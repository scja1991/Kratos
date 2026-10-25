// License: BSD License (Kratos default license: kratos/license.txt)
#pragma once
#include <cmath>
#include <limits>
#include "includes/properties.h"
#include "includes/variables.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
/// Stateless pure Voce sigma_y(p), d(sigma_y)/dp; no history or integration.
class VoceHardeningLaw
{
public:
    static bool TryCalculateStressThreshold(const Properties& rProperties, double P,
        double& rThreshold, double& rDerivative)
    {
        const double initial = rProperties[YIELD_STRESS];
        const double saturation = rProperties[EXPONENTIAL_SATURATION_YIELD_STRESS];
        const double b = rProperties[HARDENING_EXPONENT];
        if (!std::isfinite(P) || P < 0.0 || !std::isfinite(initial) || initial <= 0.0 ||
            !std::isfinite(saturation) || saturation < initial || !std::isfinite(b) || b < 0.0) return false;
        const double q = saturation - initial;
        double threshold = initial, derivative = 0.0;
        if (q > 0.0 && b > 0.0) {
            const double x = b * P;
            if (std::isinf(x)) {
                threshold = saturation;
            } else {
                const double decay = std::exp(-x);
                const double fraction = -std::expm1(-x);
                // Accurate at small bp; the second form avoids a large positive sum.
                const double initial_derivative = q * b;
                if (P > 0.0 && x < std::numeric_limits<double>::min()) {
                    // bp may underflow while Q*bp is still representable. The
                    // omitted relative correction is O(bp), below double resolution.
                    threshold = std::isnormal(initial_derivative)
                        ? std::fma(initial_derivative, P, initial)
                        : initial + std::exp(std::log(q) + std::log(b) + std::log(P));
                } else {
                    threshold = fraction <= 0.5 ? std::fma(q, fraction, initial)
                                                : saturation - q * decay;
                }
                // Avoid both overflow*0 and premature underflow of exp(-bp).
                derivative = std::isfinite(initial_derivative) && std::isnormal(decay)
                    ? initial_derivative * decay : CalculateScaledDerivative(q, b, x);
            }
        }
        if (!std::isfinite(threshold) || threshold <= 0.0 ||
            !std::isfinite(derivative) || derivative < 0.0) return false;
        rThreshold = threshold;
        rDerivative = derivative;
        return true;
    }

    static void CalculateStressThreshold(const Properties& rProperties, double P,
        double& rThreshold, double& rDerivative)
    {
        KRATOS_ERROR_IF_NOT(TryCalculateStressThreshold(rProperties, P, rThreshold, rDerivative))
            << "Voce hardening is not representable at p=" << P;
    }

    static int Check(const Properties& rProperties)
    {
        for (const auto* p_variable : {&YIELD_STRESS, &EXPONENTIAL_SATURATION_YIELD_STRESS, &HARDENING_EXPONENT}) {
            KRATOS_ERROR_IF_NOT(rProperties.Has(*p_variable)) << "Missing " << p_variable->Name();
            KRATOS_ERROR_IF_NOT(std::isfinite(rProperties[*p_variable])) << "Nonfinite " << p_variable->Name();
        }
        KRATOS_ERROR_IF(rProperties[YIELD_STRESS] <= 0.0 ||
            rProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] < rProperties[YIELD_STRESS] ||
            rProperties[HARDENING_EXPONENT] < 0.0) << "Invalid Voce hardening parameters";
        if (rProperties.Has(ISOTROPIC_HARDENING_MODULUS)) {
            KRATOS_ERROR_IF(!std::isfinite(rProperties[ISOTROPIC_HARDENING_MODULUS]) ||
                rProperties[ISOTROPIC_HARDENING_MODULUS] != 0.0)
                << "Pure Voce does not support a nonzero ISOTROPIC_HARDENING_MODULUS";
        }
        double threshold = 0.0, derivative = 0.0;
        CalculateStressThreshold(rProperties, 0.0, threshold, derivative);
        return 0;
    }

private:
    static double CalculateScaledDerivative(double Q, double B, double X)
    {
        // Q,B < 2^1024, so beyond 4096 even their product times exp(-X)
        // is far below half the smallest subnormal double. No overflow clamp.
        if (X > 4096.0) return 0.0;
        int q_exponent = 0, b_exponent = 0;
        double mantissa = std::frexp(Q, &q_exponent) * std::frexp(B, &b_exponent);
        int exponent = q_exponent + b_exponent;
        const double chunk_decay = std::exp(-512.0);
        // Each product remains normal. Carry its binary exponent separately;
        // this avoids both Q*B overflow and premature exponential underflow.
        while (X > 512.0) {
            int shift = 0;
            mantissa = std::frexp(mantissa * chunk_decay, &shift);
            exponent += shift;
            X -= 512.0;
        }
        return std::scalbn(mantissa * std::exp(-X), exponent);
    }
};
} // namespace Kratos
