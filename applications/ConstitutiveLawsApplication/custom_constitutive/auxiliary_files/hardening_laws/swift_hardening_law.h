// License: BSD License (Kratos default license: kratos/license.txt)
#pragma once
#include <cmath>
#include "includes/properties.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
/// Stateless sigma_y(p), d(sigma_y)/dp; no integration or history ownership.
class SwiftHardeningLaw
{
public:
    static bool TryCalculateStressThreshold(const Properties& rProperties, double P,
        double& rThreshold, double& rDerivative)
    {
        const double base = rProperties[SWIFT_INITIAL_STRAIN] + P;
        if (!std::isfinite(P) || P < 0.0 || !std::isfinite(base) || base <= 0.0) return false;
        const double n = rProperties[SWIFT_HARDENING_EXPONENT];
        rThreshold = rProperties[SWIFT_COEFFICIENT] * std::pow(base, n);
        rDerivative = n == 0.0 ? 0.0 : n * (rThreshold / base);
        return std::isfinite(rThreshold) && rThreshold > 0.0 &&
            std::isfinite(rDerivative) && rDerivative >= 0.0;
    }

    static void CalculateStressThreshold(const Properties& rProperties, double P,
        double& rThreshold, double& rDerivative)
    {
        KRATOS_ERROR_IF_NOT(TryCalculateStressThreshold(rProperties, P, rThreshold, rDerivative))
            << "Swift hardening is not representable at p=" << P;
    }

    static int Check(const Properties& rProperties)
    {
        for (const auto* p_variable : {&SWIFT_COEFFICIENT, &SWIFT_INITIAL_STRAIN, &SWIFT_HARDENING_EXPONENT}) {
            KRATOS_ERROR_IF_NOT(rProperties.Has(*p_variable)) << "Missing " << p_variable->Name();
            KRATOS_ERROR_IF_NOT(std::isfinite(rProperties[*p_variable])) << "Nonfinite " << p_variable->Name();
        }
        KRATOS_ERROR_IF(rProperties[SWIFT_COEFFICIENT] <= 0.0 || rProperties[SWIFT_INITIAL_STRAIN] <= 0.0 ||
            rProperties[SWIFT_HARDENING_EXPONENT] < 0.0) << "Invalid Swift hardening parameters";
        double threshold = 0.0, derivative = 0.0;
        CalculateStressThreshold(rProperties, 0.0, threshold, derivative);
        return 0;
    }
};
} // namespace Kratos
