// License: BSD License (Kratos default license: kratos/license.txt)
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include "includes/stream_serializer.h"
#include "testing/testing.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_isotropic_hardening_3d.h"

namespace Kratos::Testing::J2HardeningTestUtilities
{
inline Vector make_vector(std::initializer_list<double> Values)
{
    Vector result(Values.size());
    std::copy(Values.begin(), Values.end(), result.begin());
    return result;
}
inline Matrix elastic_reference(const Properties& rProperties)
{
    const double e = rProperties[YOUNG_MODULUS], nu = rProperties[POISSON_RATIO];
    const double shear = e / (2.0 * (1.0 + nu));
    const double lame = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    Matrix c = ZeroMatrix(6,6);
    for (SizeType i = 0; i < 3; ++i) {
        for (SizeType j = 0; j < 3; ++j) c(i,j) = lame + (i == j ? 2.0 * shear : 0.0);
        c(i+3,i+3) = shear;
    }
    return c;
}

inline double q_reference(const Vector& rStress)
{
    return std::sqrt(0.5 * (std::pow(rStress[0]-rStress[1],2) +
        std::pow(rStress[1]-rStress[2],2) + std::pow(rStress[2]-rStress[0],2)) +
        3.0 * (rStress[3]*rStress[3] + rStress[4]*rStress[4] + rStress[5]*rStress[5]));
}

inline std::string state_trace(ConstitutiveLaw& rLaw)
{
    StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ERROR);
    serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
    serializer.save("Law", rLaw);
    return serializer.GetStringRepresentation();
}
struct TangentSweepResult
{
    int BranchFaithfulSteps = 0, AcceptedSteps = 0, MaxConsecutive = 0;
    double MaxMatrixError = 0.0, MaxDirectionError = 0.0;
};

template<class TPoint, class THardening>
TangentSweepResult check_tangent_sweep(TPoint& rPoint, bool Plastic, THardening Hardening)
{
    const Vector strain = rPoint.strain;
    const Matrix c = elastic_reference(rPoint.properties);
    const Vector old_plastic = rPoint.Plastic();
    const double threshold = Hardening(rPoint.properties, rPoint.P());
    // Independent elastic predictor and scalar yield formula determine the branch.
    // The guard exceeds residual/roundoff resolution; no derivative at the switch.
    const double guard = std::max(1.0e-9 * threshold,
        512.0 * std::numeric_limits<double>::epsilon() * norm_inf(c));
    auto on_branch = [&](const Vector& rStrain) {
        const Vector predictor = prod(c, rStrain - old_plastic);
        const double margin = q_reference(predictor) - threshold;
        return Plastic ? margin > guard : margin < -guard;
    };
    KRATOS_EXPECT_TRUE(on_branch(strain));
    rPoint.Evaluate();
    const Vector saved_stress = rPoint.stress;
    const Matrix tangent = rPoint.tangent;
    const std::string state = state_trace(rPoint.law);
    const double tangent_norm = norm_inf(tangent);
    const double bound = 2.0e-6 * tangent_norm + 1.0e-8 * norm_inf(c);
    const Vector direction = make_vector({0.8,-0.3,0.5,0.4,-0.2,0.1});
    auto stress_at = [&](const Vector& rStrain) {
        rPoint.strain = rStrain;
        Vector stress;
        rPoint.law.CalculateValue(rPoint.values, CAUCHY_STRESS_VECTOR, stress);
        KRATOS_EXPECT_VECTOR_NEAR(rPoint.strain, rStrain, 0.0);
        return stress;
    };
    TangentSweepResult result;
    int consecutive = 0;
    for (double h : {1.0e-5,3.0e-6,1.0e-6,3.0e-7,1.0e-7,3.0e-8,1.0e-8,3.0e-9,1.0e-9}) {
        bool faithful = on_branch(Vector(strain + h * direction)) && on_branch(Vector(strain - h * direction));
        for (SizeType j = 0; j < 6; ++j) {
            Vector plus = strain, minus = strain;
            plus[j] += h; minus[j] -= h;
            faithful = faithful && on_branch(plus) && on_branch(minus);
        }
        if (!faithful) { consecutive = 0; continue; }
        ++result.BranchFaithfulSteps;
        Matrix numerical = ZeroMatrix(6,6);
        for (SizeType j = 0; j < 6; ++j) {
            Vector plus = strain, minus = strain;
            plus[j] += h; minus[j] -= h;
            const Vector difference = (stress_at(plus) - stress_at(minus)) / (2.0 * h);
            for (SizeType i = 0; i < 6; ++i) numerical(i,j) = difference[i];
        }
        const double matrix_error = norm_inf(Matrix(numerical - tangent));
        const Vector plus = stress_at(Vector(strain + h * direction));
        const Vector minus = stress_at(Vector(strain - h * direction));
        const double direction_error = norm_inf(Vector((plus - minus) / (2.0 * h) - prod(tangent, direction)));
        const bool passed = matrix_error <= bound && direction_error <= bound;
        consecutive = passed ? consecutive + 1 : 0;
        result.MaxConsecutive = std::max(result.MaxConsecutive, consecutive);
        if (passed) {
            ++result.AcceptedSteps;
            result.MaxMatrixError = std::max(result.MaxMatrixError, matrix_error / tangent_norm);
            result.MaxDirectionError = std::max(result.MaxDirectionError, direction_error / tangent_norm);
        }
        KRATOS_EXPECT_EQ(state_trace(rPoint.law), state);
        KRATOS_EXPECT_VECTOR_NEAR(rPoint.stress, saved_stress, 0.0);
        KRATOS_EXPECT_MATRIX_NEAR(rPoint.tangent, tangent, 0.0);
    }
    KRATOS_EXPECT_GE(result.MaxConsecutive, 2);
    KRATOS_EXPECT_MATRIX_NEAR(tangent, Matrix(trans(tangent)), 1.0e-8);
    rPoint.strain = strain;
    Matrix queried;
    rPoint.law.CalculateValue(rPoint.values, CONSTITUTIVE_MATRIX, queried);
    KRATOS_EXPECT_MATRIX_NEAR(queried, tangent, 0.0);
    KRATOS_EXPECT_EQ(state_trace(rPoint.law), state);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.stress, saved_stress, 0.0);
    KRATOS_EXPECT_MATRIX_NEAR(rPoint.tangent, tangent, 0.0);
    return result;
}
} // namespace Kratos::Testing::J2HardeningTestUtilities
