// License: BSD License (Kratos default license: kratos/license.txt)
#include <array>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include "containers/model.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"
#include "tests/cpp_tests/j2_hardening_test_utilities.h"
#include "custom_constitutive/auxiliary_files/hardening_laws/swift_hardening_law.h"

namespace Kratos::Testing
{
namespace
{
using namespace J2HardeningTestUtilities;

class InspectableSwift : public SmallStrainSwiftJ2Plasticity3D
{
public:
    using SmallStrainJ2IsotropicHardening3D::EvaluateTrial;
    double Cache() { return GetThreshold(); }
    double Dissipation() { return GetPlasticDissipation(); }
};

class ScaledElasticSwift : public SmallStrainSwiftJ2Plasticity3D
{
public:
    int ElasticCalls = 0;
protected:
    void CalculateElasticMatrix(Matrix& rMatrix, ConstitutiveLaw::Parameters& rValues) override
    {
        ++ElasticCalls;
        ElasticIsotropic3D::CalculateElasticMatrix(rMatrix, rValues);
        rMatrix *= 1.25;
    }
};

// Reconstructed writers, NOT archived historical files. ConstitutiveLaw's
// Flags/InitialState chain is unchanged between the cited commits and this base.
// ab168b6f962: standalone Swift -> ConstitutiveLaw -> plastic vector -> p.
class ReconstructedStandaloneSwiftWriter : public ConstitutiveLaw
{
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.save("PlasticStrain", make_vector({0.002,-0.001,-0.001,0.0,0.0,0.0}));
        rSerializer.save("AccumulatedPlasticStrain", 0.002);
    }
};

// b2dde3d44cc: wrapper -> historical generic base -> ConstitutiveLaw ->
// dissipation, threshold, plastic vector, p. Keep BOTH base-class trace levels.
class ReconstructedGenericSwiftBase : public ConstitutiveLaw
{
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
        rSerializer.save("PlasticDissipation", 0.0);
        rSerializer.save("Threshold", 400.0);
        rSerializer.save("PlasticStrain", make_vector({0.002,-0.001,-0.001,0.0,0.0,0.0}));
        rSerializer.save("AccumulatedPlasticStrain", 0.002);
    }
};

template<class TBase>
class ReconstructedSwiftWrapper : public TBase
{
public:
    using TBase::TBase;
private:
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBase);
    }
};

// Immediately preceding local prototype: wrapper -> family -> unchanged master
// generic base -> marker, version 1, p. This format is deliberately not migrated.
class ReconstructedPrototypeV1Family : public SmallStrainJ2IsotropicHardening3D::BaseType
{
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const override
    {
        using BaseType = SmallStrainJ2IsotropicHardening3D::BaseType;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("HardeningStateMarker", std::uint64_t{0x7ff84a3248415244ULL});
        rSerializer.save("HardeningStateVersion", 1);
        rSerializer.save("AccumulatedPlasticStrain", 0.002);
    }
};

// Deliberately no base payload: the requested error must occur at the header.
class HeaderOnlyFamilyWriter : public ConstitutiveLaw
{
public:
    HeaderOnlyFamilyWriter(std::uint64_t Marker, int Version) : mMarker(Marker), mVersion(Version) {}
private:
    std::uint64_t mMarker;
    int mVersion;
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const override
    {
        rSerializer.save("HardeningStateMarker", mMarker);
        rSerializer.save("HardeningStateVersion", mVersion);
    }
};

struct SwiftPoint
{
    Model model;
    ModelPart& part;
    Tetrahedra3D4<Node>::Pointer geometry;
    Properties properties;
    Vector strain = ZeroVector(6), stress = ZeroVector(6), n = Vector(4, 0.25);
    Matrix tangent = ZeroMatrix(6,6);
    ConstitutiveLaw::Parameters values;
    InspectableSwift law;

    SwiftPoint() : part(model.CreateModelPart("Main"))
    {
        geometry = Kratos::make_shared<Tetrahedra3D4<Node>>(
            part.CreateNewNode(1,0,0,0), part.CreateNewNode(2,1,0,0),
            part.CreateNewNode(3,0,1,0), part.CreateNewNode(4,0,0,1));
        properties[YOUNG_MODULUS] = 210000.0; // MPa throughout these tests.
        properties[POISSON_RATIO] = 0.3;
        properties[SWIFT_COEFFICIENT] = 1000.0;
        properties[SWIFT_INITIAL_STRAIN] = 0.01;
        properties[SWIFT_HARDENING_EXPONENT] = 0.2;
        part.GetProcessInfo()[STEP] = 1;
        part.GetProcessInfo()[NL_ITERATION_NUMBER] = 1;
        values.SetElementGeometry(*geometry);
        values.SetMaterialProperties(properties);
        values.SetProcessInfo(part.GetProcessInfo());
        values.SetShapeFunctionsValues(n);
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        law.InitializeMaterial(properties, *geometry, n);
    }
    double P(ConstitutiveLaw* pLaw = nullptr)
    {
        double p = 0.0;
        (pLaw ? pLaw : &law)->GetValue(ACCUMULATED_PLASTIC_STRAIN, p);
        return p;
    }
    Vector Plastic(ConstitutiveLaw* pLaw = nullptr)
    {
        Vector ep;
        (pLaw ? pLaw : &law)->GetValue(PLASTIC_STRAIN_VECTOR, ep);
        return ep;
    }
    void Evaluate() { law.CalculateMaterialResponseCauchy(values); }
    void Commit() { law.FinalizeMaterialResponseCauchy(values); }
};

Vector flow_reference(const Vector& rStress)
{
    Matrix tensor = MathUtils<double>::StressVectorToTensor(rStress);
    const double pressure = (tensor(0,0) + tensor(1,1) + tensor(2,2)) / 3.0;
    for (SizeType i=0; i<3; ++i) tensor(i,i) -= pressure;
    tensor *= 1.5 / q_reference(rStress);
    return MathUtils<double>::StrainTensorToVector(tensor, 6);
}

double yield_reference(const Properties& rProperties, double P)
{
    return rProperties[SWIFT_COEFFICIENT] * std::pow(rProperties[SWIFT_INITIAL_STRAIN]+P,
        rProperties[SWIFT_HARDENING_EXPONENT]);
}

struct Reference { Vector Stress; Vector Plastic; double P; double Lambda; };
Reference reference_step(SwiftPoint& rPoint)
{
    const Matrix c = elastic_reference(rPoint.properties);
    Vector mechanical = rPoint.strain;
    Vector initial_stress = ZeroVector(6);
    if (rPoint.law.HasInitialState()) {
        mechanical -= rPoint.law.GetInitialState().GetInitialStrainVector();
        initial_stress = rPoint.law.GetInitialState().GetInitialStressVector();
    }
    const Vector old_plastic = rPoint.Plastic();
    const double old_p = rPoint.P();
    const Vector predictor = prod(c, mechanical-old_plastic) + initial_stress;
    const double q = q_reference(predictor);
    if (q <= yield_reference(rPoint.properties, old_p)) return {predictor,old_plastic,old_p,0.0};
    const double shear = rPoint.properties[YOUNG_MODULUS] / (2.0*(1.0+rPoint.properties[POISSON_RATIO]));
    double lower = 0.0, upper = (q-yield_reference(rPoint.properties,old_p))/(3.0*shear);
    // Independent scalar bisection: no production solver, gradient, Hessian or H.
    for (int i=0; i<120; ++i) {
        const double mid = (lower+upper)/2.0;
        if (q-3.0*shear*mid-yield_reference(rPoint.properties,old_p+mid) > 0.0) lower=mid;
        else upper=mid;
    }
    const double lambda=(lower+upper)/2.0;
    const Vector increment=lambda*flow_reference(predictor);
    return {Vector(predictor-prod(c,increment)), Vector(old_plastic+increment), old_p+lambda, lambda};
}

void compare_response(SwiftPoint& rPoint, const Reference& rReference)
{
    rPoint.Evaluate();
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.stress, rReference.Stress, 1.0e-7+1.0e-9*norm_inf(rReference.Stress));
    const auto trial = rPoint.law.EvaluateTrial(rPoint.values, true);
    KRATOS_EXPECT_TRUE(trial.Candidate.has_value());
    KRATOS_EXPECT_LE(trial.StressResidual, trial.StressTolerance);
    if (rReference.Lambda > 0.0) {
        KRATOS_EXPECT_LE(std::abs(trial.YieldResidual), trial.StressTolerance);
    }
    KRATOS_EXPECT_NEAR(trial.Candidate->P, rReference.P, 1.0e-12+1.0e-9*rReference.P);
    KRATOS_EXPECT_VECTOR_NEAR(trial.Candidate->PlasticStrain, rReference.Plastic, 1.0e-12);
}

void check_committed_yield(SwiftPoint& rPoint, InspectableSwift& rLaw)
{
    const std::string before = state_trace(rLaw);
    const Vector stress = rPoint.stress;
    const Matrix tangent = rPoint.tangent;
    double threshold = -1.0;
    rLaw.CalculateValue(rPoint.values, YIELD_STRESS, threshold);
    const double reference = yield_reference(rPoint.properties, rPoint.P(&rLaw));
    KRATOS_EXPECT_NEAR(threshold, reference, 1.0e-12 * reference);
    KRATOS_EXPECT_EQ(state_trace(rLaw), before);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.stress, stress, 0.0);
    KRATOS_EXPECT_MATRIX_NEAR(rPoint.tangent, tangent, 0.0);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointScalarFunction, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    for (double n : {0.0,0.2,1.0}) {
        point.properties[SWIFT_HARDENING_EXPONENT]=n;
        KRATOS_EXPECT_EQ(SwiftHardeningLaw::Check(point.properties),0);
        for (double p : {0.0,1.0e-9,0.003,0.1,2.0}) {
            double y=0.0,h=0.0;
            SwiftHardeningLaw::CalculateStressThreshold(point.properties,p,y,h);
            const double base=0.01+p;
            KRATOS_EXPECT_NEAR(y,1000.0*std::pow(base,n),1.0e-12*y);
            const double reference_h=n==0.0 ? 0.0 : 1000.0*n*std::pow(base,n-1.0);
            KRATOS_EXPECT_NEAR(h,reference_h,1.0e-12*std::max(1.0,reference_h));
            if (p>0.001) {
                const double step=base*1.0e-5;
                const double numerical=(yield_reference(point.properties,p+step)-yield_reference(point.properties,p-step))/(2*step);
                KRATOS_EXPECT_NEAR(h,numerical,1.0e-7*std::max(1.0,h));
            }
        }
    }
    point.properties[SWIFT_INITIAL_STRAIN]=0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SwiftHardeningLaw::Check(point.properties),"Invalid Swift");
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointEngineeringNorm, KratosConstitutiveLawsFastSuite)
{
    for (const Vector& increment : {make_vector({0.02,-0.01,-0.01,0.0,0.0,0.0}),
        make_vector({0.0,0.0,0.0,0.03,0.0,0.0}),make_vector({0.02,-0.003,0.006,0.01,-0.04,0.05})}) {
        const Matrix tensor=MathUtils<double>::StrainVectorToTensor(increment);
        double contraction=0.0;
        for (SizeType i=0;i<3;++i) for (SizeType j=0;j<3;++j) contraction+=tensor(i,j)*tensor(i,j);
        KRATOS_EXPECT_NEAR(AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(increment),
            std::sqrt(2.0*contraction/3.0),1.0e-15);
    }
    const Vector shear = make_vector({0.0,0.0,0.0,0.03,0.0,0.0});
    KRATOS_EXPECT_NEAR(AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(shear),0.03/std::sqrt(3.0),1.0e-15);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(ZeroVector(3)),"complete 3D");
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointPotentialLinearization, KratosConstitutiveLawsFastSuite)
{
    for (const Vector& stress : {make_vector({500.0,0.0,0.0,0.0,0.0,0.0}),make_vector({0.0,0.0,0.0,200.0,0.0,0.0}),
        make_vector({421.0,-93.0,82.0,100.0,50.0,-70.0})}) {
        array_1d<double,6> input=stress,g;
        Matrix dg;
        KRATOS_EXPECT_TRUE(VonMisesPlasticPotential<6>::CalculatePlasticPotentialDerivativeAndJacobian(input,q_reference(stress),g,dg));
        for (double step : {0.003,0.001}) {
            for (SizeType j=0;j<6;++j) {
                Vector plus=stress,minus=stress;
                plus[j]+=step;minus[j]-=step;
                KRATOS_EXPECT_NEAR(g[j],(q_reference(plus)-q_reference(minus))/(2.0*step),2.0e-9);
                const Vector derivative=(flow_reference(plus)-flow_reference(minus))/(2.0*step);
                for (SizeType i=0;i<6;++i) {
                    KRATOS_EXPECT_NEAR(dg(i,j),derivative[i],1.0e-10);
                    KRATOS_EXPECT_NEAR(dg(i,j),dg(j,i),1.0e-15);
                }
            }
        }
        KRATOS_EXPECT_LT(norm_2(Vector(prod(dg,input))),1.0e-14);
    }
    array_1d<double,6> hydro=ScalarVector(6,0.0),g;
    Matrix dg;
    KRATOS_EXPECT_FALSE(VonMisesPlasticPotential<6>::CalculatePlasticPotentialDerivativeAndJacobian(hydro,0.0,g,dg));
    SwiftPoint point;
    point.strain[0]=point.strain[1]=point.strain[2]=0.001;
    point.Evaluate();point.Commit();
    KRATOS_EXPECT_EQ(point.P(),0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointNonlinearPaths, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    int maximum_iterations=0;
    double accumulated=0.0;
    for (const Vector& strain : {make_vector({0.0001,0.0,0.0,0.0,0.0,0.0}),make_vector({0.006,0.0,0.0,0.0,0.0,0.0}),
        make_vector({0.005,0.0,0.0,0.0,0.0,0.0}),make_vector({-0.008,0.0,0.0,0.0,0.0,0.0}),
        make_vector({0.0,0.0,0.0,0.02,0.0,0.0}),make_vector({0.006,-0.004,-0.002,-0.012,0.009,0.004})}) {
        point.strain=strain;
        const auto reference=reference_step(point);
        const Vector old_plastic=point.Plastic();
        const double old_p=point.P();
        const std::string state=state_trace(point.law);
        compare_response(point,reference);
        const auto result=point.law.EvaluateTrial(point.values,true);
        maximum_iterations=std::max(maximum_iterations,result.Iterations);
        KRATOS_EXPECT_EQ(state,state_trace(point.law));
        point.Commit();
        KRATOS_EXPECT_NEAR(point.P(),reference.P,1.0e-12);
        const Matrix increment=MathUtils<double>::StrainVectorToTensor(Vector(point.Plastic()-old_plastic));
        double contraction=0.0;
        for(SizeType i=0;i<3;++i) for(SizeType j=0;j<3;++j) contraction+=increment(i,j)*increment(i,j);
        accumulated+=std::sqrt(2.0*contraction/3.0);
        KRATOS_EXPECT_NEAR(point.P(),accumulated,1.0e-12);
        KRATOS_EXPECT_GE(point.P(),old_p);
        const double accepted_p=point.P();
        point.Commit();
        KRATOS_EXPECT_EQ(point.P(),accepted_p);
    }
    KRATOS_EXPECT_GT(maximum_iterations,1);
    KRATOS_EXPECT_GT(point.P(),AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(point.Plastic())+0.001);
    const double projection=inner_prod(point.stress,point.Plastic())/q_reference(point.stress);
    KRATOS_EXPECT_GT(std::abs(point.P()-projection),0.001);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointUniaxialStress, KratosConstitutiveLawsFastSuite)
{
    for (double sign : {-1.0, 1.0}) {
        for (double axial : {0.0001, 0.004, 0.02}) {
            SwiftPoint point;
            const double e = point.properties[YOUNG_MODULUS];
            const double nu = point.properties[POISSON_RATIO];
            double p = 0.0;
            if (e * axial > yield_reference(point.properties, 0.0)) {
                double lower = 0.0, upper = axial;
                // Independent uniaxial equation: epsilon_x = sigma/E + p.
                for (int i = 0; i < 120; ++i) {
                    const double mid = 0.5 * (lower + upper);
                    if (e * (axial - mid) > yield_reference(point.properties, mid)) lower = mid;
                    else upper = mid;
                }
                p = 0.5 * (lower + upper);
            }
            const double sigma = e * (axial - p);
            point.strain = make_vector({sign * axial, sign * (-nu * sigma / e - 0.5 * p),
                sign * (-nu * sigma / e - 0.5 * p), 0.0, 0.0, 0.0});
            point.Evaluate();
            KRATOS_EXPECT_VECTOR_NEAR(point.stress, make_vector({sign*sigma,0.0,0.0,0.0,0.0,0.0}), 1.0e-7);
            point.Commit();
            KRATOS_EXPECT_NEAR(point.P(), p, 1.0e-12);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointNewtonCorrectionsAreNotHistory, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.properties[SWIFT_COEFFICIENT] = 1.0e6;
    point.properties[SWIFT_HARDENING_EXPONENT] = 2.0;
    point.law.ResetMaterial(point.properties, *point.geometry, point.n);
    point.strain[0] = 0.02;
    const auto reference = reference_step(point);
    const double shear = point.properties[YOUNG_MODULUS] / (2.0 * (1.0 + point.properties[POISSON_RATIO]));
    const Vector predictor = prod(elastic_reference(point.properties), point.strain);
    const double initial_h = 2.0 * 1.0e6 * 0.01;
    const double first_lambda = (q_reference(predictor) - 100.0) / (3.0 * shear + initial_h);
    KRATOS_EXPECT_GT(first_lambda, reference.Lambda + 1.0e-4); // A later Newton correction must reverse.
    compare_response(point, reference);
    const auto result = point.law.EvaluateTrial(point.values, false);
    KRATOS_EXPECT_GT(result.Iterations, 1);
    point.Commit();
    KRATOS_EXPECT_NEAR(point.P(), reference.P, 1.0e-12);
    KRATOS_EXPECT_LT(point.P(), first_lambda);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointLinearLimits, KratosConstitutiveLawsFastSuite)
{
    for (double exponent : {0.0, 1.0}) {
        SwiftPoint point;
        point.properties[SWIFT_HARDENING_EXPONENT] = exponent;
        point.law.ResetMaterial(point.properties, *point.geometry, point.n);
        for (double axial : {0.002, 0.02}) {
            point.strain[0] = axial;
            compare_response(point, reference_step(point));
            point.Commit();
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointDiscreteTangent, KratosConstitutiveLawsFastSuite)
{
    enum class Path { Elastic, Tension, Compression, Shear, Unload, Reload, NonProportional, BelowYield, AboveYield };
    const std::array<std::pair<const char*,Path>,9> cases{{
        {"elastic",Path::Elastic}, {"tension",Path::Tension}, {"compression",Path::Compression},
        {"shear",Path::Shear}, {"unload",Path::Unload}, {"reload",Path::Reload},
        {"nonproportional",Path::NonProportional}, {"below_yield",Path::BelowYield}, {"above_yield",Path::AboveYield}}};
    double maximum_matrix = 0.0, maximum_direction = 0.0;
    int accepted = 0, faithful = 0, minimum_adjacent = 9;
    for (double exponent : {0.0,0.2,1.0,2.0}) for (const auto& r_case : cases) {
        std::ostringstream label;
        label << "n" << exponent << "_" << r_case.first;
        SCOPED_TRACE(label.str());
        SwiftPoint point;
        point.properties[SWIFT_HARDENING_EXPONENT] = exponent;
        point.properties[SWIFT_COEFFICIENT] = 400.0 / std::pow(0.01, exponent);
        point.law.ResetMaterial(point.properties, *point.geometry, point.n);
        const double e = point.properties[YOUNG_MODULUS], nu = point.properties[POISSON_RATIO];
        const double two_g = e / (1.0 + nu);
        bool plastic = true;
        switch (r_case.second) {
            case Path::Elastic:
                point.strain[0] = 0.0001; plastic = false;
                break;
            case Path::Tension: case Path::Compression: {
                // Independent uniaxial-stress construction, E*(eps-p)=sigma_y(p).
                double lower = 0.0, upper = 0.008;
                for (int i = 0; i < 120; ++i) {
                    const double mid = (lower + upper) / 2.0;
                    if (e * (0.008 - mid) > yield_reference(point.properties, mid)) lower = mid;
                    else upper = mid;
                }
                const double p = (lower + upper) / 2.0, stress = e * (0.008 - p);
                const double sign = r_case.second == Path::Tension ? 1.0 : -1.0;
                point.strain = make_vector({sign*0.008,sign*(-nu*stress/e-p/2.0),sign*(-nu*stress/e-p/2.0),0.0,0.0,0.0});
                break;
            }
            case Path::Shear:
                point.strain[3] = 0.01;
                break;
            case Path::Unload: case Path::Reload: {
                point.strain[0] = 0.008; point.Commit();
                point.strain = point.Plastic();
                point.strain[0] += 0.5 * yield_reference(point.properties, point.P()) / two_g;
                if (r_case.second == Path::Unload) plastic = false;
                else { point.Commit(); point.strain = make_vector({0.012,0.0,0.0,0.0,0.0,0.0}); }
                break;
            }
            case Path::NonProportional:
                point.strain = make_vector({0.007,-0.001,-0.001,0.003,0.0,0.0}); point.Commit();
                point.strain = make_vector({0.012,-0.004,0.001,-0.006,0.004,0.002});
                break;
            case Path::BelowYield: case Path::AboveYield:
                plastic = r_case.second == Path::AboveYield;
                point.strain[0] = 400.0 / two_g * (1.0 + (plastic ? 1.0 : -1.0) * 1.0e-4);
                break;
        }
        const auto result = check_tangent_sweep(point, plastic, yield_reference);
        maximum_matrix = std::max(maximum_matrix, result.MaxMatrixError);
        maximum_direction = std::max(maximum_direction, result.MaxDirectionError);
        accepted += result.AcceptedSteps; faithful += result.BranchFaithfulSteps;
        minimum_adjacent = std::min(minimum_adjacent, result.MaxConsecutive);
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(17) << "branch_faithful=" << result.BranchFaithfulSteps
            << ";accepted=" << result.AcceptedSteps << ";adjacent=" << result.MaxConsecutive
            << ";max_matrix=" << result.MaxMatrixError << ";max_direction=" << result.MaxDirectionError;
        ::testing::Test::RecordProperty(label.str(), diagnostic.str());
    }
    ::testing::Test::RecordProperty("tangent_cases", 36);
    ::testing::Test::RecordProperty("branch_faithful_steps", faithful);
    ::testing::Test::RecordProperty("accepted_steps", accepted);
    ::testing::Test::RecordProperty("minimum_adjacent_steps", minimum_adjacent);
    for (const auto& metric : {std::make_pair("max_relative_matrix_error",maximum_matrix),
                              std::make_pair("max_relative_direction_error",maximum_direction)}) {
        std::ostringstream value; value << std::setprecision(17) << metric.second;
        ::testing::Test::RecordProperty(metric.first, value.str());
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointUsesVirtualElasticity, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    ScaledElasticSwift law;
    law.InitializeMaterial(point.properties, *point.geometry, point.n);
    point.strain[0] = 0.02;
    point.properties[YOUNG_MODULUS] *= 1.25;
    const auto reference = reference_step(point);
    point.properties[YOUNG_MODULUS] /= 1.25;
    law.CalculateMaterialResponseCauchy(point.values);
    KRATOS_EXPECT_EQ(law.ElasticCalls, 1);
    KRATOS_EXPECT_VECTOR_NEAR(point.stress, reference.Stress, 1.0e-7);
    law.FinalizeMaterialResponseCauchy(point.values);
    KRATOS_EXPECT_EQ(law.ElasticCalls, 2);
    KRATOS_EXPECT_NEAR(point.P(&law), reference.P, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointResultStatus, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    using Integrator = SmallStrainJ2IsotropicHardening3D::Integrator;
    Integrator::TrialInput input;
    input.ElasticMatrix=elastic_reference(point.properties);
    input.EffectiveStrain=make_vector({0.02,-0.003,0.0,0.02,0.0,0.0});
    noalias(input.PredictorStress)=prod(input.ElasticMatrix,input.EffectiveStrain);
    input.NeedTangent=true;
    const Integrator::EvolutionContext context{point.properties,&SwiftHardeningLaw::TryCalculateStressThreshold,
        [](const array_1d<double,6>& rStress,double Q,Integrator::FlowLinearization& rFlow) {
            const bool valid=VonMisesPlasticPotential<6>::CalculatePlasticPotentialDerivativeAndJacobian(rStress,Q,rFlow.G,rFlow.DG);
            rFlow.F=rFlow.G;
            return valid;
        },
        [](double P,double Lambda,const Integrator::FlowLinearization&,Integrator::HistoryResponse& rHistory) {
            rHistory.P=P+Lambda;rHistory.DPLambda=1.0;return true;
        }};
    input.MaxIterations=1;
    const auto failed=Integrator::IntegrateStressVectorWithHistory(input,context,point.values);
    KRATOS_EXPECT_EQ(failed.Status,Integrator::IntegrationStatus::NoConvergence);
    KRATOS_EXPECT_EQ(failed.Failure,Integrator::IntegrationFailure::IterationLimit);
    KRATOS_EXPECT_FALSE(failed.Candidate.has_value());
    input.MaxIterations=100;
    const auto success=Integrator::IntegrateStressVectorWithHistory(input,context,point.values);
    KRATOS_EXPECT_EQ(success.Status,Integrator::IntegrationStatus::Converged);
    KRATOS_EXPECT_GT(success.Iterations,1);
    KRATOS_EXPECT_TRUE(success.Candidate.has_value());
    const Integrator::EvolutionContext hydro_context{point.properties,context.Hardening,
        [](const array_1d<double,6>&,double,Integrator::FlowLinearization&) { return false; },context.History};
    input.EffectiveStrain=make_vector({0.001,0.001,0.001,0.0,0.0,0.0});
    noalias(input.PredictorStress)=prod(input.ElasticMatrix,input.EffectiveStrain);
    const auto hydro=Integrator::IntegrateStressVectorWithHistory(input,hydro_context,point.values);
    KRATOS_EXPECT_EQ(hydro.Status,Integrator::IntegrationStatus::Elastic);
    KRATOS_EXPECT_TRUE(hydro.Candidate.has_value());
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointFailedFinalizeAtomicity, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.strain[0]=0.007;point.Commit();
    const std::string state=state_trace(point.law);
    point.strain=make_vector({0.03,-0.005,0.002,0.025,0.003,0.0});
    point.properties[MAX_NUMBER_NL_CL_ITERATIONS]=1;
    point.stress=ScalarVector(6,-123.0);point.tangent=ScalarMatrix(6,6,-456.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"iteration limit");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Commit(),"iteration limit");
    KRATOS_EXPECT_EQ(state,state_trace(point.law));
    KRATOS_EXPECT_VECTOR_NEAR(point.stress,Vector(ScalarVector(6,-123.0)),0.0);
    KRATOS_EXPECT_MATRIX_NEAR(point.tangent,Matrix(ScalarMatrix(6,6,-456.0)),0.0);
    point.properties[MAX_NUMBER_NL_CL_ITERATIONS]=100;
    compare_response(point,reference_step(point));
    const double old_p=point.P();point.Commit();KRATOS_EXPECT_GT(point.P(),old_p);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointFlagsAndInitialState, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    Vector initial_strain = make_vector({0.001,0.0,0.0,0.0002,0.0,0.0});
    Vector initial_stress = make_vector({40.0,-10.0,0.0,15.0,0.0,0.0});
    point.law.SetInitialState(Kratos::make_intrusive<InitialState>(initial_strain,initial_stress));
    point.strain=make_vector({0.008,-0.001,0.0,0.004,0.0,0.0});
    const Vector input=point.strain;
    const auto reference=reference_step(point);
    const std::string state=state_trace(point.law);
    point.Evaluate();
    const Matrix full_tangent = point.tangent;
    for (bool attach_unrequested : {false,true})
    for(bool stress:{false,true}) for(bool tangent:{false,true}) {
        ConstitutiveLaw::Parameters values(*point.geometry,point.properties,point.part.GetProcessInfo());
        values.SetStrainVector(point.strain);
        values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS,stress);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,tangent);
        Vector output=ScalarVector(6,-123.0);
        Matrix matrix=ScalarMatrix(6,6,-456.0);
        if(stress || attach_unrequested) values.SetStressVector(output);
        if(tangent || attach_unrequested) values.SetConstitutiveMatrix(matrix);
        point.law.CalculateMaterialResponseCauchy(values);
        if (stress) {
            KRATOS_EXPECT_VECTOR_NEAR(output,reference.Stress,1.0e-6);
        } else {
            KRATOS_EXPECT_VECTOR_NEAR(output,Vector(ScalarVector(6,-123.0)),0.0);
        }
        if (tangent) {
            KRATOS_EXPECT_MATRIX_NEAR(matrix, full_tangent, 0.0);
        } else {
            KRATOS_EXPECT_MATRIX_NEAR(matrix,Matrix(ScalarMatrix(6,6,-456.0)),0.0);
        }
        KRATOS_EXPECT_VECTOR_NEAR(point.strain,input,0.0);
        KRATOS_EXPECT_EQ(state,state_trace(point.law));
    }
    compare_response(point,reference);
    point.values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS,false);
    point.values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,false);
    const Vector stress=point.stress;const Matrix tangent=point.tangent;
    ConstitutiveLaw::Parameters finalize_values(*point.geometry, point.properties, point.part.GetProcessInfo());
    finalize_values.SetStrainVector(point.strain);
    finalize_values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    finalize_values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    finalize_values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    point.law.FinalizeMaterialResponseCauchy(finalize_values); // No output buffers, even with flags set.
    KRATOS_EXPECT_GT(point.P(),0.0); // First STEP=1/NL_ITERATION=1 must really integrate.
    KRATOS_EXPECT_VECTOR_NEAR(point.stress,stress,0.0);
    KRATOS_EXPECT_MATRIX_NEAR(point.tangent,tangent,0.0);
    KRATOS_EXPECT_VECTOR_NEAR(point.strain,input,0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointStateCopyRestartAndReset, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.strain=make_vector({0.007,-0.001,0.0,0.004,0.002,0.0});point.Commit();
    InspectableSwift copy(point.law),assigned;assigned=point.law;
    auto clone=point.law.Clone();
    KRATOS_EXPECT_EQ(clone->Info(),"SmallStrainSwiftJ2Plasticity3DLaw");
    for(auto* p_law:{static_cast<ConstitutiveLaw*>(&copy),static_cast<ConstitutiveLaw*>(&assigned),clone.get()}) {
        KRATOS_EXPECT_EQ(point.P(p_law),point.P());
        KRATOS_EXPECT_VECTOR_NEAR(point.Plastic(p_law),point.Plastic(),0.0);
    }
    InspectableSwift imported;
    imported.SetPlasticState(point.Plastic(),point.P());
    KRATOS_EXPECT_EQ(imported.Cache(),0.0);
    KRATOS_EXPECT_FALSE(imported.Has(INTERNAL_VARIABLES));
    KRATOS_EXPECT_FALSE(imported.Has(PLASTIC_DISSIPATION));
    KRATOS_EXPECT_TRUE(imported.Has(PLASTIC_STRAIN_TENSOR));
    Vector unsupported;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(imported.GetValue(INTERNAL_VARIABLES,unsupported),"Unsupported");
    Matrix tensor;
    imported.CalculateValue(point.values,PLASTIC_STRAIN_TENSOR,tensor);
    KRATOS_EXPECT_MATRIX_NEAR(tensor,MathUtils<double>::StrainVectorToTensor(point.Plastic()),0.0);
    for(auto trace:{Serializer::SERIALIZER_NO_TRACE,Serializer::SERIALIZER_TRACE_ERROR}) {
        StreamSerializer serializer(trace);
        serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
        ConstitutiveLaw::Pointer original=point.law.Clone(),restored;
        serializer.save("Law",original);
        serializer.save("Sentinel",42.125);
        serializer.load("Law",restored);
        double sentinel=0.0;serializer.load("Sentinel",sentinel);
        KRATOS_EXPECT_EQ(sentinel,42.125);
        for(const Vector& strain:{make_vector({0.004,0.0,0.0,0.002,0.0,0.0}),make_vector({0.012,-0.002,0.0,-0.007,0.005,0.0})}) {
            point.strain=strain;
            original->CalculateMaterialResponseCauchy(point.values);
            const Vector stress=point.stress;const Matrix tangent=point.tangent;
            restored->CalculateMaterialResponseCauchy(point.values);
            KRATOS_EXPECT_VECTOR_NEAR(point.stress,stress,1.0e-7);
            KRATOS_EXPECT_MATRIX_NEAR(point.tangent,tangent,1.0e-5);
            original->FinalizeMaterialResponseCauchy(point.values);restored->FinalizeMaterialResponseCauchy(point.values);
            KRATOS_EXPECT_NEAR(point.P(original.get()),point.P(restored.get()),1.0e-14);
        }
    }
    point.strain=make_vector({0.012,-0.002,0.0,-0.007,0.005,0.0});point.Evaluate();
    const Vector expected=point.stress;
    imported.CalculateMaterialResponseCauchy(point.values);
    KRATOS_EXPECT_VECTOR_NEAR(point.stress,expected,1.0e-7);
    const double original_p = point.P();
    for (auto* p_law : {static_cast<ConstitutiveLaw*>(&copy), static_cast<ConstitutiveLaw*>(&assigned), clone.get()}) {
        p_law->CalculateMaterialResponseCauchy(point.values);
        KRATOS_EXPECT_VECTOR_NEAR(point.stress, expected, 1.0e-7);
        p_law->FinalizeMaterialResponseCauchy(point.values);
        KRATOS_EXPECT_GT(point.P(p_law), original_p);
        KRATOS_EXPECT_EQ(point.P(), original_p);
    }
    point.law.ResetMaterial(point.properties,*point.geometry,point.n);
    KRATOS_EXPECT_EQ(point.P(),0.0);KRATOS_EXPECT_EQ(norm_2(point.Plastic()),0.0);
    KRATOS_EXPECT_NEAR(point.law.Cache(),yield_reference(point.properties,0.0),1.0e-12);
    KRATOS_EXPECT_EQ(point.law.Dissipation(),0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointYieldAfterImportAndRestart, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.strain = make_vector({0.007,-0.001,0.0,0.004,0.002,0.0}); point.Commit();
    InspectableSwift imported;
    imported.SetPlasticState(point.Plastic(), point.P());
    KRATOS_EXPECT_EQ(imported.Cache(), 0.0);
    check_committed_yield(point, imported);
    KRATOS_EXPECT_EQ(imported.Cache(), 0.0);
    for (auto trace : {Serializer::SERIALIZER_NO_TRACE,Serializer::SERIALIZER_TRACE_ERROR}) {
        StreamSerializer serializer(trace);
        serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
        serializer.save("Law", point.law);
        InspectableSwift restored;
        serializer.load("Law", restored);
        KRATOS_EXPECT_EQ(restored.Cache(), 0.0);
        check_committed_yield(point, restored);
        KRATOS_EXPECT_EQ(restored.Cache(), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointRejectsHistoricalArchives, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.strain[0] = 0.008; point.Commit();
    point.law.SetInitialState(Kratos::make_intrusive<InitialState>(Vector(ZeroVector(6)),Vector(ScalarVector(6,10.0))));
    const std::string state = state_trace(point.law);
    for (auto trace : {Serializer::SERIALIZER_NO_TRACE,Serializer::SERIALIZER_TRACE_ERROR}) {
        auto reject = [&](auto& rWriter) {
            StreamSerializer serializer(trace);
            serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
            serializer.save("Law", rWriter);
            // Binary rejects the marker; Trace rejects the first incompatible tag.
            const std::string message = trace == Serializer::SERIALIZER_NO_TRACE
                ? "Unsupported J2 hardening archive format" : "trace tag is not the expected one";
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(serializer.load("Law", point.law), message);
            KRATOS_EXPECT_EQ(state_trace(point.law), state);
        };
        ReconstructedStandaloneSwiftWriter standalone;
        ReconstructedSwiftWrapper<ReconstructedGenericSwiftBase> generic;
        ReconstructedSwiftWrapper<ReconstructedPrototypeV1Family> prototype_v1;
        reject(standalone); reject(generic); reject(prototype_v1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointRejectsArchiveHeaderBeforeBase, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    point.strain[0] = 0.008; point.Commit();
    point.law.SetInitialState(Kratos::make_intrusive<InitialState>(Vector(ZeroVector(6)),Vector(ScalarVector(6,10.0))));
    const std::string state = state_trace(point.law);
    constexpr std::uint64_t marker = 0x7ff84a3248415244ULL;
    // Lock the new binary family prefix independently of the production reader.
    StreamSerializer prefix;
    SmallStrainJ2IsotropicHardening3D family;
    prefix.save("Law", family);
    std::uint64_t actual_marker = 0;
    int actual_version = 0;
    prefix.load("HardeningStateMarker", actual_marker);
    prefix.load("HardeningStateVersion", actual_version);
    KRATOS_EXPECT_EQ(actual_marker, marker);
    KRATOS_EXPECT_EQ(actual_version, 2);
    for (auto trace : {Serializer::SERIALIZER_NO_TRACE,Serializer::SERIALIZER_TRACE_ERROR}) {
        for (const auto& header : {std::make_pair(marker ^ std::uint64_t{1},2),std::make_pair(marker,1),std::make_pair(marker,99)}) {
            ReconstructedSwiftWrapper<HeaderOnlyFamilyWriter> writer(header.first, header.second);
            StreamSerializer serializer(trace);
            serializer.save("Law", writer);
            const std::string message = header.first != marker
                ? "Unsupported J2 hardening archive format" : "Unsupported J2 hardening archive version";
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(serializer.load("Law", point.law), message);
            KRATOS_EXPECT_EQ(state_trace(point.law), state);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftEndpointSelectionAndScope, KratosConstitutiveLawsFastSuite)
{
    SwiftPoint point;
    const auto props_before=point.properties;
    KRATOS_EXPECT_EQ(point.law.Check(point.properties,*point.geometry,point.part.GetProcessInfo()),0);
    KRATOS_EXPECT_FALSE(point.properties.Has(HARDENING_CURVE));
    KRATOS_EXPECT_FALSE(point.properties.Has(FRACTURE_ENERGY));
    KRATOS_EXPECT_FALSE(point.properties.Has(YIELD_STRESS));
    SmallStrainJ2IsotropicHardening3D family;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(family.Check(point.properties,*point.geometry,point.part.GetProcessInfo()),"requires HARDENING_CURVE");
    point.properties[HARDENING_CURVE]=7;
    family.InitializeMaterial(point.properties,*point.geometry,point.n);
    point.strain[0]=0.01;point.Evaluate();const Vector expected=point.stress;
    family.CalculateMaterialResponsePK2(point.values);KRATOS_EXPECT_VECTOR_NEAR(point.stress,expected,0.0);
    family.CalculateMaterialResponseKirchhoff(point.values);KRATOS_EXPECT_VECTOR_NEAR(point.stress,expected,0.0);
    family.CalculateMaterialResponsePK1(point.values);KRATOS_EXPECT_VECTOR_NEAR(point.stress,expected,0.0);
    point.properties[HARDENING_CURVE]=8;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"conflicting HARDENING_CURVE");
    point.properties[HARDENING_CURVE]=9; // 8 is now the supported Voce selection.
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(family.Check(point.properties,*point.geometry,point.part.GetProcessInfo()),"Unsupported HARDENING_CURVE");
    point.properties=props_before;
    point.values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,false);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"USE_ELEMENT_PROVIDED_STRAIN");
    point.values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    point.values.GetOptions().Set(ConstitutiveLaw::U_P_LAW);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"U_P_LAW");
    point.values.GetOptions().Set(ConstitutiveLaw::U_P_LAW,false);
    point.properties[TANGENT_OPERATOR_ESTIMATION]=2;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"analytic discrete tangent");
    point.properties=props_before;
    const std::string state=state_trace(point.law);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.law.SetPlasticState(ZeroVector(3),0.0),"Invalid J2 plastic state");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.law.SetPlasticState(ZeroVector(6),-1.0),"Invalid J2 plastic state");
    KRATOS_EXPECT_EQ(state,state_trace(point.law));
    point.strain.resize(3,false);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"six engineering");
    point.strain=ZeroVector(6);point.strain[0]=std::numeric_limits<double>::quiet_NaN();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"Nonfinite strain");
    point.strain=ZeroVector(6);
    Triangle2D3<Node> triangle(point.geometry->pGetPoint(0),point.geometry->pGetPoint(1),point.geometry->pGetPoint(2));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.law.Check(point.properties,triangle,point.part.GetProcessInfo()),"3D solid");
    point.properties[THERMAL_EXPANSION_COEFFICIENT]=1.0e-5;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"Thermal J2");
    point.properties=props_before;
    Matrix initial_f=IdentityMatrix(3);initial_f(0,1)=0.1;
    point.law.SetInitialState(Kratos::make_intrusive<InitialState>(initial_f));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(),"Nontrivial initial F");
}
} // namespace Kratos::Testing
