// License: BSD License (Kratos default license: kratos/license.txt)
#include "small_strain_j2_isotropic_hardening_3d.h"
#include <cstdint>
#include "custom_constitutive/auxiliary_files/hardening_laws/swift_hardening_law.h"

namespace Kratos
{
namespace {
// Version 2 validates this prefix before interpreting any base-state payload.
constexpr std::uint64_t J2HardeningArchiveMarker = 0x7ff84a3248415244ULL;
constexpr int J2HardeningArchiveVersion = 2;
}

ConstitutiveLaw::Pointer SmallStrainJ2IsotropicHardening3D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2IsotropicHardening3D>(*this);
}

SmallStrainJ2IsotropicHardening3D::Integrator::HardeningCurveType
SmallStrainJ2IsotropicHardening3D::ResolveHardening(const Properties& rProperties) const
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(HARDENING_CURVE)) << "J2 hardening family requires HARDENING_CURVE";
    KRATOS_ERROR_IF(rProperties[HARDENING_CURVE] != static_cast<int>(Integrator::HardeningCurveType::Swift))
        << "Unsupported HARDENING_CURVE in J2 hardening family: " << rProperties[HARDENING_CURVE];
    return Integrator::HardeningCurveType::Swift;
}

SmallStrainJ2IsotropicHardening3D::Integrator::EvolutionContext
SmallStrainJ2IsotropicHardening3D::MakeEvolutionContext(const Properties& rProperties) const
{
    switch (ResolveHardening(rProperties)) {
        case Integrator::HardeningCurveType::Swift:
            return {rProperties, &SwiftHardeningLaw::TryCalculateStressThreshold, &EvaluateFlow, &EvaluateHistory};
        default:
            KRATOS_ERROR << "Unsupported hardening selection";
    }
}

int SmallStrainJ2IsotropicHardening3D::Check(const Properties& rProperties,
    const GeometryType& rGeometry, const ProcessInfo& rInfo) const
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(YOUNG_MODULUS) && rProperties.Has(POISSON_RATIO))
        << "J2 hardening requires YOUNG_MODULUS and POISSON_RATIO";
    KRATOS_ERROR_IF(!std::isfinite(rProperties[YOUNG_MODULUS]) || rProperties[YOUNG_MODULUS] <= 0.0 ||
        !std::isfinite(rProperties[POISSON_RATIO]) || !std::isfinite(rProperties[DENSITY])) << "Invalid elasticity/density";
    ElasticIsotropic3D::Check(rProperties, rGeometry, rInfo);
    KRATOS_ERROR_IF(rGeometry.WorkingSpaceDimension() != 3 || rGeometry.LocalSpaceDimension() != 3)
        << "J2 hardening requires a 3D solid geometry";
    switch (ResolveHardening(rProperties)) {
        case Integrator::HardeningCurveType::Swift:
            SwiftHardeningLaw::Check(rProperties);
            break;
        default: KRATOS_ERROR << "Unsupported hardening selection";
    }
    for (const auto* p_variable : {&YOUNG_MODULUS, &POISSON_RATIO, &SWIFT_COEFFICIENT,
        &SWIFT_INITIAL_STRAIN, &SWIFT_HARDENING_EXPONENT}) {
        KRATOS_ERROR_IF(rProperties.HasAccessor(*p_variable) || rProperties.HasTable(TEMPERATURE, *p_variable))
            << "J2 hardening supports constant isothermal properties only: " << p_variable->Name();
    }
    KRATOS_ERROR_IF(rProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "Thermal J2 hardening is not supported";
    KRATOS_ERROR_IF(rProperties.Has(TANGENT_OPERATOR_ESTIMATION) && rProperties[TANGENT_OPERATOR_ESTIMATION] != 0)
        << "J2 hardening supports the analytic discrete tangent only";
    KRATOS_ERROR_IF(rProperties.Has(MAX_NUMBER_NL_CL_ITERATIONS) && rProperties[MAX_NUMBER_NL_CL_ITERATIONS] <= 0)
        << "MAX_NUMBER_NL_CL_ITERATIONS must be positive";
    double initial_yield = 0.0, hardening = 0.0;
    const auto context = MakeEvolutionContext(rProperties);
    KRATOS_ERROR_IF_NOT(context.Hardening(rProperties, 0.0, initial_yield, hardening))
        << "Initial hardening is not representable";
    for (const auto* p_variable : {&YIELD_STRESS, &YIELD_STRESS_TENSION, &YIELD_STRESS_COMPRESSION}) {
        if (rProperties.Has(*p_variable)) {
            const double value = rProperties[*p_variable];
            KRATOS_ERROR_IF(!std::isfinite(value) || std::abs(value - initial_yield) > 1.0e-12 * initial_yield)
                << "Redundant yield stress conflicts with Swift initial yield: " << p_variable->Name();
        }
    }
    return 0;
}

void SmallStrainJ2IsotropicHardening3D::ValidatePlasticState(const Vector& rPlasticStrain, double P)
{
    KRATOS_ERROR_IF(rPlasticStrain.size() != 6 || !std::isfinite(P) || P < 0.0) << "Invalid J2 plastic state";
    double scale = 1.0;
    for (double value : rPlasticStrain) {
        KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Nonfinite J2 plastic strain";
        scale = std::max(scale, std::abs(value));
    }
    KRATOS_ERROR_IF(std::abs(rPlasticStrain[0] + rPlasticStrain[1] + rPlasticStrain[2]) > 1.0e-12 * scale)
        << "J2 plastic strain must be traceless";
}

void SmallStrainJ2IsotropicHardening3D::SetPlasticState(const Vector& rPlasticStrain, double P)
{
    ValidatePlasticState(rPlasticStrain, P);
    Vector copy = rPlasticStrain;
    GetPlasticStrain().swap(copy);
    mAccumulatedPlasticStrain = P;
    SetThreshold(0.0); // Derived cache is invalid, never used as authoritative history.
    SetPlasticDissipation(0.0);
}

void SmallStrainJ2IsotropicHardening3D::ValidateInitialState()
{
    if (!HasInitialState()) return;
    const auto& r_initial = GetInitialState();
    for (const auto* p_vector : {&r_initial.GetInitialStrainVector(), &r_initial.GetInitialStressVector()}) {
        KRATOS_ERROR_IF(p_vector->size() != 6) << "InitialState requires six 3D components";
        for (double value : *p_vector) KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Nonfinite InitialState";
    }
    const auto& r_f = r_initial.GetInitialDeformationGradientMatrix();
    bool zero = true, identity = r_f.size1() == 3 && r_f.size2() == 3;
    for (SizeType i = 0; i < r_f.size1(); ++i) {
        for (SizeType j = 0; j < r_f.size2(); ++j) {
            KRATOS_ERROR_IF_NOT(std::isfinite(r_f(i,j))) << "Nonfinite initial F";
            zero = zero && r_f(i,j) == 0.0;
            identity = identity && r_f(i,j) == (i == j ? 1.0 : 0.0);
        }
    }
    KRATOS_ERROR_IF(!zero && !identity) << "Nontrivial initial F is not supported";
}

void SmallStrainJ2IsotropicHardening3D::InitializeMaterial(const Properties& rProperties,
    const GeometryType& rGeometry, const Vector& rN)
{
    ProcessInfo info;
    Check(rProperties, rGeometry, info);
    ValidateInitialState();
    double threshold = 0.0, hardening = 0.0;
    const auto context = MakeEvolutionContext(rProperties);
    KRATOS_ERROR_IF_NOT(context.Hardening(rProperties, 0.0, threshold, hardening))
        << "Initial hardening is not representable";
    SetPlasticState(ZeroVector(6), 0.0);
    SetThreshold(threshold);
}

void SmallStrainJ2IsotropicHardening3D::ResetMaterial(const Properties& rProperties,
    const GeometryType& rGeometry, const Vector& rN)
{
    InitializeMaterial(rProperties, rGeometry, rN);
}

void SmallStrainJ2IsotropicHardening3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.clear();
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 6;
    rFeatures.mSpaceDimension = 3;
}

void SmallStrainJ2IsotropicHardening3D::ValidateResponse(Parameters& rValues)
{
    Check(rValues.GetMaterialProperties(), rValues.GetElementGeometry(), rValues.GetProcessInfo());
    KRATOS_ERROR_IF_NOT(rValues.GetOptions().Is(USE_ELEMENT_PROVIDED_STRAIN))
        << "J2 hardening requires USE_ELEMENT_PROVIDED_STRAIN; F fallback is not supported";
    KRATOS_ERROR_IF(rValues.GetOptions().Is(U_P_LAW)) << "J2 hardening does not support U_P_LAW";
    KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector()) << "Missing infinitesimal strain";
    KRATOS_ERROR_IF(rValues.GetStrainVector().size() != 6) << "Expected six engineering strain components";
    for (double value : rValues.GetStrainVector()) KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Nonfinite strain";
    ValidateInitialState();
}

bool SmallStrainJ2IsotropicHardening3D::EvaluateFlow(const array_1d<double,6>& rStress,
    double Q, Integrator::FlowLinearization& rFlow)
{
    if (!VonMisesPlasticPotential<6>::CalculatePlasticPotentialDerivativeAndJacobian(rStress, Q, rFlow.G, rFlow.DG)) return false;
    rFlow.F = rFlow.G; // Association belongs to this material family, not to the generic solver.
    return true;
}

bool SmallStrainJ2IsotropicHardening3D::EvaluateHistory(double OldP, double Lambda,
    const Integrator::FlowLinearization& rFlow, Integrator::HistoryResponse& rHistory)
{
    const double eta = AdvancedConstitutiveLawUtilities<6>::CalculateEquivalentPlasticStrainIncrement(rFlow.G);
    if (!std::isfinite(eta) || std::abs(eta - 1.0) > 1.0e-12) return false;
    rHistory.P = OldP + Lambda * eta;
    if (Lambda > 0.0 && rHistory.P == OldP) return false;
    rHistory.DPLambda = eta;
    for (SizeType j = 0; j < 6; ++j) {
        double derivative = 0.0;
        for (SizeType i = 0; i < 6; ++i)
            derivative += rFlow.DG(i,j) * (i < 3 ? 1.0 : 0.5) * rFlow.G[i];
        rHistory.DPStress[j] = Lambda * 2.0 / (3.0 * eta) * derivative;
    }
    return std::isfinite(rHistory.P);
}

SmallStrainJ2IsotropicHardening3D::Integrator::IntegrationResult
SmallStrainJ2IsotropicHardening3D::EvaluateTrial(Parameters& rValues, bool NeedTangent)
{
    ValidateResponse(rValues);
    Integrator::TrialInput input;
    input.EffectiveStrain = rValues.GetStrainVector();
    AddInitialStrainVectorContribution(input.EffectiveStrain);
    input.CommittedPlasticStrain = GetPlasticStrain();
    input.CommittedP = mAccumulatedPlasticStrain;
    CalculateElasticMatrix(input.ElasticMatrix, rValues); // Preserve the virtual elastic path.
    KRATOS_ERROR_IF(input.ElasticMatrix.size1() != 6 || input.ElasticMatrix.size2() != 6) << "Expected 6x6 elasticity";
    Matrix cholesky = ZeroMatrix(6,6);
    const double scale = norm_inf(input.ElasticMatrix);
    KRATOS_ERROR_IF(!std::isfinite(scale) || scale <= 0.0) << "Invalid elasticity";
    for (SizeType i = 0; i < 6; ++i) {
        for (SizeType j = 0; j <= i; ++j) {
            KRATOS_ERROR_IF(!std::isfinite(input.ElasticMatrix(i,j)) ||
                !std::isfinite(input.ElasticMatrix(j,i)) ||
                std::abs(input.ElasticMatrix(i,j) - input.ElasticMatrix(j,i)) > 1.0e-12 * scale)
                << "Elastic matrix must be finite and symmetric";
            double value = input.ElasticMatrix(i,j) / scale;
            for (SizeType k = 0; k < j; ++k) value -= cholesky(i,k) * cholesky(j,k);
            if (i == j) {
                KRATOS_ERROR_IF(value <= 0.0) << "Elastic matrix must be positive definite";
                cholesky(i,j) = std::sqrt(value);
            } else cholesky(i,j) = value / cholesky(j,j);
        }
    }
    const Vector elastic_strain = input.EffectiveStrain - input.CommittedPlasticStrain;
    noalias(input.PredictorStress) = prod(input.ElasticMatrix, elastic_strain);
    AddInitialStressVectorContribution(input.PredictorStress);
    const auto& r_properties = rValues.GetMaterialProperties();
    input.MaxIterations = r_properties.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_properties[MAX_NUMBER_NL_CL_ITERATIONS] : 100;
    input.NeedTangent = NeedTangent;
    const auto context = MakeEvolutionContext(r_properties);
    auto result = Integrator::IntegrateStressVectorWithHistory(input, context, rValues);
    if (!result.Candidate) {
        const char* reason = "unknown";
        switch (result.Failure) {
            case Integrator::IntegrationFailure::UnresolvedStress: reason = "unresolved stress scale"; break;
            case Integrator::IntegrationFailure::NonFiniteState: reason = "nonfinite or inadmissible state"; break;
            case Integrator::IntegrationFailure::SingularJacobian: reason = "singular Jacobian"; break;
            case Integrator::IntegrationFailure::LineSearch: reason = "line search failed"; break;
            case Integrator::IntegrationFailure::IterationLimit: reason = "iteration limit"; break;
            default: break;
        }
        KRATOS_ERROR << "J2 endpoint integration failed: " << reason << ", iterations=" << result.Iterations
            << ", stress residual=" << result.StressResidual << ", yield residual=" << result.YieldResidual
            << ", tolerance=" << result.StressTolerance;
    }
    ValidatePlasticState(result.Candidate->PlasticStrain, result.Candidate->P);
    return result;
}

void SmallStrainJ2IsotropicHardening3D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const bool stress = rValues.GetOptions().Is(COMPUTE_STRESS);
    const bool tangent = rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR);
    if (!stress && !tangent) {
        ValidateResponse(rValues);
        return;
    }
    if (stress) {
        KRATOS_ERROR_IF_NOT(rValues.IsSetStressVector()) << "Missing requested stress buffer";
        KRATOS_ERROR_IF(rValues.GetStressVector().size() != 6) << "Expected stress buffer of size six";
    }
    if (tangent) {
        KRATOS_ERROR_IF_NOT(rValues.IsSetConstitutiveMatrix()) << "Missing requested tangent buffer";
        KRATOS_ERROR_IF(rValues.GetConstitutiveMatrix().size1() != 6 || rValues.GetConstitutiveMatrix().size2() != 6)
            << "Expected 6x6 tangent buffer";
    }
    const auto result = EvaluateTrial(rValues, tangent);
    if (stress) noalias(rValues.GetStressVector()) = result.Candidate->Stress;
    if (tangent) noalias(rValues.GetConstitutiveMatrix()) = *result.Candidate->Tangent;
}

void SmallStrainJ2IsotropicHardening3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    const auto result = EvaluateTrial(rValues, false);
    const auto& r_candidate = *result.Candidate;
    // All validation/allocation precedes persistent writes. Finalize never publishes buffers.
    Vector plastic_strain = r_candidate.PlasticStrain;
    GetPlasticStrain().swap(plastic_strain);
    mAccumulatedPlasticStrain = r_candidate.P;
    SetThreshold(r_candidate.Threshold);
    SetPlasticDissipation(0.0);
}

bool SmallStrainJ2IsotropicHardening3D::Has(const Variable<double>& rVariable)
{ return rVariable == ACCUMULATED_PLASTIC_STRAIN; }
bool SmallStrainJ2IsotropicHardening3D::Has(const Variable<Vector>& rVariable)
{ return rVariable == PLASTIC_STRAIN_VECTOR; }
bool SmallStrainJ2IsotropicHardening3D::Has(const Variable<Matrix>& rVariable)
{ return rVariable == PLASTIC_STRAIN_TENSOR; }

double& SmallStrainJ2IsotropicHardening3D::GetValue(const Variable<double>& rVariable, double& rValue)
{
    KRATOS_ERROR_IF(rVariable != ACCUMULATED_PLASTIC_STRAIN) << "Unsupported J2 state variable: " << rVariable.Name();
    rValue = mAccumulatedPlasticStrain;
    return rValue;
}
Vector& SmallStrainJ2IsotropicHardening3D::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    KRATOS_ERROR_IF(rVariable != PLASTIC_STRAIN_VECTOR) << "Unsupported J2 state variable: " << rVariable.Name();
    rValue = GetPlasticStrain();
    return rValue;
}
Matrix& SmallStrainJ2IsotropicHardening3D::GetValue(const Variable<Matrix>& rVariable, Matrix& rValue)
{
    KRATOS_ERROR_IF(rVariable != PLASTIC_STRAIN_TENSOR) << "Unsupported J2 state variable: " << rVariable.Name();
    rValue = MathUtils<double>::StrainVectorToTensor(GetPlasticStrain());
    return rValue;
}
void SmallStrainJ2IsotropicHardening3D::SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rInfo)
{
    KRATOS_ERROR_IF(rVariable != ACCUMULATED_PLASTIC_STRAIN) << "Unsupported J2 state variable: " << rVariable.Name();
    SetPlasticState(GetPlasticStrain(), rValue);
}
void SmallStrainJ2IsotropicHardening3D::SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rInfo)
{
    KRATOS_ERROR_IF(rVariable != PLASTIC_STRAIN_VECTOR) << "Unsupported J2 state variable: " << rVariable.Name();
    SetPlasticState(rValue, mAccumulatedPlasticStrain);
}
void SmallStrainJ2IsotropicHardening3D::SetValue(const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rInfo)
{ KRATOS_ERROR << "Unsupported J2 state import: " << rVariable.Name() << "; use SetPlasticState"; }

double& SmallStrainJ2IsotropicHardening3D::CalculateValue(Parameters& rValues, const Variable<double>& rVariable, double& rValue)
{
    if (rVariable == YIELD_STRESS) {
        Check(rValues.GetMaterialProperties(), rValues.GetElementGeometry(), rValues.GetProcessInfo());
        double derivative = 0.0;
        const auto& r_properties = rValues.GetMaterialProperties();
        const auto context = MakeEvolutionContext(r_properties);
        KRATOS_ERROR_IF_NOT(context.Hardening(r_properties, mAccumulatedPlasticStrain, rValue, derivative))
            << "Committed hardening is not representable";
    } else if (rVariable == UNIAXIAL_STRESS) {
        rValue = EvaluateTrial(rValues, false).Candidate->EquivalentStress;
    } else return GetValue(rVariable, rValue);
    return rValue;
}
Vector& SmallStrainJ2IsotropicHardening3D::CalculateValue(Parameters& rValues, const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == STRESSES || rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR || rVariable == KIRCHHOFF_STRESS_VECTOR) {
        rValue = EvaluateTrial(rValues, false).Candidate->Stress;
    } else if (rVariable == STRAIN || rVariable == GREEN_LAGRANGE_STRAIN_VECTOR || rVariable == ALMANSI_STRAIN_VECTOR) {
        ValidateResponse(rValues);
        rValue = rValues.GetStrainVector();
    } else if (rVariable == INITIAL_STRAIN_VECTOR) {
        ValidateInitialState();
        rValue = HasInitialState() ? GetInitialState().GetInitialStrainVector() : Vector(ZeroVector(6));
    } else return GetValue(rVariable, rValue);
    return rValue;
}
Matrix& SmallStrainJ2IsotropicHardening3D::CalculateValue(Parameters& rValues, const Variable<Matrix>& rVariable, Matrix& rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX || rVariable == CONSTITUTIVE_MATRIX_PK2 || rVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        rValue = *EvaluateTrial(rValues, true).Candidate->Tangent;
    } else if (rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR) {
        const Vector stress = EvaluateTrial(rValues, false).Candidate->Stress;
        rValue = MathUtils<double>::StressVectorToTensor(stress);
    } else return GetValue(rVariable, rValue);
    return rValue;
}

void SmallStrainJ2IsotropicHardening3D::save(Serializer& rSerializer) const
{
    rSerializer.save("HardeningStateMarker", J2HardeningArchiveMarker);
    rSerializer.save("HardeningStateVersion", J2HardeningArchiveVersion);
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    rSerializer.save("AccumulatedPlasticStrain", mAccumulatedPlasticStrain);
}
void SmallStrainJ2IsotropicHardening3D::load(Serializer& rSerializer)
{
    std::uint64_t marker = 0;
    rSerializer.load("HardeningStateMarker", marker);
    KRATOS_ERROR_IF(marker != J2HardeningArchiveMarker) << "Unsupported J2 hardening archive format";
    int version = 0;
    double p = 0.0;
    rSerializer.load("HardeningStateVersion", version);
    KRATOS_ERROR_IF(version != J2HardeningArchiveVersion) << "Unsupported J2 hardening archive version";
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    rSerializer.load("AccumulatedPlasticStrain", p);
    KRATOS_ERROR_IF(GetPlasticDissipation() != 0.0 || !std::isfinite(GetThreshold()) || GetThreshold() < 0.0)
        << "Invalid J2 hardening base state";
    SetPlasticState(GetPlasticStrain(), p);
}
} // namespace Kratos
