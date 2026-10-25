// License: BSD License (Kratos default license: kratos/license.txt)
#pragma once
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "generic_small_strain_isotropic_plasticity.h"

namespace Kratos
{
/** @brief Isothermal associated 3D J2 endpoint plasticity with accumulated history.
 * The generic base archive/layout and all its legacy response paths are unchanged.
 * This family owns p; the inherited threshold is a cache (zero means invalid).
 * Requires element-provided infinitesimal engineering strain; no U-P/F fallback.
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainJ2IsotropicHardening3D
    : public GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<
        VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>
{
public:
    using Integrator = GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;
    using BaseType = GenericSmallStrainIsotropicPlasticity<Integrator>;
    using BaseType::Has;
    using BaseType::GetValue;
    using BaseType::SetValue;
    using BaseType::CalculateValue;
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainJ2IsotropicHardening3D);
    SmallStrainJ2IsotropicHardening3D() = default;
    SmallStrainJ2IsotropicHardening3D(const SmallStrainJ2IsotropicHardening3D&) = default;
    SmallStrainJ2IsotropicHardening3D& operator=(const SmallStrainJ2IsotropicHardening3D&) = default;
    ConstitutiveLaw::Pointer Clone() const override;
    std::string Info() const override { return "SmallStrainJ2IsotropicHardening3DLaw"; }
    int Check(const Properties&, const GeometryType&, const ProcessInfo&) const override;
    void InitializeMaterial(const Properties&, const GeometryType&, const Vector&) override;
    void ResetMaterial(const Properties&, const GeometryType&, const Vector&) override;
    void GetLawFeatures(Features&) override;
    void CalculateMaterialResponseCauchy(Parameters&) override;
    void FinalizeMaterialResponseCauchy(Parameters&) override;
    bool Has(const Variable<double>&) override;
    bool Has(const Variable<Vector>&) override;
    bool Has(const Variable<Matrix>&) override;
    double& GetValue(const Variable<double>&, double&) override;
    Vector& GetValue(const Variable<Vector>&, Vector&) override;
    Matrix& GetValue(const Variable<Matrix>&, Matrix&) override;
    void SetValue(const Variable<double>&, const double&, const ProcessInfo&) override;
    void SetValue(const Variable<Vector>&, const Vector&, const ProcessInfo&) override;
    void SetValue(const Variable<Matrix>&, const Matrix&, const ProcessInfo&) override;
    double& CalculateValue(Parameters&, const Variable<double>&, double&) override;
    Vector& CalculateValue(Parameters&, const Variable<Vector>&, Vector&) override;
    Matrix& CalculateValue(Parameters&, const Variable<Matrix>&, Matrix&) override;
    /// Atomic explicit import. p is never inferred from the plastic tensor.
    void SetPlasticState(const Vector& rPlasticStrain, double P);

protected:
    virtual Integrator::HardeningCurveType ResolveHardening(const Properties&) const;
    Integrator::IntegrationResult EvaluateTrial(Parameters&, bool NeedTangent);

private:
    double mAccumulatedPlasticStrain = 0.0;
    Integrator::EvolutionContext MakeEvolutionContext(const Properties&) const;
    void ValidateResponse(Parameters&);
    void ValidateInitialState();
    static void ValidatePlasticState(const Vector&, double);
    static bool EvaluateFlow(const array_1d<double,6>&, double, Integrator::FlowLinearization&);
    static bool EvaluateHistory(double, double, const Integrator::FlowLinearization&, Integrator::HistoryResponse&);
    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};
} // namespace Kratos
