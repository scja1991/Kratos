// License: BSD License (Kratos default license: kratos/license.txt)
#pragma once
#include "small_strain_j2_isotropic_hardening_3d.h"
namespace Kratos
{
/// Compatibility name: Swift is implicit; a conflicting explicit curve is rejected.
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainSwiftJ2Plasticity3D
    : public SmallStrainJ2IsotropicHardening3D
{
public:
    using BaseType = SmallStrainJ2IsotropicHardening3D;
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainSwiftJ2Plasticity3D);
    ConstitutiveLaw::Pointer Clone() const override { return Kratos::make_shared<SmallStrainSwiftJ2Plasticity3D>(*this); }
    std::string Info() const override { return "SmallStrainSwiftJ2Plasticity3DLaw"; }
protected:
    Integrator::HardeningCurveType ResolveHardening(const Properties& rProperties) const override
    {
        KRATOS_ERROR_IF(rProperties.Has(HARDENING_CURVE) &&
            rProperties[HARDENING_CURVE] != static_cast<int>(Integrator::HardeningCurveType::Swift))
            << "Swift wrapper: conflicting HARDENING_CURVE";
        return Integrator::HardeningCurveType::Swift;
    }
private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override { KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType); }
    void load(Serializer& rSerializer) override { KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType); }
};
} // namespace Kratos
