# Kratos default BSD license: kratos/license.txt
"""Second consumer of the homogeneous J2 FE/restart regression."""
import math

import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_small_strain_swift_fe as swift_fe


def voce_reference_threshold(p):
    # Independent scalar reference, no production constitutive-law calls.
    return 250.0 - 200.0 * math.expm1(-20.0 * p)


class TestSmallStrainVoceFE(swift_fe.TestSmallStrainSwiftFE):
    material_folder = "VocePlasticity"
    law_name = "SmallStrainJ2IsotropicHardening3DLaw"
    report_variable = "KRATOS_VOCE_FE_REPORT"
    material_values = ((KM.YOUNG_MODULUS, 210000.0), (KM.POISSON_RATIO, 0.3),
                       (CLA.HARDENING_CURVE, 8), (KM.YIELD_STRESS, 250.0),
                       (CLA.EXPONENTIAL_SATURATION_YIELD_STRESS, 450.0),
                       (KM.HARDENING_EXPONENT, 20.0))
    absent_properties = (KM.FRACTURE_ENERGY, CLA.SWIFT_COEFFICIENT,
                         CLA.SWIFT_INITIAL_STRAIN, CLA.SWIFT_HARDENING_EXPONENT)

    def _reference(self):
        return swift_fe.uniaxial_reference(voce_reference_threshold)


if __name__ == "__main__":
    KratosUnittest.main()
