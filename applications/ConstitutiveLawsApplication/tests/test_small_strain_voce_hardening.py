# Kratos default BSD license: kratos/license.txt
from pathlib import Path
import math

import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestSmallStrainVoceHardening(KratosUnittest.TestCase):
    def test_material_definition_selection_and_no_swift_parameters(self):
        model = KM.Model()
        part = model.CreateModelPart("Structure")
        settings = KM.Parameters('{"Parameters": {"materials_filename": ""}}')
        settings["Parameters"]["materials_filename"].SetString(
            str(Path(__file__).parent / "VocePlasticity" / "Materials.json"))
        KM.ReadMaterialsUtility(settings, model)
        properties = part.GetProperties()[1]
        law = properties[KM.CONSTITUTIVE_LAW].Clone()
        geometry = KM.Tetrahedra3D4(*[part.CreateNewNode(i+1, *xyz) for i, xyz in
            enumerate(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))])
        for variable in (CLA.SWIFT_COEFFICIENT, CLA.SWIFT_INITIAL_STRAIN,
                         CLA.SWIFT_HARDENING_EXPONENT, KM.FRACTURE_ENERGY):
            self.assertFalse(properties.Has(variable))
        self.assertEqual(properties[CLA.HARDENING_CURVE], 8)
        self.assertEqual(law.Check(properties, geometry, part.ProcessInfo), 0)
        law.InitializeMaterial(properties, geometry, KM.Vector([0.25]*4))
        self.assertEqual(law.Info(), "SmallStrainJ2IsotropicHardening3DLaw")
        properties[KM.ISOTROPIC_HARDENING_MODULUS] = 0.0
        self.assertEqual(law.Check(properties, geometry, part.ProcessInfo), 0)
        properties[KM.ISOTROPIC_HARDENING_MODULUS] = 1.0
        with self.assertRaisesRegex(RuntimeError, "ISOTROPIC_HARDENING_MODULUS"):
            law.Check(properties, geometry, part.ProcessInfo)
        properties[KM.ISOTROPIC_HARDENING_MODULUS] = 0.0
        properties[CLA.HARDENING_CURVE] = 9
        with self.assertRaisesRegex(RuntimeError, "Unsupported HARDENING_CURVE"):
            law.Check(properties, geometry, part.ProcessInfo)

    def test_import_clone_and_committed_yield_query(self):
        model = KM.Model()
        part = model.CreateModelPart("Structure")
        properties = part.CreateNewProperties(1)
        for variable, value in ((KM.YOUNG_MODULUS, 210000.0), (KM.POISSON_RATIO, 0.3),
                                (CLA.HARDENING_CURVE, 8), (KM.YIELD_STRESS, 250.0),
                                (CLA.EXPONENTIAL_SATURATION_YIELD_STRESS, 450.0),
                                (KM.HARDENING_EXPONENT, 20.0)):
            properties[variable] = value
        geometry = KM.Tetrahedra3D4(*[part.CreateNewNode(i+1, *xyz) for i, xyz in
            enumerate(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))])
        law = CLA.SmallStrainJ2IsotropicHardening3DLaw()
        plastic = KM.Vector([0.002, -0.001, -0.001, 0.0, 0.0, 0.0])
        law.SetPlasticState(plastic, 0.01)
        values = KM.ConstitutiveLawParameters(geometry, properties, part.ProcessInfo)
        expected = 250.0 + 200.0 * (1.0-math.exp(-0.2))
        for current in (law, law.Clone()):
            self.assertAlmostEqual(current.CalculateValue(values, KM.YIELD_STRESS), expected, places=11)
            self.assertEqual(current.GetValue(CLA.ACCUMULATED_PLASTIC_STRAIN), 0.01)
            self.assertVectorAlmostEqual(current.GetValue(CLA.PLASTIC_STRAIN_VECTOR), plastic)


if __name__ == "__main__":
    KratosUnittest.main()
