# Kratos default BSD license: kratos/license.txt
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ConstitutiveLawsApplication as CLA


class TestSmallStrainJ2Hardening(KratosUnittest.TestCase):
    def test_existing_swift_material_and_lifecycle(self):
        model = KM.Model()
        part = model.CreateModelPart("Structure")
        materials = Path(__file__).parent / "SwiftPlasticity" / "Materials.json"
        settings = KM.Parameters('{"Parameters": {"materials_filename": ""}}')
        settings["Parameters"]["materials_filename"].SetString(str(materials))
        KM.ReadMaterialsUtility(settings, model)
        properties = part.GetProperties()[1]
        for variable in (CLA.HARDENING_CURVE, KM.FRACTURE_ENERGY, KM.YIELD_STRESS):
            self.assertFalse(properties.Has(variable))
        law = properties[KM.CONSTITUTIVE_LAW].Clone()
        self.assertEqual(law.Info(), "SmallStrainSwiftJ2Plasticity3DLaw")
        geometry = KM.Tetrahedra3D4(*[
            part.CreateNewNode(1, 0, 0, 0), part.CreateNewNode(2, 1, 0, 0),
            part.CreateNewNode(3, 0, 1, 0), part.CreateNewNode(4, 0, 0, 1)])
        shape = KM.Vector([0.25] * 4)
        self.assertEqual(law.Check(properties, geometry, part.ProcessInfo), 0)
        law.InitializeMaterial(properties, geometry, shape)
        values = KM.ConstitutiveLawParameters(geometry, properties, part.ProcessInfo)
        strain = KM.Vector([0.008, 0, 0, 0.004, 0, 0])
        stress = KM.Vector(6)
        tangent = KM.Matrix(6, 6)
        values.SetStrainVector(strain)
        values.SetStressVector(stress)
        values.SetConstitutiveMatrix(tangent)
        values.GetOptions().Set(KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        values.GetOptions().Set(KM.ConstitutiveLaw.COMPUTE_STRESS, True)
        values.GetOptions().Set(KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        law.CalculateMaterialResponseCauchy(values)
        self.assertEqual(law.GetValue(CLA.ACCUMULATED_PLASTIC_STRAIN), 0.0)
        law.FinalizeMaterialResponseCauchy(values)
        p = law.GetValue(CLA.ACCUMULATED_PLASTIC_STRAIN)
        self.assertGreater(p, 0.0)
        law.FinalizeMaterialResponseCauchy(values)
        self.assertEqual(law.GetValue(CLA.ACCUMULATED_PLASTIC_STRAIN), p)
        self.assertFalse(properties.Has(CLA.HARDENING_CURVE))

    def test_family_selection_and_import(self):
        law = CLA.SmallStrainJ2IsotropicHardening3DLaw()
        self.assertEqual(law.Clone().Info(), "SmallStrainJ2IsotropicHardening3DLaw")
        law.SetPlasticState(KM.Vector([0.002, -0.001, -0.001, 0, 0, 0]), 0.01)
        self.assertEqual(law.GetValue(CLA.ACCUMULATED_PLASTIC_STRAIN), 0.01)
        self.assertFalse(law.Has(KM.INTERNAL_VARIABLES))
        self.assertTrue(law.Has(CLA.PLASTIC_STRAIN_VECTOR))


if __name__ == "__main__":
    KratosUnittest.main()
