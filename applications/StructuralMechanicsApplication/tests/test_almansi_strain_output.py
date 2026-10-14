"""Independent finite-strain output regression; Almansi fails on upstream."""

import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication  # noqa: F401
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    import KratosMultiphysics.ConstitutiveLawsApplication as CLA


def create_fixture(element_name, law_name="HyperElastic3DLaw"):
    model = KM.Model()
    model_part = model.CreateModelPart("cube")
    model_part.SetBufferSize(2)
    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
    model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
    model_part.ProcessInfo[KM.STEP] = 0
    coordinates = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    for node_id, xyz in enumerate(coordinates, 1):
        node = model_part.CreateNewNode(node_id, *xyz)
        for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                         KM.DISPLACEMENT_Z):
            node.AddDof(variable)
            node.Fix(variable)
    properties = model_part.GetProperties()[1]
    properties[KM.YOUNG_MODULUS] = 210000.0
    properties[KM.POISSON_RATIO] = 0.3
    properties[KM.DENSITY] = 1.0
    properties[KM.VOLUME_ACCELERATION] = [0.0, 0.0, 0.0]
    properties[KM.CONSTITUTIVE_LAW] = getattr(CLA, law_name)()
    element = model_part.CreateNewElement(
        element_name, 1, list(range(1, 9)), properties)
    element.Initialize(model_part.ProcessInfo)
    return model, model_part, element


def prescribe(model_part, deformation_gradient):
    """Move the mesh and set total displacement, preserving buffer slot 1."""
    for node in model_part.Nodes:
        reference = np.array([node.X0, node.Y0, node.Z0])
        current = deformation_gradient @ reference
        node.SetSolutionStepValue(KM.DISPLACEMENT, 0,
                                  (current - reference).tolist())
        node.X, node.Y, node.Z = current.tolist()


def analytical_strains(deformation_gradient):
    """Reference equations use no Kratos strain-output functions."""
    f = deformation_gradient
    almansi = 0.5 * (np.eye(3) - np.linalg.inv(f @ f.T))
    green_lagrange = 0.5 * (f.T @ f - np.eye(3))
    return almansi, green_lagrange


def engineering_voigt(tensor):
    return np.array([tensor[0, 0], tensor[1, 1], tensor[2, 2],
                     2 * tensor[0, 1], 2 * tensor[1, 2],
                     2 * tensor[0, 2]])


def serialized_state(model_part):
    serializer = KM.StreamSerializer()
    serializer.Save("cube", model_part)
    return serializer.__getstate__()[0]


def initialize_solution_step(model_part):
    KM.ResidualBasedIncrementalUpdateStaticScheme().InitializeSolutionStep(
        model_part, KM.CompressedMatrix(), KM.Vector(), KM.Vector())


class TestAlmansiStrainOutput(KratosUnittest.TestCase):
    def setUp(self):
        self.skipTestIfApplicationsNotAvailable("ConstitutiveLawsApplication")

    def _check_strain(self, formulation, measure):
        for name, f in (("axial", np.diag([1.004, 1.0, 1.0])),
                        ("shear", np.array([[1.0, 0.15, 0.0],
                                            [0.0, 1.0, 0.0],
                                            [0.0, 0.0, 1.0]]))):
            model, model_part, element = create_fixture(
                formulation + "LagrangianElement3D8N")
            model_part.CloneTimeStep(1.0)
            model_part.ProcessInfo[KM.STEP] = 1
            prescribe(model_part, f)
            initialize_solution_step(model_part)
            e, green = analytical_strains(f)
            expected = e if measure == "ALMANSI" else green
            for suffix in ("VECTOR", "TENSOR"):
                variable = getattr(KM, measure + "_STRAIN_" + suffix)
                actual = np.array(element.CalculateOnIntegrationPoints(
                    variable, model_part.ProcessInfo))
                reference = (engineering_voigt(expected)
                             if suffix == "VECTOR" else expected)
                with self.subTest(case=name, output=variable.Name()):
                    self.assertEqual(actual.shape[0], 8)
                    np.testing.assert_allclose(
                        actual, np.broadcast_to(reference, actual.shape),
                        rtol=0.0, atol=1e-12)

    def _check_repeatability(self, formulation):
        model, model_part, element = create_fixture(
            formulation + "LagrangianElement3D8N")
        model_part.CloneTimeStep(1.0)
        model_part.ProcessInfo[KM.STEP] = 1
        f = np.array([[1.0, 0.15, 0.0], [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])
        prescribe(model_part, f)
        initialize_solution_step(model_part)
        before = serialized_state(model_part)
        variables = (KM.ALMANSI_STRAIN_VECTOR, KM.ALMANSI_STRAIN_TENSOR,
                     KM.GREEN_LAGRANGE_STRAIN_VECTOR,
                     KM.GREEN_LAGRANGE_STRAIN_TENSOR,
                     KM.DEFORMATION_GRADIENT, KM.CAUCHY_STRESS_VECTOR)
        first = {v.Name(): np.array(element.CalculateOnIntegrationPoints(
            v, model_part.ProcessInfo)) for v in variables}
        for variable in reversed(variables):
            repeated = np.array(element.CalculateOnIntegrationPoints(
                variable, model_part.ProcessInfo))
            np.testing.assert_array_equal(first[variable.Name()], repeated)
        np.testing.assert_allclose(
            first["DEFORMATION_GRADIENT"], np.broadcast_to(f, (8, 3, 3)),
            rtol=0.0, atol=1e-12)
        self.assertEqual(before, serialized_state(model_part))

    def test_tl_almansi(self):
        self._check_strain("Total", "ALMANSI")

    def test_ul_almansi(self):
        self._check_strain("Updated", "ALMANSI")

    def test_tl_green_lagrange(self):
        self._check_strain("Total", "GREEN_LAGRANGE")

    def test_ul_green_lagrange(self):
        self._check_strain("Updated", "GREEN_LAGRANGE")

    def test_tl_repeatability(self):
        self._check_repeatability("Total")

    def test_ul_repeatability(self):
        self._check_repeatability("Updated")


if __name__ == "__main__":
    KratosUnittest.main()
