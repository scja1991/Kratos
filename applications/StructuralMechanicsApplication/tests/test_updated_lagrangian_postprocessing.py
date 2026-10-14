"""Regression for stress output under unchanged prescribed finite deformation."""
import math

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication  # Register elements.
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    import KratosMultiphysics.ConstitutiveLawsApplication as CLA


def create_case(element_name):
    model = KM.Model()
    mp = model.CreateModelPart("cube", 2)
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
    mp.ProcessInfo[KM.DOMAIN_SIZE] = 3
    mp.ProcessInfo[KM.STEP] = 0
    mp.ProcessInfo[KM.IS_RESTARTED] = False
    properties = mp.CreateNewProperties(1)
    properties[KM.YOUNG_MODULUS] = 210000.0  # MPa; unit cube in mm.
    properties[KM.POISSON_RATIO] = 0.3
    properties[KM.DENSITY] = 1.0
    properties[KM.CONSTITUTIVE_LAW] = CLA.HyperElastic3DLaw()
    coordinates = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    for node_id, xyz in enumerate(coordinates, 1):
        node = mp.CreateNewNode(node_id, *xyz)
        for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                         KM.DISPLACEMENT_Z):
            node.AddDof(variable)
            node.Fix(variable)
    element = mp.CreateNewElement(element_name, 1, list(range(1, 9)), properties)
    element.Initialize(mp.ProcessInfo)
    mp.Check()
    return model, mp, element


def prescribe_stretch(mp, step, stretch):
    mp.CloneTimeStep(float(step))
    mp.ProcessInfo[KM.STEP] = step
    KM.Scheme().InitializeSolutionStep(
        mp, KM.CompressedMatrix(), KM.Vector(), KM.Vector())
    for node in mp.Nodes:
        node.SetSolutionStepValue(KM.DISPLACEMENT,
                                  [(stretch - 1.0) * node.X0, 0.0, 0.0])
        node.X = stretch * node.X0
        node.Y = node.Y0
        node.Z = node.Z0


def finalize_solution_step(mp):
    KM.Scheme().FinalizeSolutionStep(
        mp, KM.CompressedMatrix(), KM.Vector(), KM.Vector())


def analytical_cauchy(stretch):
    mu = 210000.0 / (2.0 * 1.3)
    lame_lambda = 210000.0 * 0.3 / (1.3 * 0.4)
    transverse = lame_lambda * math.log(stretch) / stretch
    return [mu * (stretch * stretch - 1.0) / stretch + transverse,
            transverse, transverse, 0.0, 0.0, 0.0]


def assembled_cauchy(mp, element, stretch):
    """Recover homogeneous Cauchy stress from assembled current face forces."""
    lhs, rhs = KM.Matrix(), KM.Vector()
    element.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
    stress = []
    for axis in range(3):
        face_force = sum(rhs[3 * i + axis] for i, node in enumerate(mp.Nodes)
                         if (node.X0, node.Y0, node.Z0)[axis] == 1.0)
        area = 1.0 if axis == 0 else stretch
        stress.append(-face_force / area)
    return stress


class TestUpdatedLagrangianPostprocessing(KratosUnittest.TestCase):
    def _check_stress_output(self, element_name):
        self.skipTestIfApplicationsNotAvailable("ConstitutiveLawsApplication")
        model, mp, element = create_case(element_name)
        for step, stretch in enumerate((1.002, 1.004), 1):
            prescribe_stretch(mp, step, stretch)
            expected = analytical_cauchy(stretch)
            self.assertVectorAlmostEqual(
                assembled_cauchy(mp, element, stretch), expected[:3], places=7)
            for stress in element.CalculateOnIntegrationPoints(
                    KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo):
                self.assertVectorAlmostEqual(stress, expected, places=7)
            finalize_solution_step(mp)
            for repeat in range(3):
                with self.subTest(step=step, repeat=repeat):
                    for stress in element.CalculateOnIntegrationPoints(
                            KM.CAUCHY_STRESS_VECTOR, mp.ProcessInfo):
                        self.assertVectorAlmostEqual(stress, expected, places=7)

    def test_updated_lagrangian_stress_after_finalization(self):
        self._check_stress_output("UpdatedLagrangianElement3D8N")

    def test_total_lagrangian_stress_after_finalization(self):
        self._check_stress_output("TotalLagrangianElement3D8N")


if __name__ == "__main__":
    KratosUnittest.main()
