# Kratos default BSD license: kratos/license.txt
"""Homogeneous one-element lifecycle/restart regression; mm, N, MPa."""
import json
import math
import os
from pathlib import Path
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication  # Element registration
import KratosMultiphysics.ConstitutiveLawsApplication as CLA


# Fixed before FE evaluation. No productive hardening/integration calls in reference.
YOUNG = 210000.0
POISSON = 0.3
COEFFICIENT = 1000.0
INITIAL_STRAIN = 0.01
EXPONENT = 0.2
AXIAL_PATH = (0.0005, 0.0015, 0.002, 0.003, 0.0045, 0.006,
              0.005, 0.0045, 0.0055, 0.006, 0.0065, 0.008)
RESTART_STEP = 6
RELATIVE_TOLERANCE = 1.0e-9
STRESS_ABSOLUTE = 1.0e-6       # MPa: local integration floor ~1.5e-8 MPa plus global error.
REACTION_ABSOLUTE = 1.0e-6     # N: unit initial area, same stress-to-force scale.
STATE_ABSOLUTE = 1.0e-11      # Strain or mm: stress error / E and roundoff allowance.
FREE_RESIDUAL_LIMIT = 1.0e-7  # N: independently assembled free-DOF L2 norm.
DISPLACEMENT_DOFS = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)


def uniaxial_reference():
    """Signed axial plastic strain and accumulated p are independent histories."""
    p = axial_plastic = 0.0
    result = []
    for strain in AXIAL_PATH:
        trial = YOUNG * (strain - axial_plastic)
        yield_stress = COEFFICIENT * (INITIAL_STRAIN + p)**EXPONENT
        delta_p = 0.0
        if abs(trial) > yield_stress:
            lower, upper = 0.0, (abs(trial) - yield_stress) / YOUNG
            for _ in range(120):
                mid = 0.5 * (lower + upper)
                residual = (abs(trial) - YOUNG * mid
                            - COEFFICIENT * (INITIAL_STRAIN + p + mid)**EXPONENT)
                if residual > 0.0:
                    lower = mid
                else:
                    upper = mid
            delta_p = 0.5 * (lower + upper)
            axial_plastic += math.copysign(delta_p, trial)
            p += delta_p
        stress = YOUNG * (strain - axial_plastic)
        transverse = -POISSON * stress / YOUNG - axial_plastic / 2.0
        result.append({"axial_strain": strain, "stress": stress, "p": p,
                       "axial_plastic_strain": axial_plastic, "delta_p": delta_p,
                       "transverse_strain": transverse})
    return result


class TestSmallStrainSwiftFE(KratosUnittest.TestCase):
    def _create_model(self):
        model = KM.Model()
        part = model.CreateModelPart("Structure")
        part.SetBufferSize(2)
        for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION):
            part.AddNodalSolutionStepVariable(variable)
        part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        part.ProcessInfo[KM.IS_RESTARTED] = False
        part.ProcessInfo[KM.STEP] = 0
        part.ProcessInfo[KM.TIME] = 0.0
        settings = KM.Parameters('{"Parameters": {"materials_filename": ""}}')
        settings["Parameters"]["materials_filename"].SetString(
            str(Path(__file__).parent / "SwiftPlasticity" / "Materials.json"))
        KM.ReadMaterialsUtility(settings, model)
        properties = part.GetProperties()[1]
        for variable, value in ((KM.YOUNG_MODULUS, YOUNG), (KM.POISSON_RATIO, POISSON),
                                (CLA.SWIFT_COEFFICIENT, COEFFICIENT),
                                (CLA.SWIFT_INITIAL_STRAIN, INITIAL_STRAIN),
                                (CLA.SWIFT_HARDENING_EXPONENT, EXPONENT)):
            self.assertEqual(properties[variable], value)
        self.assertEqual(properties[KM.CONSTITUTIVE_LAW].Info(), "SmallStrainSwiftJ2Plasticity3DLaw")
        for variable in (CLA.HARDENING_CURVE, KM.FRACTURE_ENERGY, KM.YIELD_STRESS):
            self.assertFalse(properties.Has(variable))
        coordinates = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                       (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
        for node_id, xyz in enumerate(coordinates, 1):
            node = part.CreateNewNode(node_id, *xyz)
            for displacement, reaction in zip(DISPLACEMENT_DOFS,
                                               (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z)):
                node.AddDof(displacement, reaction)
            if node.X0 == 0.0 or node.X0 == 1.0:
                node.Fix(KM.DISPLACEMENT_X)
            if node.Y0 == 0.0:
                node.Fix(KM.DISPLACEMENT_Y)
            if node.Z0 == 0.0:
                node.Fix(KM.DISPLACEMENT_Z)
        part.CreateNewElement("SmallDisplacementElement3D8N", 1, list(range(1, 9)), properties)
        self.assertEqual(sum(not n.IsFixed(v) for n in part.Nodes for v in DISPLACEMENT_DOFS), 8)
        return model

    def _create_solver(self, part):
        scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = KM.SkylineLUFactorizationSolver()
        builder = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
        criterion = KM.ResidualCriteria(1.0e-10, 1.0e-9)
        echo_level = int(bool(os.environ.get("KRATOS_SWIFT_FE_REPORT")))
        criterion.SetEchoLevel(echo_level)
        strategy = KM.ResidualBasedNewtonRaphsonStrategy(
            part, scheme, criterion, builder, 25, True, False, False)
        strategy.SetEchoLevel(echo_level)
        strategy.Initialize()
        strategy.Check()
        return strategy

    def _history(self, part):
        element = part.Elements[1]
        return {
            "p": list(element.CalculateOnIntegrationPoints(CLA.ACCUMULATED_PLASTIC_STRAIN, part.ProcessInfo)),
            "plastic_strain": [list(v) for v in element.CalculateOnIntegrationPoints(
                CLA.PLASTIC_STRAIN_VECTOR, part.ProcessInfo)]}

    def _snapshot(self, part):
        element = part.Elements[1]
        history = self._history(part)
        stress = [list(v) for v in element.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)]
        strain = [list(v) for v in element.CalculateOnIntegrationPoints(KM.GREEN_LAGRANGE_STRAIN_VECTOR, part.ProcessInfo)]
        self.assertEqual(self._history(part), history)  # Element output remains a trial query.
        return {**history, "stress": stress, "strain": strain,
                "displacement": [list(n.GetSolutionStepValue(KM.DISPLACEMENT)) for n in part.Nodes],
                "reaction": [list(n.GetSolutionStepValue(KM.REACTION)) for n in part.Nodes]}

    def _free_residual(self, part):
        rhs = KM.Vector()
        element = part.Elements[1]
        element.CalculateRightHandSide(rhs, part.ProcessInfo)
        values = [rhs[3*i+j] for i, node in enumerate(element.GetGeometry())
                  for j, variable in enumerate(DISPLACEMENT_DOFS) if not node.IsFixed(variable)]
        return math.sqrt(sum(v*v for v in values))

    def _near(self, actual, expected, absolute):
        self.assertLessEqual(abs(actual - expected), absolute + RELATIVE_TOLERANCE * abs(expected))

    def _check_reference(self, part, state, reference, committed):
        self.assertEqual(len(state["stress"]), 8)
        expected_stress = [reference["stress"], 0.0, 0.0, 0.0, 0.0, 0.0]
        ep = reference["axial_plastic_strain"]
        expected_plastic = [ep, -ep/2, -ep/2, 0.0, 0.0, 0.0]
        expected_strain = [reference["axial_strain"], reference["transverse_strain"],
                           reference["transverse_strain"], 0.0, 0.0, 0.0]
        for j, expected in enumerate(expected_stress):
            values = [v[j] for v in state["stress"]]
            self._near(max(values) - min(values), 0.0,
                       STRESS_ABSOLUTE + RELATIVE_TOLERANCE * abs(expected))
        for field in ("strain", "plastic_strain") if committed else ("strain",):
            for j in range(6):
                values = [v[j] for v in state[field]]
                self._near(max(values) - min(values), 0.0, STATE_ABSOLUTE)
        if committed:
            self._near(max(state["p"]) - min(state["p"]), 0.0, STATE_ABSOLUTE)
        for point in range(8):
            for actual, expected in zip(state["stress"][point], expected_stress):
                self._near(actual, expected, STRESS_ABSOLUTE)
            for actual, expected in zip(state["strain"][point], expected_strain):
                self._near(actual, expected, STATE_ABSOLUTE)
            if committed:
                self._near(state["p"][point], reference["p"], STATE_ABSOLUTE)
                for actual, expected in zip(state["plastic_strain"][point], expected_plastic):
                    self._near(actual, expected, STATE_ABSOLUTE)
        for node, displacement in zip(part.Nodes, state["displacement"]):
            expected = [reference["axial_strain"] * node.X0,
                        reference["transverse_strain"] * node.Y0,
                        reference["transverse_strain"] * node.Z0]
            for actual, value in zip(displacement, expected):
                self._near(actual, value, STATE_ABSOLUTE)
        # KM.REACTION = internal - external. Positive tension: +sigma*A at x=1,
        # -sigma*A at x=0; A is the ORIGINAL unit area, never the contracted area.
        left = sum(n.GetSolutionStepValue(KM.REACTION_X) for n in part.Nodes if n.X0 == 0.0)
        right = sum(n.GetSolutionStepValue(KM.REACTION_X) for n in part.Nodes if n.X0 == 1.0)
        self._near(left, -reference["stress"], REACTION_ABSOLUTE)
        self._near(right, reference["stress"], REACTION_ABSOLUTE)
        self._near(left + right, 0.0, REACTION_ABSOLUTE)
        for node in part.Nodes:
            reaction = node.GetSolutionStepValue(KM.REACTION)
            for j, variable in enumerate(DISPLACEMENT_DOFS):
                if not node.IsFixed(variable):
                    self._near(reaction[j], 0.0, REACTION_ABSOLUTE)
        return left, right

    def _solve_step(self, part, strategy, step, reference):
        part.CloneTimeStep(float(step))
        part.ProcessInfo[KM.STEP] = step
        for node in part.Nodes:
            if node.X0 == 1.0:
                node.SetSolutionStepValue(KM.DISPLACEMENT_X, reference["axial_strain"])
        strategy.InitializeSolutionStep()
        strategy.Predict()
        before = self._history(part)
        initial_residual = self._free_residual(part)
        converged = strategy.SolveSolutionStep()
        self.assertTrue(converged, f"Step {step} did not converge; no finalize is allowed")
        final_residual = self._free_residual(part)
        self.assertLessEqual(final_residual, FREE_RESIDUAL_LIMIT)
        trial = self._snapshot(part)
        self.assertEqual({key: trial[key] for key in before}, before)
        left, right = self._check_reference(part, trial, reference, committed=False)
        iterations = part.ProcessInfo[KM.NL_ITERATION_NUMBER]
        criterion_norm = part.ProcessInfo[KM.RESIDUAL_NORM]
        criterion_ratio = part.ProcessInfo[KM.CONVERGENCE_RATIO]
        strategy.FinalizeSolutionStep()
        committed = self._snapshot(part)
        self._check_reference(part, committed, reference, committed=True)
        if step in (7, 8, 9, 10):
            self.assertEqual(committed["p"], before["p"])
        if step in (11, 12):
            self.assertTrue(all(new > old for new, old in zip(committed["p"], before["p"])))
        return {"step": step, "time": part.ProcessInfo[KM.TIME], "reference": reference,
                "converged": bool(converged), "iterations": iterations,
                "initial_free_residual_l2": initial_residual, "final_free_residual_l2": final_residual,
                "criterion_norm": criterion_norm, "criterion_ratio": criterion_ratio,
                "left_reaction": left, "right_reaction": right,
                "trial": trial, "committed": committed}

    def _restart_signature(self, part):
        return {"state": self._snapshot(part), "buffer_size": part.GetBufferSize(),
                "step": part.ProcessInfo[KM.STEP], "time": part.ProcessInfo[KM.TIME],
                "dt": part.ProcessInfo[KM.DELTA_TIME],
                "nodes": [{"id": n.Id, "coordinates": [n.X0, n.Y0, n.Z0],
                           "fixity": [n.IsFixed(v) for v in DISPLACEMENT_DOFS],
                           "displacement_buffer": [list(n.GetSolutionStepValue(KM.DISPLACEMENT, i)) for i in range(2)],
                           "reaction_buffer": [list(n.GetSolutionStepValue(KM.REACTION, i)) for i in range(2)]}
                          for n in part.Nodes],
                "properties": [part.GetProperties()[1][v] for v in (
                    KM.YOUNG_MODULUS, KM.POISSON_RATIO, CLA.SWIFT_COEFFICIENT,
                    CLA.SWIFT_INITIAL_STRAIN, CLA.SWIFT_HARDENING_EXPONENT)],
                "law": part.GetProperties()[1][KM.CONSTITUTIVE_LAW].Info()}

    def _error_maxima(self, part, records):
        errors = {key: 0.0 for key in ("axial_stress_mpa", "other_stress_mpa", "stress_spread_mpa",
                  "p", "plastic_strain", "strain", "displacement_mm", "axial_reaction_n", "free_residual_n")}
        for record in records:
            reference = record["reference"]
            ep = reference["axial_plastic_strain"]
            expected_plastic = [ep, -ep/2, -ep/2, 0.0, 0.0, 0.0]
            expected_strain = [reference["axial_strain"], reference["transverse_strain"],
                               reference["transverse_strain"], 0.0, 0.0, 0.0]
            for phase in ("trial", "committed"):
                state = record[phase]
                for point in range(8):
                    errors["axial_stress_mpa"] = max(errors["axial_stress_mpa"],
                        abs(state["stress"][point][0] - reference["stress"]))
                    errors["other_stress_mpa"] = max(errors["other_stress_mpa"],
                        max(abs(v) for v in state["stress"][point][1:]))
                    errors["strain"] = max(errors["strain"],
                        max(abs(x-y) for x, y in zip(state["strain"][point], expected_strain)))
                    if phase == "committed":
                        errors["p"] = max(errors["p"], abs(state["p"][point] - reference["p"]))
                        errors["plastic_strain"] = max(errors["plastic_strain"],
                            max(abs(x-y) for x, y in zip(state["plastic_strain"][point], expected_plastic)))
                errors["stress_spread_mpa"] = max(errors["stress_spread_mpa"],
                    max(max(v[j] for v in state["stress"]) - min(v[j] for v in state["stress"]) for j in range(6)))
                for node, displacement in zip(part.Nodes, state["displacement"]):
                    expected = [reference["axial_strain"]*node.X0, reference["transverse_strain"]*node.Y0,
                                reference["transverse_strain"]*node.Z0]
                    errors["displacement_mm"] = max(errors["displacement_mm"],
                        max(abs(x-y) for x, y in zip(displacement, expected)))
            errors["axial_reaction_n"] = max(errors["axial_reaction_n"],
                abs(record["right_reaction"] - reference["stress"]),
                abs(record["left_reaction"] + reference["stress"]))
            errors["free_residual_n"] = max(errors["free_residual_n"], record["final_free_residual_l2"])
        return errors

    def test_uniaxial_path_and_model_restart(self):
        reference = uniaxial_reference()  # The entire reference exists before FE.
        model = self._create_model()
        part = model["Structure"]
        strategy = self._create_solver(part)
        uninterrupted = []
        with tempfile.TemporaryDirectory(prefix="swift_fe_restart_") as folder:
            archive = str(Path(folder) / "checkpoint")
            for step, expected in enumerate(reference, 1):
                uninterrupted.append(self._solve_step(part, strategy, step, expected))
                if step == RESTART_STEP:
                    checkpoint = self._restart_signature(part)
                    self.assertTrue(all(p > 0.0 for p in checkpoint["state"]["p"]))
                    KM.FileSerializer(archive, KM.SerializerTraceType.SERIALIZER_NO_TRACE).Save("Model", model)
            restored_model = KM.Model()
            KM.FileSerializer(archive, KM.SerializerTraceType.SERIALIZER_NO_TRACE).Load("Model", restored_model)
            restored = restored_model["Structure"]
            self.assertEqual(self._restart_signature(restored), checkpoint)
            restored.ProcessInfo[KM.IS_RESTARTED] = True
            restored.ProcessInfo[KM.LOAD_RESTART] = RESTART_STEP + 1
            restored_strategy = self._create_solver(restored)
            self.assertEqual(self._restart_signature(restored), checkpoint)  # Initialize must not reset IP laws.
            continued = []
            for step in range(RESTART_STEP + 1, len(reference) + 1):
                record = self._solve_step(restored, restored_strategy, step, reference[step - 1])
                continued.append(record)
                original = uninterrupted[step - 1]
                for phase in ("trial", "committed"):
                    for key, absolute in (("stress", STRESS_ABSOLUTE), ("strain", STATE_ABSOLUTE),
                                          ("plastic_strain", STATE_ABSOLUTE), ("p", STATE_ABSOLUTE),
                                          ("displacement", STATE_ABSOLUTE), ("reaction", REACTION_ABSOLUTE)):
                        a, b = original[phase][key], record[phase][key]
                        if key == "p":
                            for x, y in zip(a, b):
                                self._near(x, y, absolute)
                        else:
                            for row_a, row_b in zip(a, b):
                                for x, y in zip(row_a, row_b):
                                    self._near(x, y, absolute)
        restart_errors = {key: 0.0 for key in ("stress", "strain", "plastic_strain", "p", "displacement", "reaction")}
        for record in continued:
            original = uninterrupted[record["step"] - 1]
            for phase in ("trial", "committed"):
                for key in restart_errors:
                    a, b = original[phase][key], record[phase][key]
                    if key != "p":
                        a = [x for row in a for x in row]
                        b = [x for row in b for x in row]
                    restart_errors[key] = max(restart_errors[key], max(abs(x-y) for x, y in zip(a, b)))
        output = os.environ.get("KRATOS_SWIFT_FE_REPORT")
        if output:
            Path(output).write_text(json.dumps({"element": "SmallDisplacementElement3D8N",
                "units": {"length": "mm", "force": "N", "stress": "MPa"},
                "tolerances": {"relative": RELATIVE_TOLERANCE, "stress": STRESS_ABSOLUTE,
                    "reaction": REACTION_ABSOLUTE, "state": STATE_ABSOLUTE,
                    "free_residual_l2": FREE_RESIDUAL_LIMIT},
                "phase_semantics": {"trial": "Converged stress/strain; p and plastic_strain still committed input",
                                    "committed": "After strategy.FinalizeSolutionStep"},
                "restart_step": RESTART_STEP, "checkpoint": checkpoint,
                "max_errors": self._error_maxima(part, uninterrupted + continued),
                "max_restart_errors": restart_errors,
                "uninterrupted": uninterrupted, "restarted": continued}, indent=2) + "\n")
        restored_strategy.Clear()
        strategy.Clear()


if __name__ == "__main__":
    KratosUnittest.main()
