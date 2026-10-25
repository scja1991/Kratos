# Local small-strain Swift prototype

Based on public upstream master `8b7484c2ae79451b8c83e6831865f6641dc62018`.
This is a local architecture proposal, without maintainer approval.

`SmallStrainSwiftJ2Plasticity3DLaw` accepts the unchanged `Materials.json`:
`YOUNG_MODULUS`, `POISSON_RATIO`, `SWIFT_COEFFICIENT`,
`SWIFT_INITIAL_STRAIN`, and `SWIFT_HARDENING_EXPONENT`.
`HARDENING_CURVE` is optional for this wrapper; if supplied it must be 7.
The shared `SmallStrainJ2IsotropicHardening3DLaw` requires that explicit selection.
This branch also supports Voce through selection 8 in the shared family; see
[VocePlasticity/README.md](../VocePlasticity/README.md). Other selections are rejected.

The shared family owns accumulated plastic strain. Its explicit opt-in integrator
solves six stress equations and yield consistency with a total plastic multiplier.
The analytic tangent differentiates these same endpoint equations with six linear
right-hand sides. Elasticity is supplied through the existing virtual elastic path.
The legacy dissipation-based integrator and base archive format are unchanged.

The first scope requires element-provided infinitesimal 3D engineering strain,
constant isothermal material parameters, associated J2 flow, and the analytic tangent
(default, or TANGENT_OPERATOR_ESTIMATION=0). U-P, F-based input, reduced dimensions,
thermal coupling, damage, kinematic hardening and automatic substeps are unsupported.
Initial stress and strain are supported through InitialState; repeated trials do not
modify either the input strain or persistent history. Finalize recomputes and commits
only a converged state, independently of output flags, without publishing buffers.

ACCUMULATED_PLASTIC_STRAIN and PLASTIC_STRAIN_VECTOR expose committed state;
PLASTIC_STRAIN_TENSOR is a read-only view. SetPlasticState imports both independent
quantities atomically. Neither legacy INTERNAL_VARIABLES nor PLASTIC_DISSIPATION,
projected EQUIVALENT_PLASTIC_STRAIN or an unspecified energy are exposed. Import/load
invalidates the derived threshold cache (zero marker). Restart requires matching
material properties and must not be followed by InitializeMaterial (which resets).

The new family archive format is **version 2**. Its payload order is:

1. `HardeningStateMarker`: fixed `uint64_t`, `0x7ff84a3248415244`.
2. `HardeningStateVersion`: integer `2`.
3. Unchanged generic-base serialization: `ConstitutiveLaw` (Flags and InitialState),
   `PlasticDissipation`, `Threshold`, and the engineering `PlasticStrain` vector.
4. Independent `AccumulatedPlasticStrain`.

The Swift wrapper retains its base-class serialization frame and adds no fields.
Marker and version are checked before the first base payload is loaded. New binary
and precise Trace archives roundtrip with restart continuation. A wrong marker or
unsupported version leaves an already plasticized target unchanged; header-only
fixtures with no base payload verify that rejection does not depend on reading it.

The following formats are explicitly rejected, without automatic migration:

| Reconstructed writer / source | Historical serialization chain |
| --- | --- |
| Original standalone Swift, commit `ab168b6f962`, `small_strain_swift_j2_plasticity_3d.cpp`, save at line 283 | Swift -> ConstitutiveLaw -> plastic vector -> p |
| Later generic Swift, commit `b2dde3d44cc`, generic header save at line 477 and wrapper save | Wrapper -> historical generic base -> ConstitutiveLaw -> dissipation, threshold, plastic vector, p |
| Immediately preceding uncommitted local prototype, version 1 | Wrapper -> family -> master generic base -> marker, version 1, p |

These are reconstructed writers from verified save sequences and complete base
chains, **not actual previously saved archive files**. Both binary and Trace
rejection are tested. Trace rejects the first incompatible tag; binary rejects
the missing prefix marker. A prefix carrying version 1 or 99 is also rejected.
The preceding local version-1 prototype archives must be regenerated; they are
not accepted as version 2. Existing master material-law serializers, including
both plane-stress chains fixed by #14895, remain unchanged.

The early rejection guarantee covers these specified layouts and incompatible
headers. It is not general validation of arbitrarily damaged archives. A matching
header followed by a damaged base payload still uses the existing core reader;
transactional recovery from every possible payload failure is not promised.

The C++ tests use an independent scalar bisection and full tensor contractions.
The parameterized tangent table contains nine paths for each exponent
`n=0, 0.2, 1, 2`: elastic, uniaxial-stress tension/compression, pure shear,
unloading, renewed yielding, nonproportional history, and either side of initial
yield. `K=400/0.01^n` gives the same initial yield across these comparisons.
The uniaxial-stress endpoints are constructed by independent scalar bisection.

Every case compares the complete matrix from external stress differences and an
independently evaluated directional difference at fixed committed history. The
independent elastic predictor, scalar Swift formula, and a roundoff guard reject
perturbations crossing the branch. No unique two-sided derivative is required at
the switch. Nine h values from `1e-5` to `1e-9` are considered; near initial yield
only the five branch-faithful values `1e-7, 3e-8, 1e-8, 3e-9, 1e-9` are accepted.
A step passes only when BOTH matrix and direction errors pass, and at least two
adjacent steps must pass per case. Committed state, input copies and foreign
stress/tangent buffers are checked for invariance during all evaluations.

The unchanged tangent bound is `2e-6*||C_alg||_inf + 1e-8*||C_el||_inf`.
Recorded relative matrix and directional errors use `||C_alg||_inf` as denominator.
The records report maxima over ALL accepted steps, not the best single error.
Stress comparisons retain `1e-7 MPa + 1e-9` of reference magnitude; history uses
`1e-12` absolute. No acceptance tolerance was relaxed for this correction.

Additional lifecycle checks query committed yield after explicit import and after
binary/Trace restart with an invalid cache, against the independent Swift formula.
Queries leave history and foreign buffers unchanged. All four response-flag
combinations are checked both with missing unrequested buffers and with attached
unrequested buffers. Tangent-only output equals the full response from the same
committed state. Existing failed-finalization, clone and continuation checks remain.

Validation after the archive correction (2026-10-24, Release, GCC 13.3, C++20,
one OpenMP thread):

| Check | Result |
| --- | --- |
| Focused Swift C++ tests, including archive/tangent/lifecycle checks | 17/17 passed |
| Full ConstitutiveLaws C++ suite (includes those 17) | 85/85 passed |
| Both existing plane-stress serialization regressions from #14895 | Passed within the 85 |
| Python registration and unchanged material-definition tests | 2/2 passed |
| ConstitutiveLaws Python nightly (includes those 2) | 30/30 passed |
| Tangent case/step coverage | 36 cases, 292/292 branch-faithful steps passed |
| Minimum adjacent accepted steps per case | 5 |
| Maximum relative matrix error over accepted steps | 1.622901e-7 |
| Maximum relative directional error over accepted steps | 8.345052e-8 |
| `git diff --check` and legacy source audit | Passed |

The small roundoff/residual plateau near yield and the larger truncation error at
coarser h remain below the original acceptance bound; no tolerance change was
needed. Per-case maxima and counts are retained in the GTest XML and validation
summary. Both review findings are addressed by the early header correction and
the permanent case table; these tests do not depend on temporary review files.

Build with `bash build/configure.sh KratosConstitutiveLawsCoreTest KratosConstitutiveLawsApplication`,
then `python3 build/stage_runtime.py`. Run `bash build/run_checks.sh focused`,
`bash build/run_checks.sh cpp_all`, and `bash build/run_checks.sh python`.
The existing helpers stage libraries only from this worktree; `runtime_manifest.json`
and `runtime_provenance.json` identify actual source/destination hashes and module
paths. Build log: `build/logs/archive_fix_build.log`; test logs:
`archive_fix_focused.log`, `archive_fix_cpp_all.log`, `archive_fix_python.log`.
The XML files and `build/logs/validation_summary.json` describe this correction run.

No broad Core/StructuralMechanics suites, FE benchmark or performance harness was
run in that correction. Voce and further physics remain outside this package.
The separate small FE lifecycle run below was added afterwards.
The architecture, selection value and public archive contract remain proposals
subject to maintainer review; no ABI or universal zero-cost claim is made.


## Homogeneous 3D FE lifecycle and model restart (2026-10-25)

`test_small_strain_swift_fe.py` is included in the small/nightly/all suites. It
uses one registered `SmallDisplacementElement3D8N`, a straight 1 mm cube, and the
existing `Materials.json` through `ReadMaterialsUtility`. Units are mm, N, MPa.
All eight Gauss points are checked. Element registration is in StructuralMechanics
`structural_mechanics_application.cpp:650`; `SmallDisplacement::UseElementProvidedStrain`
returns true. The material definition adds no fracture energy or separate initial
yield stress.

Symmetry conditions prescribe ux=0 on x=0, uy=0 on y=0, uz=0 on z=0; ux is driven
on x=1. The eight remaining transverse DOFs are solved by the standard
`ResidualBasedNewtonRaphsonStrategy`, incremental static scheme, block builder,
and Skyline LU. Mesh movement is disabled. The solver rebuilds the material
stiffness during Newton. Only a converged SolveSolutionStep may be finalized.
There is no prescribed analytical transverse displacement field.

The full reference path is computed before the FE solve. Its independent scalar
bisection uses `abs(E*(eps-ep_axial_old)) - E*dp - K*(eps0+p_old+dp)^n = 0`.
The modulus here is E because the stress state is uniaxial; the reference does not
use the 3G equation of the 3D elastic predictor. Accumulated p and signed axial
plastic strain are distinct stored reference quantities, updated with dp and
sign(trial)*dp respectively. Lateral strain is `-nu*sigma/E-ep_axial/2`; the plastic
engineering vector is `[ep_axial,-ep_axial/2,-ep_axial/2,0,0,0]`.

Positive tension gives a positive summed KM.REACTION_X on x=1 and a negative sum
on x=0. Both signed sums are compared with +/-sigma times the ORIGINAL 1 mm^2
area; their sum and all free-DOF reactions are checked. Independently assembling
the element RHS checks the L2 equilibrium residual of the eight free DOFs.

| Step | Axial strain | Expected phase | Global Newton iterations, A/B |
| --- | --- | --- | --- |
| 1 | 0.0005 | elastic | 1/— |
| 2 | 0.0015 | elastic | 1/— |
| 3 | 0.0020 | first plastic flow | 3/— |
| 4 | 0.0030 | plastic loading | 2/— |
| 5 | 0.0045 | plastic loading | 2/— |
| 6 | 0.0060 | plastic loading / checkpoint | 2/— |
| 7 | 0.0050 | elastic unloading | 1/1 |
| 8 | 0.0045 | elastic unloading | 1/1 |
| 9 | 0.0055 | elastic reloading | 1/1 |
| 10 | 0.0060 | return to previous yield | 1/1 |
| 11 | 0.0065 | renewed plastic flow | 2/2 |
| 12 | 0.0080 | further plastic flow | 2/2 |

The JSON distinguishes the converged trial response before FinalizeSolutionStep
from committed history afterwards: `trial.stress/strain` describe the converged
configuration, while `trial.p/plastic_strain` are still the committed inputs.
These inputs remain unchanged through Newton and output queries. After normal
strategy finalization all IP histories match the independent reference. During
unloading/reloading p remains exactly unchanged; steps 11 and 12 increase it.
Stress, strain and plastic state are checked for homogeneity at every IP.

The binary restart uses the existing `FileSerializer` on the entire `Model` after
step 6, at nonzero p approximately 0.00397306854. It includes nodes, DOFs/fixities,
two solution-step buffers, properties, ProcessInfo, element data and all eight IP
laws. The file lives in a fresh temporary directory for this test only. A new
Model is loaded and its signature, including previous nodal buffers and nonzero
history, is compared immediately. IS_RESTARTED and LOAD_RESTART are then set as
in the core RestartUtility, and a NEW solver is initialized normally.
`BaseSolidElement::Initialize` respects IS_RESTARTED; the second signature check
verifies that loaded material history was not reset. No history is copied back
manually. Steps 7-12 then compare both continuations, including trial/committed
stress, strain, history, displacement and reaction. Iteration-count identity is
not required, and no convergence-order claim is inferred from these few iterations.

Acceptance limits were fixed before FE evaluation: absolute 1e-6 MPa for stress,
1e-6 N for reaction, 1e-11 for strain/history and 1e-11 mm for displacement, plus
1e-9 times the reference magnitude. Local stress resolution is about 1.5e-8 MPa
for this C, so the stress/force allowances also cover global equilibrium and
roundoff. Strain/history allowances cover stress-to-strain conversion and rounding.
The normal ResidualCriteria uses relative 1e-10 and absolute 1e-9; the independent
free-DOF L2 residual must additionally be <=1e-7 N. No tolerance was relaxed.

Newly executed checks: standalone FE test 1/1, existing registration/material tests
2/2, ConstitutiveLaws Python nightly 31/31 (includes these three), and diff --check.
The 17/85 C++ results above are earlier evidence, not rerun for this Python-only
addition. No native sources changed and no native rebuild was needed. Runtime
paths and SHA256 hashes match the existing build/staged libraries of this worktree.
Their build metadata records 8b7484c2ae7, when the prototype was uncommitted;
the user subsequently committed that prototype as f282677a864 before this FE task.
The provenance JSON records both revisions explicitly.

| Maximum over all accepted A/B steps | Observed absolute error |
| --- | ---: |
| Axial stress [MPa] | 4.213291e-09 |
| Transverse/shear stress [MPa] | 6.531309e-11 |
| IP stress spread [MPa] | 8.810730e-13 |
| Accumulated p | 2.006208e-14 |
| Plastic strain, any component | 2.006208e-14 |
| Total strain, any component | 4.013283e-15 |
| Displacement [mm] | 4.013283e-15 |
| Signed axial face reaction [N] | 4.212893e-09 |
| Free-DOF L2 residual [N] | 4.611346e-11 |

All 12 A steps and all 6 restarted B steps converged, with 1-3 global Newton
iterations. Maximum A/B continuation differences in stress, strain, p, plastic
strain, displacement and reaction were all 0.0 in this run. No production-law,
integrator, solver or archive correction was required by the FE test.

Separate reproduction from the worktree root:

```bash
export PYTHONPATH="$PWD/bin/Release"
export LD_LIBRARY_PATH="$PWD/bin/Release/libs"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
KRATOS_SWIFT_FE_REPORT="$PWD/build/logs/swift_fe_validation.json" \
  /usr/bin/python3 applications/ConstitutiveLawsApplication/tests/test_small_strain_swift_fe.py
```

The optional report is output only. The test never reads previous build/logs or
files from earlier temporary runs. `swift_fe_validation.json` contains the fixed
reference path, all IP/node outputs, signed reactions, global convergence records,
checkpoint and error maxima. With the optional report enabled, the standard
Newton strategy/criterion also prints its per-iteration residuals to the test log. `swift_fe_runtime_provenance.json` verifies native
paths/hashes. `validation_summary.json` keeps previous results and this separate
FE run. Test logs: `swift_fe_test.log`, `swift_fe_registration.log`,
`swift_fe_python_suite.log` under build/logs.

This homogeneous one-element case establishes the intended FE lifecycle and
uniaxial response. It provides no general mesh, element-family or structural
validation. Voce, mesh studies and performance/benchmark campaigns were not added.
