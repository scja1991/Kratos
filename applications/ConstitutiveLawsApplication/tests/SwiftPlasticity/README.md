# Local small-strain Swift prototype

Based on public upstream master `8b7484c2ae79451b8c83e6831865f6641dc62018`.
This is a local architecture proposal, without maintainer approval.

`SmallStrainSwiftJ2Plasticity3DLaw` accepts the unchanged `Materials.json`:
`YOUNG_MODULUS`, `POISSON_RATIO`, `SWIFT_COEFFICIENT`,
`SWIFT_INITIAL_STRAIN`, and `SWIFT_HARDENING_EXPONENT`.
`HARDENING_CURVE` is optional for this wrapper; if supplied it must be 7.
The shared `SmallStrainJ2IsotropicHardening3DLaw` requires that explicit selection.
Only Swift is implemented. Other selections are rejected.

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
run in this correction. Voce, further physics and FE implementation remain outside
this package. A small 3D FE demonstration with the new Swift law is the next step.
The architecture, selection value and public archive contract remain proposals
subject to maintainer review; no ABI or universal zero-cost claim is made.
