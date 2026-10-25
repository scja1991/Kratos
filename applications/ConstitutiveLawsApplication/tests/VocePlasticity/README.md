# Local small-strain Voce extension

Base: completed Swift commit `5b13b5748a1d2df7167d0937bccf7d8ff3d9c900`
(`Validate small-strain Swift plasticity in a 3D FE model`).
Branch `feature/small-strain-voce-hardening`, isolated worktree
`/home/localadmin/Kratos-small-strain-voce-hardening`.
The common architecture is a local proposal, without maintainer approval.

Use registered `SmallStrainJ2IsotropicHardening3DLaw` with `HARDENING_CURVE=8`.
No Voce wrapper or new material variable is introduced. `Materials.json` is a
minimal working definition; its existing variables mean:

| Variable | Meaning | Illustrative value |
| --- | --- | ---: |
| YOUNG_MODULUS | E, MPa | 210000 |
| POISSON_RATIO | nu | 0.3 |
| YIELD_STRESS | initial yield sigma_y0, MPa | 250 |
| EXPONENTIAL_SATURATION_YIELD_STRESS | saturation sigma_infinity, MPa | 450 |
| HARDENING_EXPONENT | b, inverse plastic strain | 20 |
| HARDENING_CURVE | explicit Voce selection | 8 |

The pure Voce function is `sigma_y(p)=sigma_y0+Q*(1-exp(-b*p))`,
`H(p)=Q*b*exp(-b*p)`, with `Q=sigma_infinity-sigma_y0`.
Strains/p are dimensionless, so H uses the same stress units as sigma_y.
Finite sigma_y0>0, sigma_infinity>=sigma_y0, b>=0 and p>=0 are required.
Q=0 and b=0 give perfect plasticity. Optional ISOTROPIC_HARDENING_MODULUS=0
is accepted; any nonzero/nonfinite value is rejected. This is not a Voce-plus-linear
model. Optional tension/compression yields must equal the initial yield.

Small bp uses expm1, with reordered products if bp itself underflows but its
strength contribution is representable. Large bp approaches the saturation
stress. The derivative fallback uses frexp/scalbn and normalized products with exp(-512)
chunks for overflowing Q*b or subnormal/underflowing exponential factors. It does
not reconstruct a large derivative from a rounded sum of logarithms. Mantissas
stay normal; the final binary scaling detects overflow without clamping or a
rejection margin. For bp>4096 the derivative rounds to zero: Q,b<2^1024 implies
Q*b*exp(-4096)<2^-3860, far below half the smallest subnormal double.
A nonrepresentable result returns false without changing either output. Check
also requires a representable initial derivative. Derivative underflow to zero
in saturation is allowed. Iteration callbacks hold no material state.
This is not a claim of correctly rounded results for every double input.

Only MakeEvolutionContext, ResolveHardening and model-specific checks are extended.
Voce uses the existing history callback, potential/Hessian, coupled endpoint solve,
analytic discrete tangent, trial/finalization, import, copy/clone and restart.
The generic integrator changes only by appending enum value 8; its algorithms are
unchanged. Swift still selects 7 implicitly in its wrapper and needs no Voce fields.
Existing enum values and dissipation-based paths retain their meanings.

The common archive remains version 2, with its existing early header check and
payload. No added history, archive field or migration is needed. Restart requires
the matching model parameters and explicit hardening selection in Properties.
Historical Swift layouts retain their previous rejection contract.

Scope remains isothermal, rate-independent associated 3D infinitesimal J2,
element-provided engineering strain and analytic tangent. Active thermal tables,
accessors, thermal expansion, U-P and F-input are rejected. Finite strain, damage,
kinematic/combined hardening, reduced dimensions and substepping are outside scope.

The C++ tests use independent long-double function/derivative references and scalar
bisection of `q_trial-3G*lambda-sigma_y(p_old+lambda)=0`; they do not use the separate
SmallStrainJ2Plasticity3D as a reference. They cover parameter/domain failures,
small/large/overflowing bp, Q=0/b=0, signed loading, shear, unloading/re-yielding,
nonproportional history, trials, failed/successful commit, reset, import, concrete
family clone and binary/precise Trace restart continuation with invalid cache.

A small shared tangent test helper preserves the existing Swift stencils/tolerances.
Voce tests nine paths for nonlinear, Q=0 and b=0 configurations plus a saturating
state at bp=10. Both matrix and independent direction errors must pass the same
steps, with at least two adjacent accepted h. Branch fidelity comes from an
independent elastic predictor and yield formula. The unchanged bound is
`2e-6*||C_alg||_inf + 1e-8*||C_el||_inf`. No two-sided derivative is demanded at
the yield switch. Per-case maxima over accepted steps are written into GTest XML.

The Voce FE subclass reuses the existing Swift cube/Newton/model-restart lifecycle.
Only material definition, model metadata and independent scalar reference change.
Its uniaxial bisection uses E, not 3G. Accumulated p and signed axial plastic strain
remain separate. The reference is computed before FE; Voce starts yielding from
250/E, without a hard-coded Swift threshold. The same 12-step path up to 0.008,
free transverse DOFs, eight Gauss points, signed reactions and full binary Model
restart after step 6 are retained. Both Swift and Voce keep separate automated tests.
The shared FE acceptance limits are unchanged: stress/reaction 1e-6 absolute,
state/displacement 1e-11 absolute, plus relative 1e-9; free-DOF residual L2<=1e-7 N.

Build through `bash build/configure.sh`, then `python3 build/stage_runtime.py`.
Only local dependency SOURCES were copied from the Swift workspace; all native
libraries are built in this worktree. Use its bin/Release and libs paths at runtime.
Focused and suite commands and measured results are recorded after validation in
build/logs/validation_summary.json; raw XML/logs/FE comparison JSON stay under build/logs.
No benchmarks or general mesh/structure validation are claimed for the homogeneous
one-element proof. Properties are not supplemented to satisfy legacy requirements.


## Original validation before review corrections (2026-10-25)

| Check | Result |
| --- | ---: |
| Focused Voce C++ (function, extremes, selection, paths, lifecycle/restart, tangent) | 6/6 passed |
| Existing Swift material-point/tangent/archive C++ | 17/17 passed |
| Full ConstitutiveLaws C++ (includes both groups) | 91/91 passed |
| Existing small-/finite-strain plane-stress serialization from #14895 | Both passed within the 91 |
| Voce selection/import Python | 2/2 passed |
| Swift registration/material Python | 2/2 passed |
| Separate Swift and Voce FE/model-restart tests | 1/1 each passed |
| Relevant ConstitutiveLaws Python nightly (includes the six tests above) | 34/34 passed |
| Diff whitespace/source reuse checks | Passed |

Voce tangent table: 28 cases, all 228 branch-faithful steps accepted. Matrix and
independent direction passed the same h; at least five adjacent steps passed per
case. Maximum relative matrix error 1.394537e-7; maximum normalized directional
error 8.297511e-8. The denominator is ||C_alg||_inf, and maxima cover ALL accepted
steps. Near initial yield h=1e-7 down to 1e-9 is branch-faithful; the other cases
accept all nine h from 1e-5 to 1e-9. Swift retains 36 cases / 292 accepted steps,
maximum matrix error 1.622901e-7 and directional error 8.345052e-8. No tolerance
was relaxed. The shared sweep body/stencils were checked byte-identical after
substituting only the point type and scalar reference function.

Extreme derivative checks additionally use rounded constants calculated with
80-digit Decimal arithmetic, avoiding dependence on extended long-double range.
The initial standalone Python selection invocation revealed a missing main entry
point in the new test file; it was added, and the two tests were then actually
executed and passed. No solver, history, tangent, lifecycle or archive correction
was needed after running the tests.

Both FE cases converge for all 12 uninterrupted and all 6 restarted continuation
steps. Voce uses 1-2 global Newton iterations; Swift uses 1-3. The independently
computed reference yields first at step 2 for Voce versus step 3 for Swift. These
counts are not evidence of a convergence order. Measured maxima:

| Quantity | Voce | Swift |
| --- | ---: | ---: |
| Axial stress error [MPa] | 4.184756e-09 | 4.213291e-09 |
| Transverse/shear stress [MPa] | 6.110668e-13 | 6.531309e-11 |
| IP stress spread [MPa] | 7.958079e-13 | 8.810730e-13 |
| p error | 1.992720e-14 | 2.006208e-14 |
| Plastic strain error | 1.992677e-14 | 2.006208e-14 |
| Total strain error | 3.986395e-15 | 4.013283e-15 |
| Displacement error [mm] | 3.986828e-15 | 4.013283e-15 |
| Signed axial face reaction error [N] | 4.184528e-09 | 4.212893e-09 |
| Free-DOF L2 residual [N] | 2.392213e-13 | 4.611346e-11 |

Maximum A/B restart differences in stress, total/plastic strain, p, displacement
and reaction are 0.0 for both laws in this run. The complete model, including
Properties selection 8 and every IP law, is restored with a new solver; no history
is restored manually. All existing Swift bounds and load steps remain unchanged.

Reproduce from this worktree root:

```bash
bash build/configure.sh
python3 build/stage_runtime.py
bash build/run_checks.sh voce
bash build/run_checks.sh focused
bash build/run_checks.sh cpp_all
bash build/run_checks.sh fe
bash build/run_checks.sh python
```

For just Voce Python, set PYTHONPATH=$PWD/bin/Release,
LD_LIBRARY_PATH=$PWD/bin/Release/libs, PYTHONDONTWRITEBYTECODE=1 and
OMP_NUM_THREADS=1, then run
`/usr/bin/python3 applications/ConstitutiveLawsApplication/tests/test_small_strain_voce_hardening.py`
or `test_small_strain_voce_fe.py` in the same folder. The optional
KRATOS_VOCE_FE_REPORT output path writes FE data; the automated tests read no old
build/logs files and depend on no previous temporary files.

All native Core, StructuralMechanics, LinearSolvers and ConstitutiveLaws libraries
were freshly built here; no Swift-only library was reused. The staged runtime
matches its build-file SHA256 hashes. `runtime_provenance.json`,
`native_hash_provenance.json` and `cpp_runtime_provenance.log` verify actual paths;
`validation_summary.json`, GTest XML and both `*_fe_validation.json` files contain
results. Build logs are `bootstrap_build.log` and `voce_build.log`.

The reuse audit verifies unchanged integrator algorithms and all common methods
from state validation through serialization. Only active-model checks/context
selection and enum value 8 changed in production. One existing unsupported-choice
Swift test now uses 9, because 8 is supported here. The prior Props assignment
behavior retains accessors; new accessor tests therefore use fresh fixtures,
without changing that unrelated Core behavior.

Limits: no Finite Strain, damage, kinematic/combined hardening, reduced dimensions,
substepping or wider mesh/structure/performance claim. The homogeneous FE case
checks the intended lifecycle and uniaxial response. Maintainer acceptance of the
local common architecture remains open.


## Review corrections: representability and header independence

The previous logarithmic derivative fallback could return true just beyond the
upper double range. For binary inputs sigma_y0=1, sigma_infinity=1e308, b=10 and
p=0.17160808417760906, the independent reference is
H=1.797693134862367906...e308, which rounds to infinity. The corrected function
returns false and preserves both output arguments. This concerns an individual
scalar evaluation. The material was already rejected by Check at p=0; this
finding does not invalidate the existing FE proof.

Only the derivative fallback and the missing direct includes/variables.h include
changed in production. Its scaled mantissa/exponent calculation has no DBL_MAX
clamp and no acceptance/rejection safety margin. The direct normal product path,
small-bp strength calculation, ideal-plastic limits and permitted saturation
underflow are retained. The common family's parameter checks and all solver,
history, potential, tangent, lifecycle and archive methods remain unchanged.

VoceEndpointRepresentabilityTransitions freezes twelve cases as hexadecimal
literals: the review failure, finite/overflowing H close to DBL_MAX, both sides of
bp=DBL_MIN and bp=log(2), both sides of exp(-bp) becoming subnormal, and normal,
subnormal and zero H deep in saturation. Each side is checked against its own
reference, so a jump inconsistent with the smooth mathematical function is
bounded by the reference tolerances. No assertion depends on a private branch
number or the particular exponential chunk size.

References were generated with Python Decimal at precision 400, using
Decimal.from_float for EVERY input, q=saturation-initial at that precision and
H=q*b*exp(-b*p). Tiny 1-exp(-bp) was evaluated through its convergent series.
The test needs neither Python/Decimal nor extra long-double exponent range:
only the independently rounded constants are compiled. Bounds are 8*epsilon
relative for yield, 32*epsilon for normal H, and an absolute floor of two smallest
subnormal quanta. These cover ordinary exp/product rounding, including several
normalized products; the success flag is checked exactly, separately from value
tolerances. Existing function, tangent and FE tolerances were not relaxed.
There is no claim of correctly rounded evaluation for every double input.

VoceEndpointParameterIsolation checks that inactive Swift fields (including
invalid scalars and thermal metadata) do not affect Voce, and inactive Voce
saturation/b fields do not affect Swift. Stress and tangent remain identical.
YIELD_STRESS is NOT inactive for Swift: optional yield/tension/compression values
retain their redundancy contract. Accessors and temperature tables are rejected
on those redundant values and on Voce's optional zero linear modulus. The latter
has no corresponding active Swift contract.

The pre-fix regression was built and run before editing the production header;
its review row fails. No git reset or restoration of repository files was used.
Evidence is under build/logs/voce_review_fixes/numeric_before.{log,xml}; an exact
copy of the prior header is retained there for diagnosis only.

The isolated translation unit contains exactly:

```cpp
#include "custom_constitutive/auxiliary_files/hardening_laws/voce_hardening_law.h"
```

It is compiled by /usr/bin/g++ -fsyntax-only with CXX_DEFINES, CXX_INCLUDES and
CXX_FLAGS read from the existing KratosConstitutiveLawsCore.dir/flags.make.
It uses no unity source, prior project header or precompiled-header injection.
The local reproduction command is
`python3 build/logs/voce_review_fixes/header_check.py after`;
header_after.json records the complete expanded compiler command. Before the fix,
the same check failed on YIELD_STRESS; afterwards it exits 0.


### Validation after the corrections (2026-10-25)

| Newly executed check | Result |
| --- | ---: |
| New boundary and parameter-isolation tests | 2/2 passed |
| All focused Voce C++ | 8/8 passed |
| All focused Swift C++ | 17/17 passed |
| Full ConstitutiveLaws C++ including Swift archives and both #14895 serialization cases | 93/93 passed |
| Relevant Swift/Voce Python material/lifecycle and separate FE-restart cases | 6/6 passed |
| Isolated header syntax check, outside unity | Passed (pre-fix failure also retained) |
| git diff --check | Passed |

For the new frozen table, the largest relative yield error is 0.0 against its
rounded constants, and the largest relative normal-H error is 2.728533e-16.
Subnormal H is checked with the stated absolute quantum bound. This is evidence
for these inputs, not a guarantee of exact rounding for all double inputs.
The two new tests are VoceEndpointRepresentabilityTransitions and
VoceEndpointParameterIsolation. Their independent constants require no files
under /tmp and no high-precision dependency during test execution.

The correction build used
`bash build/configure.sh KratosConstitutiveLawsCoreTest KratosConstitutiveLawsApplication`,
followed by `python3 build/stage_runtime.py`. Native build/staged hashes match for
all nine recorded libraries, and actual module/ldd paths resolve to this worktree.
Raw current logs, XML, reference data and FE JSON are in
build/logs/voce_review_fixes; build/logs/validation_summary.json is the current
summary. Its previous version was preserved as previous_validation_summary.json
in that subdirectory. The full Python nightly above is historical; this task
reran only the six relevant Python tests. No broad Core/StructuralMechanics suite,
mesh/structure study or performance campaign was run.

Both FE cases again converge for every 12-step path and six-step restart
continuation. Voce maxima remain: axial stress 4.184756e-9 MPa, p 1.992720e-14,
displacement 3.986828e-15 mm, signed axial reaction 4.184528e-9 N. Swift axial
stress error is 4.213291e-9 MPa. Maximum restart differences remain 0.0 in all
recorded fields for both laws. The correction required no new FE scope and no
change to the FE reference or acceptance bounds.

To rerun the new C++ cases separately with this worktree's runtime libraries:

```bash
export LD_LIBRARY_PATH="$PWD/bin/Release/libs"
export OMP_NUM_THREADS=1
build/Release/applications/ConstitutiveLawsApplication/KratosConstitutiveLawsCoreTest --gtest_filter='*VoceEndpointRepresentabilityTransitions:*VoceEndpointParameterIsolation'
```

The earlier focused/full-suite/FE commands remain applicable. The relevant Python
selection is `python3 -m unittest -v test_small_strain_j2_hardening
test_small_strain_voce_hardening test_small_strain_swift_fe
test_small_strain_voce_fe` from the tests folder, with this worktree's PYTHONPATH,
LD_LIBRARY_PATH and PYTHONDONTWRITEBYTECODE=1.
