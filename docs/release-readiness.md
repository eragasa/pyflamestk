# Release-readiness audit

This audit is bounded to static inspection and an isolated wheel-build attempt
of PyFlamestk commit `5b8368cc88d91bc56f9cc1c8a7fa8d9ea3d6b359` on branch
`dev/release-readiness`. It does not import PyFlamestk, run tests, execute
calculator scripts, submit jobs, evaluate source expressions, or establish
numerical or scientific validity.

## Definition of ready

A source baseline is release-ready when it has an explicit package boundary,
reproducible metadata and dependency declarations, a buildable wheel, a bounded
non-calculator test suite, CI enforcement, documented limitations, and no
restricted or generated calculator artifacts in the maintained tip or release
artifacts. An immutable annotated tag should identify the accepted baseline;
a moving branch is not provenance.

## Evidence summary

- Git branch before audit: `master`.
- Exact source commit: `5b8368cc88d91bc56f9cc1c8a7fa8d9ea3d6b359`.
- Exact Git tree: `02e20f61a9b554ed0dcbda15bb24adc21e942007`.
- Tracked inventory before cleanup: 955 files and 290,287,120 bytes.
- Directory sizes before cleanup: `dev/` 233,625,305 bytes, `examples/`
  46,149,972 bytes, `tests/` 10,153,017 bytes, and `pyflamestk/` 352,234
  bytes.
- Python inventory: 99 tracked Python files and 42 test-like paths.
- CI configuration: absent.
- `.gitignore` at the audited commit: absent.
- Isolated wheel build: failed during package discovery because both `dev` and
  `pyflamestk` were interpreted as top-level packages.

## Findings

### Package boundary is not buildable

Disposition: `MUST_FIX`

`setup.py` declares project metadata but does not declare packages,
dependencies, supported Python versions, package data, entry points, or a build
backend. Current setuptools aborts with:

```text
Multiple top-level packages discovered in a flat-layout: ['dev', 'pyflamestk'].
```

The smallest correction is explicit modern build metadata that includes only
the intended `pyflamestk` package. Examples, development material, tests,
offline artifacts, and calculator output must not enter the wheel.

### Maintained-package syntax is incomplete

Disposition: `MUST_FIX`

Static AST parsing found an indentation error in
`pyflamestk/dakota_interface.py:187`. Three additional source-bearing example
or development scripts also fail parsing:

- `dev/2016_MRS_spring_pareto/buckingham_pareto_clean.py:348`
- `dev/2016_MRS_spring_pareto/buckingham_pareto_iterate.py:349`
- `examples/Ni_eam/eam_potential.py:17`

The package module must be repaired or explicitly excluded from the supported
package boundary. Invalid development and example scripts must be classified as
unsupported evidence or corrected in their owning scope; they must not silently
pass a release gate.

### Runtime dependencies are undeclared

Disposition: `MUST_FIX`

Static import inspection found external roots including NumPy, SciPy,
Matplotlib, pandas, PyYAML, mpi4py, seaborn, and scikit-learn. `setup.py`
declares none. The supported package surface must be selected before separating
required dependencies from optional calculator, plotting, MPI, and development
dependencies.

### Test and CI contract is absent

Disposition: `MUST_FIX`

The repository has no CI workflow or modern test configuration. Existing tests
include calculator outputs, pseudopotentials, scheduler output, and scripts that
can launch external programs. Release verification must begin with static and
pure in-memory tests. Calculator and scheduler execution must remain explicit,
opt-in, and outside ordinary CI.

### Generated and restricted artifacts were tracked

Disposition: `MUST_FIX` — addressed in the current working tree but not yet
committed.

The audited commit tracks VASP pseudopotentials and generated outputs, LAMMPS
logs and restarts, large optimization datasets, scheduler output, editor state,
and interpreter caches. The release-readiness working tree moved 161 files
(250,836,932 bytes) to ignored `.offline/release-readiness/` staging.
`OFFLINE_ARTIFACT_SHA256SUMS` preserves their original path identities, and
`OFFLINE_ARTIFACTS.md` defines the boundary. Every staged byte was verified
against the original `HEAD` content.

Existing Git history is unchanged. No release should redistribute a
pseudopotential without separately established authority.

### Machine-specific examples and effectful APIs remain

Disposition: `SAFE_TO_DEFER`

Examples contain absolute paths, executable selections, scheduler commands, and
calculator launch scripts. Package modules contain subprocess boundaries and
historical `eval()` use. These prevent examples from being advertised as
portable supported workflows, but they need not block a source baseline if:

- examples and development trees are excluded from distributions;
- ordinary tests never execute those paths;
- the README clearly states the effect and portability boundary; and
- supported APIs are narrowed before publication.

Revisit this finding before any calculator-backed feature is presented as
supported.

## Recommended sequence

1. Commit the reviewed offline-artifact boundary without rewriting history.
2. Add explicit package metadata and build only the intended package.
3. Define the supported module surface and resolve or exclude
   `dakota_interface.py`.
4. Declare core and optional dependencies from that supported surface.
5. Add bounded static and pure in-memory tests with sanitized fixtures.
6. Add CI for tests, package-boundary checks, syntax, and wheel isolation.
7. Document execution, portability, numerical-verification, and
   scientific-validation limits.
8. Build and inspect the wheel in isolation.
9. Select a release version and create an annotated tag only after all blocking
   findings are resolved.

## Review result

Review outcome: `CHANGES_REQUIRED`

The external-artifact boundary is now staged and verified. Packaging, syntax,
dependency, test, CI, and support-surface work remains before a release tag is
appropriate.
