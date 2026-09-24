# Release-readiness audit

The initial audit was bounded to static inspection and an isolated wheel-build
attempt of PyFlamestk commit
`5b8368cc88d91bc56f9cc1c8a7fa8d9ea3d6b359`. Remediation validation on
`dev/release-readiness` additionally builds distributions, imports packaged
modules with declared dependencies, and runs sanitized infrastructure tests. It
does not execute calculator scripts, submit jobs, evaluate source expressions,
or establish numerical or scientific validity.

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

## Remediation evidence

- Artifact-boundary commit: `64e7b50`.
- Explicit packaging and CI commit: `ed0f248`.
- Tracked inventory after cleanup: 808 files and 39,504,965 bytes.
- Removed collection: 161 files and 250,836,932 bytes, verified against the
  original commit and preserved in a private manifest-bound archive.
- Package syntax: all 17 packaged Python modules parse and compile without a
  `SyntaxWarning`.
- Distribution boundary: 22 wheel files and 40 source-distribution files;
  neither artifact contains `dev/`, `examples/`, offline bytes, or calculator
  output.
- Isolated imports: the base wheel imports without dependencies, and all 16
  packaged submodules import with declared core and optional dependencies.
- Bounded tests: six tests cover package selection, safe YAML loading, wheel
  boundaries, and the PyFlamestk-specific offline staging policy.
- CI: Python 3.11 and 3.14 passed tests, compilation, style checks, distribution
  builds, boundary verification, dependency installation, and isolated module
  imports in [run 35957687040](https://github.com/eragasa/pyflamestk/actions/runs/35957687040).

## Findings

### Package boundary is explicit and buildable

Disposition: `NO_ACTION_REQUIRED`

`pyproject.toml` selects only `pyflamestk`, defines Python 3.11 or newer, and
uses a metadata-only `setup.py` compatibility shim. `MANIFEST.in`, wheel and
source-distribution verifiers, and CI enforce the source and binary boundaries.
Both distributions build in isolated environments. Examples, development
material, offline bytes, and calculator output are absent from release
artifacts.

### Maintained-package syntax is complete

Disposition: `NO_ACTION_REQUIRED`

The indentation defect in `pyflamestk/dakota_interface.py` is repaired. Unsafe
unqualified YAML loading was replaced with `safe_load` and `safe_dump`, and
invalid regular-expression escapes were corrected. Every packaged module now
parses, compiles, and imports in CI.

Three development or example scripts still fail static parsing. They are
excluded from both distributions and classified with the unsupported examples
below; no release gate silently treats them as maintained package code.

### Runtime dependencies are declared by supported surface

Disposition: `NO_ACTION_REQUIRED`

NumPy and SciPy are core requirements. Matplotlib and PyYAML are explicit
`plot` and `dakota` extras and are combined by the `all` extra. Imports that
occur only under excluded development and example trees do not become wheel
requirements. CI installs all declared extras and imports every packaged
module.

### Test and CI contract is bounded

Disposition: `NO_ACTION_REQUIRED`

Ordinary tests use only synthetic temporary repositories, manifests, and
archives. CI compiles maintained source, runs the bounded tests, validates
wheel and source-distribution contents, installs declared dependencies, and
imports modules outside the checkout on Python 3.11 and 3.14. It never invokes
a calculator or scheduler.

### Generated and restricted artifacts are outside the maintained tip

Disposition: `NO_ACTION_REQUIRED`

Commit `64e7b50` removes VASP pseudopotentials and generated outputs, LAMMPS
logs and restarts, large optimization datasets, scheduler output, editor state,
and interpreter caches from the maintained tip. The 161 original files remain
in ignored local staging and a private offline archive.
`OFFLINE_ARTIFACT_SHA256SUMS` preserves original paths and identities, while
`OFFLINE_ARTIFACTS.md` records the archive identity and boundary. Existing Git
history is unchanged. No release redistributes the removed pseudopotential.

### Machine-specific examples and effectful APIs remain

Disposition: `SAFE_TO_DEFER`

Examples contain absolute paths, executable selections, scheduler commands, and
calculator launch scripts. Package modules contain subprocess boundaries and
historical `eval()` use. These prevent examples from being advertised as
portable supported workflows, but they need not block a source baseline if:

- examples and development trees are excluded from distributions;
- ordinary tests never execute those paths;
- the README clearly states the effect and portability boundary; and
- this source baseline is described as alpha research software without a
  calculator, numerical-verification, or scientific-validation claim.

Three excluded scripts remain syntactically invalid. Revisit this finding and
repair or retire the relevant material before any example or calculator-backed
feature is presented as supported.

## Recommended sequence

The artifact boundary, packaging metadata, dependency declarations, bounded
tests, CI, limitation documentation, and isolated distribution checks are
complete. The remaining release operation is to merge or fast-forward the
validated branch and create an annotated `v0.1.0` tag at the accepted commit.
Tagging is a publication action, not a technical-review step.

## Review result

Review outcome: `NO_BLOCKING_FINDINGS`

The bounded source-baseline definition is met. This result does not claim that
historical examples work, that calculators are available, or that any numerical
or scientific result is verified or validated.
