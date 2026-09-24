# Offline artifact boundary

The `dev/release-readiness` branch removes generated calculator output,
restricted calculator input, large optimization result datasets, editor state,
and interpreter caches from the maintained repository tip.

The original bytes are staged locally under the ignored directory:

```text
.offline/release-readiness/
```

`OFFLINE_ARTIFACT_SHA256SUMS` records each original repository-relative path and
its SHA-256 identity. The ignored directory contains `MANIFEST.tsv` with the
same digest, byte size, classification reason, and original path. The staging
directory is operator-managed material intended for offline storage; it is not a
repository artifact, package input, test fixture, or supported runtime path.

## Selection boundary

The initial move contains 161 files totaling 250,836,932 bytes. It includes:

- VASP pseudopotential and generated calculator files such as `POTCAR`,
  `WAVECAR`, `CHGCAR`, `OUTCAR`, and `vasprun.xml`;
- LAMMPS logs, restart files, and generated output;
- generated optimization parameter and result datasets;
- scheduler output; and
- Spyder, AppleDouble, and Python cache files.

Authored Python, shell scripts, configuration, structures, compact calculator
inputs, and potential definitions remain in the repository unless separately
classified. A filename extension alone was not used to remove every `.dat` or
`.out` file because those extensions may also contain authored input.

## Preserved procedure

The PyFlamestk-specific selection and safe-copy procedure is versioned under
[`tools/`](tools/README.md). `tools/stage_offline_artifacts.py` must not be
rerun against the current staging destination.

Reusable verification and deterministic tar creation are owned by the
[Project Koios bootstrap offline-artifact candidate](https://github.com/eragasa/projectkoios-bootstrap/blob/master/docs/candidates/offline-artifact-archive.md).
The first private archive has this identity:

```text
archive=pyflamestk-offline-artifacts-git-5b8368cc88d9.tar
sha256=b3c7f89ef88564a04ca1cb02bf85cb66c314452812cb4827db24a58dcc3c0e0a
```

## Restoration and verification

Offline restoration is manual and explicit. Before using a staged file, verify
its bytes against both `MANIFEST.tsv` and `OFFLINE_ARTIFACT_SHA256SUMS`. Do not
silently copy offline artifacts back into a release, test, wheel, or source
distribution. In particular, no repository release should redistribute a
pseudopotential without separately established redistribution authority.

Existing Git history remains unchanged and continues to identify earlier
committed bytes. This boundary changes the maintained tip; it is not a history
rewrite or a claim that the removed calculations are reproducible, correct, or
scientifically validated.
