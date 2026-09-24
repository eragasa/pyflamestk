# Offline artifact staging policy

`stage_offline_artifacts.py` preserves PyFlamestk's repository-specific
selection and safe-copy procedure. It does not run PyFlamestk, calculators,
schedulers, or retained scripts.

From a clean checkout whose `.gitignore` ignores `.offline/`:

```bash
python3 tools/stage_offline_artifacts.py \
  --repository-root . \
  --destination .offline/release-readiness \
  --checksum-manifest OFFLINE_ARTIFACT_SHA256SUMS
```

The tool uses an explicit PyFlamestk filename and generated-result policy. It
first copies selected tracked regular files into a temporary sibling directory,
verifies every copied digest and size, publishes the completed offline
directory, writes the tracked checksum manifest, and only then removes original
working-tree paths. It refuses dirty repositories, symlinks, a reused
destination, an existing checksum manifest, or a destination outside
`.offline/`.

The policy deliberately does not remove every `.dat` or `.out` file because
those extensions can contain authored inputs. Extend the reviewed selection
policy rather than broadening it implicitly.

## Generic verification and archival

Reusable manifest verification and deterministic tar creation belong to the
Project Koios bootstrap harness, not PyFlamestk. See:

- <https://github.com/eragasa/projectkoios-bootstrap/blob/master/docs/candidates/offline-artifact-archive.md>
- `projectkoios.bootstrap.harness.offline_artifacts`

Given an explicit bootstrap checkout, verify this staging tree with:

```bash
PYTHONPATH=/path/to/projectkoios-bootstrap/python python3.14 -m \
  projectkoios.bootstrap.harness.offline_artifacts verify \
  --staging-root .offline/release-readiness \
  --manifest .offline/release-readiness/MANIFEST.tsv \
  --source-checkout . \
  --source-revision 5b8368cc88d91bc56f9cc1c8a7fa8d9ea3d6b359 \
  --source-repository-url https://github.com/eragasa/pyflamestk
```

The ignored directory is operator-managed staging. Archive it to independent
private storage before deleting the local staging copy. The tracked checksum
manifest remains useful after the bytes leave the workstation.
