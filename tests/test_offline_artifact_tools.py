from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
STAGE_TOOL = REPOSITORY_ROOT / "tools/stage_offline_artifacts.py"


class OfflineArtifactStagingToolTest(unittest.TestCase):
    def test_stages_selected_artifacts_without_moving_authored_source(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository = Path(directory)
            self._git(repository, "init", "-q", "-b", "main")
            self._git(repository, "config", "user.name", "Tool Test")
            self._git(
                repository,
                "config",
                "user.email",
                "tool-test@example.invalid",
            )
            (repository / ".gitignore").write_text(".offline/\n", encoding="utf-8")
            (repository / "example").mkdir()
            source = repository / "example/input.config"
            artifact = repository / "example/results_001.out"
            payload = b"generated output\n"
            source.write_text("authored source\n", encoding="utf-8")
            artifact.write_bytes(payload)
            self._git(repository, "add", ".")
            self._git(repository, "commit", "-qm", "baseline")

            subprocess.run(
                [
                    sys.executable,
                    str(STAGE_TOOL),
                    "--repository-root",
                    str(repository),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            staged = repository / ".offline/release-readiness/example/results_001.out"
            self.assertTrue(source.is_file())
            self.assertFalse(artifact.exists())
            self.assertEqual(staged.read_bytes(), payload)
            checksum = hashlib.sha256(payload).hexdigest()
            self.assertEqual(
                (repository / "OFFLINE_ARTIFACT_SHA256SUMS").read_text(
                    encoding="utf-8"
                ),
                f"{checksum}  example/results_001.out\n",
            )
            self.assertIn(
                f"{checksum}\t{len(payload)}\t",
                (repository / ".offline/release-readiness/MANIFEST.tsv").read_text(
                    encoding="utf-8"
                ),
            )

    @staticmethod
    def _git(repository: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", str(repository), *arguments],
            check=True,
            capture_output=True,
            text=True,
        )


if __name__ == "__main__":
    unittest.main()
