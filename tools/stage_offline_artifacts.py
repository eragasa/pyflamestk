from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path, PurePosixPath

_EXACT_NAMES = frozenset(
    {
        ".spyderproject",
        ".spyderworkspace",
        "CHG",
        "CHGCAR",
        "CONTCAR",
        "DOSCAR",
        "EIGENVAL",
        "IBZKPT",
        "OSZICAR",
        "OUTCAR",
        "PCDAT",
        "POTCAR",
        "WAVECAR",
        "XDATCAR",
        "job.log",
        "log.lammps",
        "out.dat",
        "param_out.dat",
        "params.dat",
        "pareto.out",
        "pyposmat.out",
        "regress_lammps.out",
        "restart.equil",
        "results.out",
        "vasprun.xml",
    }
)
_GENERATED_PATTERNS = (
    re.compile(r"^\._"),
    re.compile(r".*\.pyc$"),
    re.compile(r".*\.o\d+$"),
    re.compile(r".*\.e\d+$"),
    re.compile(r"(?:results.*|sim_results.*|params_\d+|.*_results_\d+)\.(?:out|dat)$"),
)


def classification_for(relative_path: PurePosixPath) -> str | None:
    """Return the bounded offline classification for one tracked path."""
    if relative_path.name in _EXACT_NAMES:
        return "known-generated-calculator-or-editor-artifact"
    if "__pycache__" in relative_path.parts:
        return "python-cache"
    if any(pattern.fullmatch(relative_path.name) for pattern in _GENERATED_PATTERNS):
        return "generated-result-or-cache-pattern"
    return None


def sha256_file(path: Path) -> tuple[str, int]:
    """Hash one non-symlink regular file without loading it all into memory."""
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"artifact must be a regular non-symlink file: {path}")
    digest = hashlib.sha256()
    byte_size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
            byte_size += len(chunk)
    return digest.hexdigest(), byte_size


def tracked_paths(repository_root: Path) -> tuple[PurePosixPath, ...]:
    result = subprocess.run(
        ["git", "-C", str(repository_root), "ls-files", "-z"],
        check=True,
        capture_output=True,
    )
    return tuple(
        PurePosixPath(value.decode("utf-8", errors="strict"))
        for value in result.stdout.split(b"\0")
        if value
    )


def stage(
    repository_root: Path,
    destination_relative: PurePosixPath,
    checksum_manifest_relative: PurePosixPath,
) -> tuple[int, int]:
    """Copy selected artifacts to ignored staging, verify, then remove originals."""
    root = repository_root.resolve(strict=True)
    _require_normalized_relative(destination_relative, "destination")
    _require_normalized_relative(checksum_manifest_relative, "checksum manifest")
    if destination_relative.parts[0] != ".offline":
        raise ValueError("destination must be below the ignored .offline directory")
    destination = root.joinpath(*destination_relative.parts)
    checksum_manifest = root.joinpath(*checksum_manifest_relative.parts)
    status = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain"],
        check=True,
        capture_output=True,
    )
    if status.stdout:
        raise ValueError("repository must be clean before staging artifacts")
    ignored = subprocess.run(
        ["git", "-C", str(root), "check-ignore", "--quiet", str(destination)],
        check=False,
    )
    if ignored.returncode != 0:
        raise ValueError("destination must be ignored by repository policy")
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"destination already exists: {destination}")
    if checksum_manifest.exists() or checksum_manifest.is_symlink():
        raise ValueError(f"checksum manifest already exists: {checksum_manifest}")

    records: list[tuple[str, PurePosixPath, int, str]] = []
    for relative_path in tracked_paths(root):
        reason = classification_for(relative_path)
        if reason is None:
            continue
        source = root.joinpath(*relative_path.parts)
        digest, byte_size = sha256_file(source)
        records.append((digest, relative_path, byte_size, reason))
    if not records:
        raise ValueError("selection policy found no tracked artifacts")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix="offline-artifacts-", dir=destination.parent)
    )
    try:
        for digest, relative_path, byte_size, _reason in records:
            source = root.joinpath(*relative_path.parts)
            copied = temporary.joinpath(*relative_path.parts)
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, copied, follow_symlinks=False)
            copied_digest, copied_size = sha256_file(copied)
            if (copied_digest, copied_size) != (digest, byte_size):
                raise ValueError(f"copied artifact identity changed: {relative_path}")
        temporary.joinpath("MANIFEST.tsv").write_text(
            "sha256\tbyte_size\treason\toriginal_path\n"
            + "".join(
                f"{digest}\t{byte_size}\t{reason}\t{relative_path.as_posix()}\n"
                for digest, relative_path, byte_size, reason in records
            ),
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    checksum_manifest.write_text(
        "".join(
            f"{digest}  {relative_path.as_posix()}\n"
            for digest, relative_path, _byte_size, _reason in records
        ),
        encoding="utf-8",
    )
    for _digest, relative_path, _byte_size, _reason in records:
        root.joinpath(*relative_path.parts).unlink()
    return len(records), sum(record[2] for record in records)


def _require_normalized_relative(path: PurePosixPath, field_name: str) -> None:
    value = path.as_posix()
    if (
        not value
        or value == "."
        or path.is_absolute()
        or "." in path.parts
        or ".." in path.parts
    ):
        raise ValueError(f"{field_name} must be a normalized relative POSIX path")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Stage selected tracked artifacts for operator-managed offline storage."
        )
    )
    result.add_argument("--repository-root", type=Path, required=True)
    result.add_argument(
        "--destination",
        type=PurePosixPath,
        default=PurePosixPath(".offline/release-readiness"),
    )
    result.add_argument(
        "--checksum-manifest",
        type=PurePosixPath,
        default=PurePosixPath("OFFLINE_ARTIFACT_SHA256SUMS"),
    )
    return result


def main(arguments: list[str] | None = None) -> int:
    args = parser().parse_args(arguments)
    count, byte_size = stage(
        args.repository_root,
        args.destination,
        args.checksum_manifest,
    )
    print(f"staged {count} files ({byte_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
