from __future__ import annotations

import argparse
import zipfile
from pathlib import Path, PurePosixPath


class WheelBoundaryError(ValueError):
    """Raised when a wheel violates the PyFlamestk release boundary."""


def verify_wheel(repository_root: Path, wheel: Path) -> tuple[str, ...]:
    """Verify that a wheel contains only package and distribution metadata."""
    root = repository_root.resolve(strict=True)
    wheel_path = wheel.resolve(strict=True)
    expected_modules = {
        f"pyflamestk/{path.name}" for path in (root / "pyflamestk").glob("*.py")
    }
    with zipfile.ZipFile(wheel_path) as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise WheelBoundaryError("wheel contains duplicate paths")
        for info in infos:
            path = PurePosixPath(info.filename)
            if path.is_absolute() or "." in path.parts or ".." in path.parts:
                raise WheelBoundaryError(
                    f"wheel contains an unsafe path: {info.filename}"
                )
            if info.is_dir():
                continue
            if not (
                info.filename.startswith("pyflamestk/")
                or ".dist-info/" in info.filename
            ):
                raise WheelBoundaryError(
                    f"wheel contains content outside its boundary: {info.filename}"
                )
            if info.filename.endswith((".pyc", ".pyo")):
                raise WheelBoundaryError(
                    f"wheel contains interpreter output: {info.filename}"
                )
    packaged_modules = {
        name
        for name in names
        if name.startswith("pyflamestk/") and name.endswith(".py")
    }
    if packaged_modules != expected_modules:
        raise WheelBoundaryError(
            "wheel module set does not match the explicit source package"
        )
    return tuple(sorted(names))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Verify the PyFlamestk wheel boundary."
    )
    result.add_argument("wheel", type=Path)
    result.add_argument("--repository-root", type=Path, default=Path.cwd())
    return result


def main(arguments: list[str] | None = None) -> int:
    args = parser().parse_args(arguments)
    for name in verify_wheel(args.repository_root, args.wheel):
        print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
