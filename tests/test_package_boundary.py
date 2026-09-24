from __future__ import annotations

import ast
import tomllib
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPOSITORY_ROOT / "pyflamestk"


class PackageBoundaryTest(unittest.TestCase):
    def test_all_packaged_python_modules_parse(self) -> None:
        failures = []
        for path in sorted(PACKAGE_ROOT.glob("*.py")):
            try:
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            except SyntaxError as error:
                failures.append(f"{path.name}:{error.lineno}: {error.msg}")

        self.assertEqual(failures, [])

    def test_build_metadata_selects_only_pyflamestk(self) -> None:
        configuration = tomllib.loads(
            (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )

        self.assertEqual(
            configuration["tool"]["setuptools"]["packages"],
            ["pyflamestk"],
        )
        self.assertEqual(configuration["project"]["license"], "BSD-2-Clause")
        self.assertEqual(
            configuration["project"]["dependencies"],
            ["numpy>=1.24", "scipy>=1.10"],
        )
        self.assertEqual(
            configuration["project"]["optional-dependencies"]["dakota"],
            ["PyYAML>=6"],
        )
        self.assertEqual(
            configuration["project"]["optional-dependencies"]["plot"],
            ["matplotlib>=3.7"],
        )

    def test_dakota_yaml_loading_is_safe(self) -> None:
        module = ast.parse(
            (PACKAGE_ROOT / "dakota_interface.py").read_text(encoding="utf-8")
        )
        yaml_calls = {
            node.func.attr
            for node in ast.walk(module)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "yaml"
        }

        self.assertIn("safe_load", yaml_calls)
        self.assertIn("safe_dump", yaml_calls)
        self.assertNotIn("load", yaml_calls)


if __name__ == "__main__":
    unittest.main()
