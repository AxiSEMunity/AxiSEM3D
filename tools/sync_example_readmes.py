#!/usr/bin/env python3

from __future__ import annotations

import shutil
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    examples_root = repo_root / "examples"
    docs_root = repo_root / "doc" / "sphinx" / "source" / "user" / "examples" / "examples"

    synced = []
    skipped = []

    for example_dir in sorted(path for path in examples_root.iterdir() if path.is_dir()):
        source_readme = example_dir / "README.md"
        target_readme = docs_root / example_dir.name / "README.md"

        if not source_readme.is_file() or not target_readme.parent.is_dir():
            skipped.append(example_dir.name)
            continue

        shutil.copyfile(source_readme, target_readme)
        synced.append(example_dir.name)

    print(f"Synced {len(synced)} example README.md files into the Sphinx tree.")
    if skipped:
        print("Skipped:", ", ".join(skipped))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())