from __future__ import annotations

import re
import subprocess
from pathlib import Path, PurePosixPath

import tomllib

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_COMMIT = "3063f031301c04e395cb9d7bc5064f2a52e1110c"
HISTORICAL_TREE = "fd9788560365e93faae596ac91f2530aa6bd92a8"
PACKAGE_SUBTREE = "54866bd536c8131b1932354d61f1ff33ae881a69"
_CHECKSUM = re.compile(r"([0-9a-f]{64})  ([^\r\n]+)")


def test_offline_selection_policy_is_declarative_and_bounded() -> None:
    policy = tomllib.loads(
        (REPOSITORY_ROOT / "devtools/offline-artifacts.toml").read_text(
            encoding="utf-8"
        )
    )

    assert set(policy) == {"schema_version", "rules"}
    assert policy["schema_version"] == 1
    assert policy["rules"]
    assert all(set(rule) == {"kind", "reason", "values"} for rule in policy["rules"])


def test_offline_checksum_evidence_is_canonical_and_absent() -> None:
    lines = (
        (REPOSITORY_ROOT / "OFFLINE_ARTIFACT_SHA256SUMS")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    paths: list[str] = []
    for line in lines:
        match = _CHECKSUM.fullmatch(line)
        assert match is not None
        raw_path = match.group(2)
        path = PurePosixPath(raw_path)
        assert not path.is_absolute()
        assert "." not in path.parts
        assert ".." not in path.parts
        paths.append(raw_path)

    assert len(paths) == 3_191
    assert paths == sorted(paths)
    assert len(paths) == len(set(paths))
    assert sum(PurePosixPath(path).name == "POTCAR" for path in paths) == 30
    assert not any((REPOSITORY_ROOT / path).exists() for path in paths)


def test_historical_commit_and_package_subtree_are_distinguished() -> None:
    repository_available = (
        subprocess.run(
            ["git", "-C", str(REPOSITORY_ROOT), "cat-file", "-e", HISTORICAL_COMMIT],
            check=False,
            capture_output=True,
        ).returncode
        == 0
    )
    if not repository_available:
        return

    commit_tree = subprocess.run(
        [
            "git",
            "-C",
            str(REPOSITORY_ROOT),
            "rev-parse",
            f"{HISTORICAL_COMMIT}^{{tree}}",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    package_subtree = subprocess.run(
        [
            "git",
            "-C",
            str(REPOSITORY_ROOT),
            "rev-parse",
            f"{HISTORICAL_COMMIT}:src/pymatmc2",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    assert commit_tree == HISTORICAL_TREE
    assert package_subtree == PACKAGE_SUBTREE
    assert commit_tree != package_subtree
