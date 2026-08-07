"""The vendored unit layer must stay byte-identical to the caller pipeline.

WES-PON-smk builds the normals that WES-snakemake normalises its tumours
against. If the two resolve read groups differently, duplicate marking differs,
coverage shifts, and every PON artifact is subtly wrong for the tumours it is
applied to -- while both pipelines keep running and reporting success. Nothing
downstream would catch it, so it is caught here.

The comparison is against the caller's published `main`, not against whatever
happens to sit in a sibling working tree: a local checkout can be mid-edit, on a
feature branch, or simply stale, and agreeing with it proves nothing. Whole
files are compared rather than individual objects, because unlike
`fastq_header.py` in the caller -- which vendors part of wesingest and adds to
it -- these two files are vendored entire and add nothing.

Upstream resolution, in order:

    1. `git fetch --depth 1 <WES_CALLER_URL> <WES_CALLER_REF>`, needs network
    2. `git show <ref>:<path>` in a local clone at `WES_CALLER_REPO`, needs the
       clone to have been fetched recently
    3. skip

Skipping keeps a vendored deployment installable with no access to the caller
repository, which is the same policy as the caller's own test_vendored_parser.py.
Set WES_CALLER_REF to a tag (`v1.1.0`) to pin the vendored copy to a caller
release; the default tracks `main` so a divergence goes red as soon as it lands.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import pytest

CALLER_URL = os.environ.get(
    "WES_CALLER_URL", "https://github.com/dmkv1/WES-snakemake.git")
CALLER_REF = os.environ.get("WES_CALLER_REF", "main")
CALLER_REPO = os.environ.get("WES_CALLER_REPO", "")
OFFLINE = os.environ.get("WES_CALLER_OFFLINE", "") not in ("", "0")

# Vendored entire, no local additions. Paths are the same in both repositories.
VENDORED = [
    "workflow/scripts/units.py",
    "workflow/scripts/fastq_header.py",
]

FETCH_TIMEOUT = 120


def _git(*args: str, cwd: Path | str | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                          text=True, timeout=FETCH_TIMEOUT)


class _Upstream:
    """A resolved upstream tree: where it came from, and how to read a file out."""

    def __init__(self, kind: str, cwd: Path, rev: str, origin: str):
        self.kind, self.cwd, self.rev, self.origin = kind, cwd, rev, origin
        self.commit = (_git("rev-parse", self.rev, cwd=cwd).stdout.strip()
                       or "unknown")

    def read(self, relpath: str) -> str:
        done = _git("show", f"{self.rev}:{relpath}", cwd=self.cwd)
        if done.returncode != 0:
            raise AssertionError(
                f"{relpath} is not present in {self.origin} at "
                f"{self.commit[:10]}: {done.stderr.strip()}")
        return done.stdout

    def __str__(self) -> str:
        return f"{self.origin} ({self.commit[:10]})"


def _from_remote(tmp: Path) -> _Upstream | None:
    if OFFLINE:
        return None
    if _git("init", "--bare", "--quiet", str(tmp)).returncode != 0:
        return None
    # A depth-1 fetch of one ref: the whole tree, none of the history.
    done = _git("fetch", "--depth", "1", "--quiet", CALLER_URL, CALLER_REF, cwd=tmp)
    if done.returncode != 0:
        return None
    return _Upstream("remote", tmp, "FETCH_HEAD", f"{CALLER_URL}@{CALLER_REF}")


def _from_local() -> _Upstream | None:
    """A sibling clone, read at its remote-tracking ref rather than its worktree."""
    default = Path(__file__).resolve().parent.parent.parent / "WES-snakemake-dev"
    repo = Path(CALLER_REPO) if CALLER_REPO else default
    if not (repo / ".git").exists():
        return None
    for rev in (f"origin/{CALLER_REF}", CALLER_REF):
        if _git("rev-parse", "--verify", "--quiet", rev, cwd=repo).returncode == 0:
            return _Upstream("local", repo, rev, f"{repo}@{rev}")
    return None


@pytest.fixture(scope="session")
def upstream(tmp_path_factory) -> _Upstream:
    tmp = tmp_path_factory.mktemp("caller")
    resolved = _from_remote(tmp) or _from_local()
    if resolved is None:
        pytest.skip(
            f"caller pipeline unreachable: could not fetch {CALLER_REF} from "
            f"{CALLER_URL}, and no clone at WES_CALLER_REPO")
    return resolved


@pytest.mark.parametrize("relpath", VENDORED)
def test_vendored_file_matches_the_caller(upstream, repo_root, relpath):
    ours = (repo_root / relpath).read_text(encoding="utf-8")
    theirs = upstream.read(relpath)
    assert ours == theirs, (
        f"{relpath} has diverged from {upstream}.\n"
        f"The caller's copy is canonical -- re-copy it here rather than editing "
        f"this file, then re-run the alignment comparison. Both pipelines must "
        f"resolve read groups identically or the panel of normals stops "
        f"matching the tumors it normalises."
    )


def test_the_upstream_ref_actually_resolved(upstream):
    """A silently-skipped identity test protects nothing; make the source visible."""
    assert upstream.commit != "unknown"
    print(f"\ncompared against {upstream} via {upstream.kind}")
