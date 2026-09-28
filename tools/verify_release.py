#!/usr/bin/env python3
"""Verify a release archive from a fresh extraction, with nothing borrowed from this tree.

Running checks in the tree that produced an archive proves only that the bytes on disk are
the bytes that were hashed; it says nothing about whether someone who downloads the archive
can build and reproduce anything. So the archive is extracted into a scratch directory, the
build is done there, and every command runs with that directory as its root.

Steps, each reported pass/fail/skip and none of them skipped silently:

  1. the archive's own SHA-256 matches its published `.sha256`, if one is alongside;
  2. every entry matches `ARCHIVE_MANIFEST.sha256` inside it;
  3. the library compiles as strict C++11 under every available compiler;
  4. its regression suite passes under every available compiler;
  5. the Python validator's self-test passes (needs `uv`; skipped otherwise);
  6. Fig. 2 of the paper and its caption and manifest regenerate byte-identically from the
     archived results of the finite-precision study (`make figures` in
     experiments/exp4_finite_precision).

Usage:
    tools/make_release_archive.py --ref v3.0.0 --out dist/
    tools/verify_release.py dist/bi-kappa-v3.0.0-<sha>.tar.gz
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile

CHECK, CROSS = "PASS", "FAIL"


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str]] = []

    def add(self, ok: bool | None, name: str, detail: str = "") -> None:
        status = "SKIP" if ok is None else (CHECK if ok else CROSS)
        self.rows.append((status, name, detail))
        print(f"  [{status}] {name}" + (f" -- {detail}" if detail else ""), flush=True)

    @property
    def failed(self) -> bool:
        return any(s == CROSS for s, _, _ in self.rows)


def sh(cmd: list[str], cwd: str, timeout: int = 3600) -> tuple[int, str]:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout + p.stderr)[-4000:]
    except FileNotFoundError as exc:
        return 127, str(exc)
    except subprocess.TimeoutExpired:
        return 124, "timed out"


def tree_digests(root: str) -> dict[str, str]:
    """SHA-256 of every file under `root`, keyed by path relative to it."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            with open(full, "rb") as fh:
                out[os.path.relpath(full, root)] = hashlib.sha256(fh.read()).hexdigest()
    return out


def compilers() -> list[tuple[str, list[str]]]:
    found = []
    if shutil.which("clang++"):
        found.append(("clang++/libc++", ["clang++", "-stdlib=libc++"]))
    for g in ("g++-15", "g++-14", "g++-13", "g++"):
        if shutil.which(g):
            found.append((f"{g}/libstdc++", [g]))
            break
    return found


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Extract a release archive into a scratch directory and verify it there.")
    ap.add_argument("archive", help="the .tar.gz written by tools/make_release_archive.py")
    ap.add_argument("--keep", action="store_true", help="keep the extraction directory")
    ap.add_argument("--skip-python", action="store_true",
                    help="skip the Python validator's self-test")
    ap.add_argument("--skip-experiment", action="store_true",
                    help="skip regenerating Fig. 2 of the finite-precision study")
    args = ap.parse_args()

    archive = os.path.abspath(args.archive)
    rep = Report()
    print(f"verifying {archive}\n")

    # 1 -- archive digest
    digest = hashlib.sha256(open(archive, "rb").read()).hexdigest()
    side = archive + ".sha256"
    if os.path.exists(side):
        published = open(side).read().split()[0]
        rep.add(digest == published, "archive SHA-256 matches its published digest",
                f"{digest[:16]}...")
    else:
        rep.add(None, "archive SHA-256 matches its published digest",
                f"no {os.path.basename(side)} alongside; computed {digest[:16]}...")

    tmp = tempfile.mkdtemp(prefix="bikappa-verify-")
    try:
        with tarfile.open(archive, "r:gz") as tar:
            names = tar.getnames()
            root = os.path.commonprefix(names).split("/")[0]
            if hasattr(tarfile, "data_filter"):
                tar.extractall(tmp, filter="data")
            else:
                tar.extractall(tmp)
        top = os.path.join(tmp, root)
        rep.add(os.path.isdir(top), "archive extracts to a single root",
                f"{root} ({len(names)} entries)")

        # 2 -- internal manifest
        man = os.path.join(top, "ARCHIVE_MANIFEST.sha256")
        if os.path.exists(man):
            bad, n = [], 0
            for line in open(man):
                if line.startswith("#") or not line.strip():
                    continue
                want, rel = line.split(maxsplit=1)
                rel = rel.strip()
                full = os.path.join(top, rel)
                n += 1
                if not os.path.exists(full):
                    bad.append(f"{rel} missing")
                elif hashlib.sha256(open(full, "rb").read()).hexdigest() != want:
                    bad.append(f"{rel} differs")
            rep.add(not bad, "every archived file matches ARCHIVE_MANIFEST.sha256",
                    f"{n} entries" + (f"; {len(bad)} bad: {bad[:3]}" if bad else ""))
        else:
            rep.add(False, "every archived file matches ARCHIVE_MANIFEST.sha256",
                    "manifest absent from the archive")

        cc = compilers()
        rep.add(bool(cc), "a C++ compiler is available",
                ", ".join(n for n, _ in cc) or "none found")

        # 3 & 4 -- the released library, built and tested from the extraction
        cpp = os.path.join(top, "cpp")
        for label, base in cc:
            rc, out = sh(base + ["-std=c++11", "-pedantic-errors", "-Wall", "-Wextra",
                                 "-fsyntax-only", "-I", ".", "main.cpp"], cwd=cpp)
            rep.add(rc == 0, f"strict C++11 compile [{label}]", out.strip().splitlines()[-1]
                    if rc else "")
        if os.path.exists(os.path.join(cpp, "GNUmakefile")):
            for label, base in cc:
                env_cxx = base[0]
                rc, out = sh(["make", "-s", "test", f"CXX={env_cxx}"], cwd=cpp)
                if rc == 127 or "No rule to make target" in out:
                    rep.add(None, f"regression suite [{label}]", "no `test` target")
                else:
                    tail = [l for l in out.splitlines() if "PASS" in l or "FAIL" in l]
                    rep.add(rc == 0, f"regression suite [{label}]",
                            tail[-1] if tail else out.strip()[-200:])

        have_uv = shutil.which("uv") is not None

        # 5 -- the Python validator
        if args.skip_python:
            rep.add(None, "Python validator self-test", "--skip-python")
        elif not have_uv:
            rep.add(None, "Python validator self-test", "uv not found")
        else:
            rc, out = sh(["uv", "run", "python", "test_bikappa_validate.py"],
                         cwd=os.path.join(top, "python"))
            rep.add(rc == 0, "Python validator self-test", "" if rc == 0 else out.strip()[-300:])

        # 6 -- Fig. 2 from the archived results of the finite-precision study
        exp = os.path.join(top, "experiments", "exp4_finite_precision")
        name = "Fig. 2 regenerates byte-identically from the archived results"
        if args.skip_experiment:
            rep.add(None, name, "--skip-experiment")
        elif not os.path.isdir(exp):
            rep.add(None, name, "experiments/exp4_finite_precision not archived")
        elif not have_uv:
            rep.add(None, name, "uv not found")
        else:
            before = tree_digests(exp)
            rc, out = sh(["make", "-s", "figures"], cwd=exp)
            after = tree_digests(exp)
            changed = sorted(p for p in before if after.get(p) != before[p])
            if rc != 0:
                rep.add(False, name, out.strip()[-300:])
            else:
                rep.add(not changed, name,
                        f"{len(changed)} archived file(s) differ: {changed[:4]}" if changed
                        else f"{len(before)} archived files unchanged")
    finally:
        if args.keep:
            print(f"\nextraction kept at {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    print("\n" + ("VERIFICATION FAILED" if rep.failed else "VERIFICATION PASSED"))
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
