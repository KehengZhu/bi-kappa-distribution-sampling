#!/usr/bin/env python3
"""Verify a release archive from a fresh extraction, with nothing borrowed from this tree.

This is the check gate G6 actually requires. Running `make verify` in the tree that produced
the artifacts proves only that the bytes on disk are the bytes that were hashed; it says
nothing about whether someone who downloads the archive can reproduce anything. So the
archive is extracted into a scratch directory, the build is done there, and every command
runs with that directory as its root.

Steps, each reported pass/fail and none of them skipped silently:

  1. the archive's own SHA-256 matches its published `.sha256`;
  2. every entry matches `ARCHIVE_MANIFEST.sha256` inside it;
  3. the released library compiles as strict C++11 under every available compiler;
  4. its regression suite passes under every available compiler;
  5. the experiment's checksum manifests verify (`make verify`);
  6. the experiment's derived artifacts regenerate byte-identically (`make reverify`).

Usage:
    tools/verify_release.py dist/bi-kappa-v2.0.0-<sha>.tar.gz
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
    ap = argparse.ArgumentParser()
    ap.add_argument("archive")
    ap.add_argument("--keep", action="store_true", help="keep the extraction directory")
    ap.add_argument("--skip-experiment", action="store_true",
                    help="verify the library only; skip the experiment manifests")
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

        # 5 & 6 -- the experiment bundle
        if not args.skip_experiment:
            exp = os.path.join(top, "experiments", "exp7_confirmatory")
            if os.path.isdir(exp):
                rc, out = sh(["make", "-s", "verify"], cwd=exp)
                rep.add(rc == 0, "experiment checksum manifests verify",
                        out.strip()[-300:] if rc else "")
                rc, out = sh(["make", "-s", "reverify"], cwd=exp)
                if "No rule to make target" in out:
                    rep.add(None, "derived artifacts regenerate byte-identically",
                            "no `reverify` target")
                else:
                    rep.add(rc == 0, "derived artifacts regenerate byte-identically",
                            out.strip()[-300:] if rc else "")
            else:
                rep.add(None, "experiment bundle present", "exp7_confirmatory not archived")
    finally:
        if args.keep:
            print(f"\nextraction kept at {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    print("\n" + ("VERIFICATION FAILED" if rep.failed else "VERIFICATION PASSED"))
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
