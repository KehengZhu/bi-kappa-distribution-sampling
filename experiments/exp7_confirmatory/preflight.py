"""Experiment 6 preflight: record the environment, and refuse a dirty production run.

Run with:  uv run --project ../../python python preflight.py [--allow-dirty]

Experiment 4's saved result records a dirty tree whose runtime state cannot be
reconstructed from the commit it names, which is why its numbers are a historical baseline
rather than reproducible evidence.  This script exists so that Experiment 6 cannot repeat
that.  It writes ``raw/environment.json`` and exits non-zero unless every file whose
content determines a result is tracked in git and identical to its committed object.

"Dirty" is deliberately narrower than ``git status``: unrelated edits elsewhere in the
working tree do not affect what this experiment computes, and the plan requires them to be
preserved rather than stashed.  What must be clean is the dependency set below.  The
overall repository state is recorded either way, honestly, in ``git.repo_dirty``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import locale
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(HERE, "raw")

# Files whose content determines the *raw draws*.  Nothing else can change a number in
# raw/; editing an analysis script cannot retroactively alter a sample that was already
# written.  Keeping this set separate is what lets the manifest say precisely which
# sources a given raw file came from.
SAMPLING_DEPENDENCIES = [
    "../../cpp/bi_kappa_distribution.H",          # the candidate under test
    "src/legacy/bi_kappa_distribution_v1.H",      # the frozen comparator
    "src/exp7_common.H",
    "src/exp7_loaders.H",
    "src/exp7_digest.H",
    "src/exp7_probe.cpp",
    "src/exp7_protocol.H",                        # generated from the protocol
    "PROTOCOL.md",                                # the frozen rules
    "config/protocol.json",
]

# Files that determine what is derived from those draws.
ANALYSIS_DEPENDENCIES = [
    "src/exp7_oracle.cpp",
    "src/gen_protocol_header.py",
    "src/check_legacy.py",
    "GNUmakefile",
    "analyze.py",
    "make_figures.py",
    "preflight.py",
    "exp7_stats.py",
    "exp7_io.py",
    "exp7_families.py",
    "exp7_gates.py",
    "exp7_portability.py",
    "config/make_protocol.py",
]

# Every file whose content determines a number this experiment reports.
DEPENDENCIES = SAMPLING_DEPENDENCIES + ANALYSIS_DEPENDENCIES

# The frozen protocol, and the hash that ties this run to it.  Loaded once, at import, so a
# missing or malformed protocol stops the run before any data is produced rather than after.
with open(os.path.join(HERE, "config", "protocol.json"), "rb") as _fh:
    _PROTOCOL_BYTES = _fh.read()
PROTOCOL_SHA256 = hashlib.sha256(_PROTOCOL_BYTES).hexdigest()
PROTOCOL = json.loads(_PROTOCOL_BYTES)


def sh(*args: str, cwd: str | None = None) -> str:
    try:
        return subprocess.run(
            args, cwd=cwd or HERE, capture_output=True, text=True, check=False
        ).stdout.strip()
    except OSError:
        return ""


def sha256(path: str) -> str | None:
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return None


def git_blob_hash(path: str) -> str | None:
    """The SHA-256 of the committed content of ``path`` at HEAD, or None if untracked."""
    rel = os.path.relpath(os.path.abspath(os.path.join(HERE, path)), ROOT)
    blob = subprocess.run(
        ["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, check=False
    )
    if blob.returncode != 0:
        return None
    return hashlib.sha256(blob.stdout).hexdigest()


def probe_env(tag: str) -> dict | None:
    exe = os.path.join(HERE, f"exp7_probe_{tag}.exe")
    if not os.path.exists(exe):
        return None
    out = subprocess.run([exe, "env"], capture_output=True, text=True, check=False).stdout
    try:
        return json.loads(out)
    except json.JSONDecodeError:
        return None


def pkg_versions() -> dict:
    out = {}
    for name in ("numpy", "scipy", "matplotlib", "pandas", "mpmath"):
        try:
            mod = __import__(name)
            out[name] = getattr(mod, "__version__", "unknown")
        except ImportError:
            out[name] = None
    return out


def boost_version() -> str | None:
    for inc in (os.environ.get("BOOST_INC"), "/opt/homebrew/include", "/usr/local/include"):
        if not inc:
            continue
        path = os.path.join(inc, "boost", "version.hpp")
        if not os.path.exists(path):
            continue
        for line in open(path, encoding="utf-8", errors="replace"):
            if line.startswith("#define BOOST_LIB_VERSION"):
                return line.split('"')[1]
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--allow-dirty",
        action="store_true",
        help="permit a development run against an uncommitted dependency; the manifest is "
        "marked exploratory and must not feed a manuscript figure or table",
    )
    args = ap.parse_args()

    deps = {}
    unclean = []
    for rel in DEPENDENCIES:
        path = os.path.join(HERE, rel)
        disk = sha256(path)
        head = git_blob_hash(rel)
        state = (
            "missing"
            if disk is None
            else "untracked"
            if head is None
            else "clean"
            if head == disk
            else "modified"
        )
        deps[rel] = {"sha256": disk, "head_sha256": head, "state": state}
        if state != "clean":
            unclean.append(f"{rel} [{state}]")

    # Continuity with the previous preflight.  Provenance is about content, not about when
    # git happened to record it: if the dependency hashes are unchanged since the preflight
    # that a set of raw data was produced under, that data was produced from this content.
    prev_path = os.path.join(RAW, "environment.json")
    prev_hashes, prev_time = None, None
    if os.path.exists(prev_path):
        try:
            prev = json.load(open(prev_path, encoding="utf-8"))
            prev_hashes = {k: v.get("sha256") for k, v in (prev.get("dependencies") or {}).items()}
            prev_time = prev.get("generated_utc")
        except (OSError, json.JSONDecodeError):
            pass
    now_hashes = {k: v["sha256"] for k, v in deps.items()}
    unchanged = bool(prev_hashes is not None and prev_hashes == now_hashes)
    sampling_unchanged = bool(
        prev_hashes is not None and
        all(prev_hashes.get(k) == now_hashes.get(k) for k in SAMPLING_DEPENDENCIES))

    commit = sh("git", "rev-parse", "HEAD", cwd=ROOT)
    repo_dirty = bool(sh("git", "status", "--porcelain", cwd=ROOT))

    env = {
        "experiment": "exp7_confirmatory",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "exploratory": bool(unclean),
        "git": {
            "commit": commit,
            "repo_dirty": repo_dirty,
            "dependencies_clean": not unclean,
            "unclean_dependencies": unclean,
        },
        "dependencies": deps,
        "continuity": {
            "previous_preflight_utc": prev_time,
            "dependency_hashes_unchanged_since_previous_preflight": unchanged,
            "sampling_dependency_hashes_unchanged_since_previous_preflight":
                sampling_unchanged,
            "note": "Raw draws depend only on the sampling dependencies. If those hashes are "
                    "unchanged since the preflight a raw file was produced under, that file "
                    "came from this content, whatever else was edited afterwards.",
        },
        "sampling_dependencies": SAMPLING_DEPENDENCIES,
        "analysis_dependencies": ANALYSIS_DEPENDENCIES,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": sh("sysctl", "-n", "machdep.cpu.brand_string") or platform.processor(),
            "cpu_count": os.cpu_count(),
            "system": platform.system(),
            "release": platform.release(),
            "python": platform.python_version(),
            "locale": str(locale.setlocale(locale.LC_ALL)),
        },
        "toolchains": {
            "clang": sh("clang++", "--version").splitlines()[:2],
            "gcc": sh("g++-15", "--version").splitlines()[:1],
            "cxxflags": os.environ.get("EXP7_CXXFLAGS",
                        "-Wall -Wextra -std=c++11 -O2 -ffp-contract=off"),
            "oracle_flags": "-Wall -Wextra -std=c++14 -O2 -ffp-contract=off",
            "boost": boost_version(),
        },
        "floating_point": {tag: probe_env(tag) for tag in ("libcxx", "libstdcxx")},
        "python_packages": pkg_versions(),
        "commands": {
            "selftest": "make selftest",
            "preflight": "make preflight",
            "pilot": "make pilot",
            "baseline": "make baseline",
            "mechanism": "make mechanism",
            "oracle": "make oracle",
            "conditioning": "make conditioning",
            "loader": "make loader",
            "portability": "make portability",
            "performance": "make performance",
            "analyze": "make analyze",
            "figures": "make figures",
            "checksums": "make checksums",
            "verify": "make verify",
        },
        # Read from the frozen protocol rather than mirrored here.  Experiment 6 kept the
        # seed block in six places -- the probe, preflight, a config file, the analysis
        # prose, the README and the plan -- and its performance phase then derived five more
        # by arithmetic, so 4006-4010 appear in its manifest and in no declaration at all.
        # One source, hashed, is the fix.
        "frozen_sizes": {
            "n_scalar": PROTOCOL["matrix"]["n_scalar"],
            "n_mechanism": PROTOCOL["matrix"]["n_mechanism"],
            "n_conditioning": PROTOCOL["matrix"]["n_conditioning"],
            "n_loader": PROTOCOL["matrix"]["n_loader"],
            "production_seeds": PROTOCOL["seeds"]["production"],
            "performance_seeds": PROTOCOL["seeds"]["performance"],
            "kappa_ladder": PROTOCOL["matrix"]["kappa_ladder"],
        },
        "protocol": {
            "version": PROTOCOL["protocol_version"],
            "sha256": PROTOCOL_SHA256,
            "path": "config/protocol.json",
        },
    }

    os.makedirs(RAW, exist_ok=True)
    out = os.path.join(RAW, "environment.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(env, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"  wrote {os.path.relpath(out, HERE)}")

    if unclean:
        print("  dependency set is NOT clean:")
        for item in unclean:
            print(f"    - {item}")
        if not args.allow_dirty:
            print(
                "  refusing a production run.  Commit the dependency set, or pass\n"
                "  ALLOW_DIRTY_DEV=1 to make preflight for an exploratory run whose output\n"
                "  must not feed a manuscript figure or table."
            )
            return 1
        print("  ALLOW_DIRTY_DEV set: continuing, output marked exploratory.")
    else:
        print(f"  dependency set clean at {commit[:12]}"
              f"{' (unrelated working-tree edits present)' if repo_dirty else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
