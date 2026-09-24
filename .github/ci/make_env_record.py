#!/usr/bin/env python3
"""Assemble the machine-readable environment record for one Experiment 7 P5 toolchain.

PROTOCOL.md 5.3 decides G4 by comparing counters across standard libraries and rates
across architectures.  Both comparisons are meaningless without the arithmetic the
numbers came from, and the plan (E5, item 3) lists what has to be on file: compiler,
standard-library version, target triple, optimization flags, floating-point contraction,
rounding mode, and FTZ/DAZ state.  This script writes exactly that, merging three
sources:

  * ``env_probe`` -- what only a compiled program can see (rounding mode, subnormal
    flushing, LDBL_MANT_DIG, the exact overflow thresholds, the standard-library version
    macros);
  * the toolchain itself -- ``--version`` and ``-dumpmachine``;
  * the host and the build log -- identity, and the command line the probe was actually
    built with.

One field is computed rather than copied: ``execution``.  It is ``native`` only when the
architecture the binary was compiled for matches the architecture the host reports and
the process is not running under a translator.  PROTOCOL.md 5.3 excludes anything else
from the F7 decision, so it must never be possible to record a translated run as native
by forgetting to say otherwise.

Standard library only; runs on any Python 3.8+ present on a hosted runner.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone

SCHEMA = "exp7-portability-environment/1"

# The build contract the P5 counters are only comparable under.  -ffp-contract=off is not
# cosmetic: a fused multiply-add changes the rounding of the very products whose overflow
# is being counted, so a build that contracts is measuring a different envelope.
DEFAULT_REQUIRED_FLAGS = ["-std=c++11", "-O2", "-ffp-contract=off"]


def run(cmd, cwd=None):
    """Return (rc, stdout+stderr) without raising; a missing tool is data, not a crash."""
    try:
        p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, timeout=120)
        return p.returncode, p.stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        return -1, "unavailable: %s" % exc


def first_line(text):
    return text.splitlines()[0].strip() if text else ""


def sha256_of(path):
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return None


def host_facts():
    facts = {
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
    }
    if os.path.exists("/etc/os-release"):
        rel = {}
        with open("/etc/os-release", "r") as fh:
            for line in fh:
                if "=" in line:
                    k, _, v = line.partition("=")
                    rel[k.strip()] = v.strip().strip('"')
        facts["os_release"] = {k: rel.get(k) for k in ("ID", "VERSION_ID", "PRETTY_NAME")}
    rc, out = run(["sw_vers"])
    if rc == 0:
        facts["sw_vers"] = out
    rc, out = run(["sysctl", "-n", "machdep.cpu.brand_string"])
    if rc == 0 and out:
        facts["cpu_model"] = out
    elif os.path.exists("/proc/cpuinfo"):
        with open("/proc/cpuinfo", "r") as fh:
            for line in fh:
                if line.lower().startswith(("model name", "cpu part", "hardware")):
                    facts["cpu_model"] = line.split(":", 1)[1].strip()
                    break
    facts["cpu_count"] = os.cpu_count()
    return facts


def translation_state(arch_compiled_for):
    """native | translated | unknown, with the evidence that decided it.

    Rosetta 2 and Docker/Rosetta can run an x86_64 binary on this project's arm64
    development host, and PROTOCOL.md 5.3 records such a run as corroborating evidence
    only.  The distinction is therefore computed here from two independent signals rather
    than passed in on the command line.
    """
    evidence = {}
    host_machine = platform.machine().lower()
    host_arch = {"aarch64": "arm64", "arm64": "arm64", "x86_64": "x86_64",
                 "amd64": "x86_64"}.get(host_machine, host_machine)
    evidence["host_machine"] = host_machine
    evidence["host_arch"] = host_arch
    evidence["arch_compiled_for"] = arch_compiled_for

    rc, out = run(["sysctl", "-n", "sysctl.proc_translated"])
    if rc == 0 and out.strip() in ("0", "1"):
        evidence["sysctl.proc_translated"] = out.strip()
        if out.strip() == "1":
            return "translated", evidence

    if arch_compiled_for in ("arm64", "x86_64") and host_arch in ("arm64", "x86_64"):
        if arch_compiled_for != host_arch:
            return "translated", evidence
        return "native", evidence
    return "unknown", evidence


def extract_build_command(build_log, tag):
    """The compile line that produced this tag's probe, taken from the build log.

    Recorded rather than reconstructed: the flags that matter are the ones make actually
    passed, and a record that restates the flags the workflow *meant* to use would not be
    evidence of anything.
    """
    if not build_log or not os.path.exists(build_log):
        return None, "build log not captured"
    needle = "exp7_probe_%s" % tag
    found = None
    with open(build_log, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if needle in line and ("-o" in line or "-c" in line):
                found = line.strip()
    if found is None:
        return None, "no compile line naming %s in the build log" % needle
    return found, None


def flags_from_command(cmd):
    opt = [t for t in cmd.split() if re.fullmatch(r"-O[0-9a-zA-Z]*", t)]
    contract = [t for t in cmd.split() if t.startswith("-ffp-contract")]
    fast_math = [t for t in cmd.split() if t in ("-ffast-math", "-Ofast", "-funsafe-math-optimizations")]
    return opt, contract, fast_math


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--probe-json", required=True,
                    help="file holding the JSON object env_probe wrote to stdout")
    ap.add_argument("--tag", required=True, help="stdlib tag, e.g. libcxx or libstdcxx")
    ap.add_argument("--env-id", required=True, help="environment id, e.g. linux-x86_64")
    ap.add_argument("--cxx", required=True, help="the compiler command make was given")
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--build-log", default=None)
    ap.add_argument("--phase", default="P5")
    ap.add_argument("--out", required=True)
    ap.add_argument("--require-flag", action="append", default=None,
                    help="repeatable; defaults to the frozen set %s" % DEFAULT_REQUIRED_FLAGS)
    args = ap.parse_args()

    required = args.require_flag if args.require_flag else list(DEFAULT_REQUIRED_FLAGS)
    root = os.path.abspath(args.repo_root)
    exp = os.path.join(root, "experiments", "exp4_finite_precision")

    with open(args.probe_json, "r") as fh:
        probe = json.load(fh)

    cxx_argv = shlex.split(args.cxx)
    rc_v, ver = run(cxx_argv + ["--version"])
    rc_m, triple = run(cxx_argv + ["-dumpmachine"])
    rc_d, dumpver = run(cxx_argv + ["-dumpversion"])

    build_cmd, build_why = extract_build_command(args.build_log, args.tag)
    opt, contract, fast_math = flags_from_command(build_cmd) if build_cmd else ([], [], [])

    if build_cmd:
        missing = [f for f in required if f not in build_cmd.split()]
        flag_check = {"required": required, "missing": missing,
                      "verified": not missing,
                      "source": "compile line captured from the build log"}
    else:
        flag_check = {"required": required, "missing": None, "verified": False,
                      "source": "not verified: %s" % build_why}

    execution, exec_evidence = translation_state(probe.get("arch_compiled_for", "unknown"))

    protocol_md = os.path.join(exp, "PROTOCOL.md")
    protocol_json = os.path.join(exp, "config", "protocol.json")
    frozen = {}
    try:
        with open(protocol_json, "r") as fh:
            cfg = json.load(fh)
        frozen = {
            "n_scalar": cfg["matrix"]["n_scalar"],
            "production_seeds": cfg["seeds"]["production"],
            "log_ratio_margin": cfg["gates"]["G4"]["log_ratio_margin"],
            "alpha_F7": cfg["statistics"]["family_alpha"]["F7_portability"],
        }
    except (OSError, KeyError, ValueError) as exc:
        frozen = {"error": "could not read config/protocol.json: %s" % exc}

    sources = {}
    for rel in ("cpp/bi_kappa_distribution.H",
                "experiments/exp4_finite_precision/PROTOCOL.md",
                "experiments/exp4_finite_precision/config/protocol.json",
                "experiments/exp4_finite_precision/src/exp7_common.H",
                "experiments/exp4_finite_precision/src/exp7_loaders.H",
                "experiments/exp4_finite_precision/src/exp7_probe.cpp",
                "experiments/exp4_finite_precision/GNUmakefile"):
        digest = sha256_of(os.path.join(root, rel))
        if digest:
            sources[rel] = digest

    rc_c, commit = run(["git", "-C", root, "rev-parse", "HEAD"])
    rc_s, status = run(["git", "-C", root, "status", "--porcelain"])

    record = {
        "schema": SCHEMA,
        "experiment": "exp7_confirmatory",
        "phase": args.phase,
        "family": "F7",
        "gate": "G4",
        "env_id": args.env_id,
        "tag": args.tag,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),

        # The field PROTOCOL.md 5.3 turns on.  Computed, never asserted.
        "execution": execution,
        "execution_evidence": exec_evidence,

        "arch": probe.get("arch_compiled_for", "unknown"),
        "compiler": probe.get("compiler"),
        "compiler_family": probe.get("compiler_family"),
        "compiler_version_string": first_line(ver) if rc_v == 0 else "unavailable",
        "compiler_dumpversion": dumpver if rc_d == 0 else "unavailable",
        "target_triple": triple if rc_m == 0 else "unavailable",
        "stdlib": probe.get("stdlib"),
        "stdlib_version": probe.get("stdlib_version"),
        "stdlib_release": probe.get("stdlib_release"),

        "build": {
            "cxx_as_passed_to_make": args.cxx,
            "command_line": build_cmd,
            "command_line_note": build_why,
            "optimization_flags": opt,
            "fp_contract": contract if contract else ["(compiler default: not specified)"],
            "fast_math_flags": fast_math,
            "frozen_build_flags": flag_check,
        },

        "floating_point": {
            "rounding_mode": probe.get("rounding_mode"),
            "flt_eval_method": probe.get("flt_eval_method"),
            "float_subnormals_flushed": probe.get("float_subnormals_flushed"),
            "double_subnormals_flushed": probe.get("double_subnormals_flushed"),
            "ftz_daz": probe.get("ftz_daz"),
            "fp_fast_fma": probe.get("fp_fast_fma"),
            "flt_radix": probe.get("flt_radix"),
            "flt_mant_dig": probe.get("flt_mant_dig"),
            "dbl_mant_dig": probe.get("dbl_mant_dig"),
            "ldbl_mant_dig": probe.get("ldbl_mant_dig"),
            "sizeof_long_double": probe.get("sizeof_long_double"),
            "dbl_max_hex": probe.get("dbl_max_hex"),
            "flt_max_hex": probe.get("flt_max_hex"),
            "log_overflow_threshold_double": probe.get("log_overflow_threshold_double"),
            "log_overflow_threshold_float": probe.get("log_overflow_threshold_float"),
            "log_max_double": probe.get("log_max_double"),
            "log_max_float": probe.get("log_max_float"),
        },

        "rng": {
            "engine": probe.get("rng"),
            "min": probe.get("mt19937_min"),
            "max": probe.get("mt19937_max"),
        },

        "loader_version": probe.get("bi_kappa_version"),
        "cplusplus": probe.get("cplusplus"),

        "frozen_protocol": frozen,
        "source_sha256": sources,

        "git": {
            "commit": commit if rc_c == 0 else "unavailable",
            "repo_dirty": bool(status) if rc_s == 0 else None,
        },

        "ci": {
            "provider": "github-actions" if os.environ.get("GITHUB_ACTIONS") else "local",
            "workflow": os.environ.get("GITHUB_WORKFLOW"),
            "run_id": os.environ.get("GITHUB_RUN_ID"),
            "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
            "repository": os.environ.get("GITHUB_REPOSITORY"),
            "ref": os.environ.get("GITHUB_REF"),
            "sha": os.environ.get("GITHUB_SHA"),
            "runner_os": os.environ.get("RUNNER_OS"),
            "runner_arch": os.environ.get("RUNNER_ARCH"),
            "image_os": os.environ.get("ImageOS"),
            "image_version": os.environ.get("ImageVersion"),
        },

        "host": host_facts(),
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=False)
        fh.write("\n")

    print("wrote %s" % args.out)
    print("  execution      : %s (%s)" % (execution, exec_evidence))
    print("  arch / stdlib  : %s / %s" % (record["arch"], record["stdlib"]))
    print("  target triple  : %s" % record["target_triple"])
    print("  LDBL_MANT_DIG  : %s" % record["floating_point"]["ldbl_mant_dig"])
    print("  build flags    : %s" % flag_check)

    if execution != "native":
        sys.stderr.write(
            "make_env_record: execution is %r, not 'native'.  PROTOCOL.md 5.3 excludes it "
            "from the F7 decision; it is corroborating evidence only.\n" % execution)
    if flag_check["missing"]:
        sys.stderr.write("make_env_record: frozen build flags missing from the compile "
                         "line: %s\n" % flag_check["missing"])
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
