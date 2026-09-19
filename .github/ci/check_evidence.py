#!/usr/bin/env python3
"""Coverage bookkeeping for Experiment 7's G4 evidence.  Not a statistical test.

PROTOCOL.md 5.3 and 6 say the same thing twice: an environment that did not run leaves
G4 *open*, is never recorded as a pass, and is never replaced by a zero.  A green CI
check is read as a pass, so this script exists to make a run that is missing an
environment fail rather than succeed quietly.

It answers only questions of presence and shape -- is the counter file there, is it
non-empty, does its environment record say the run was native, were the frozen sizes and
seeds used.  Every decision that involves a distribution belongs to
``analyze.py --portability-ingest``, which reads ``config/protocol.json``; nothing here
duplicates it.

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

SCHEMA = "exp7-portability-coverage/1"

# Fields a counter row may carry.  The probe's schema belongs to the experiment, not to
# CI, so a field that is absent is reported as unverified rather than treated as wrong --
# with one exception: `execution`, which lives in the environment record and is checked
# hard, because a translated run must never reach the F7 decision.
ROW_N_KEYS = ("n_attempted", "n_intended", "n")
ROW_SEED_KEYS = ("seed",)
ROW_DIGEST_KEYS = ("stream_sha256", "digest_sha256", "sha256", "sample_digest",
                   "stream_digest", "digest", "draw_digest", "stream_fnv1a", "fnv1a")


def read_jsonl(path, limit=None):
    rows, bad = [], 0
    with open(path, "r", errors="replace") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except ValueError:
                bad += 1
            if limit and len(rows) >= limit:
                break
    return rows, bad


def first_present(row, keys):
    for k in keys:
        if k in row:
            return k, row[k]
    return None, None


def tags_present(art_dir):
    """Every toolchain whose environment record is in this bundle.

    A hand run can contribute a standard library the CI matrix cannot: macOS has no
    libstdc++ toolchain on a hosted runner but does on the development host.  The bundle
    is therefore read for what it holds, not only for what was asked of it.  An extra
    toolchain cannot make an environment complete; it can supply a cross-library pair the
    matrix on its own does not.
    """
    if not os.path.isdir(art_dir):
        return []
    return [name[len("environment_"):-len(".json")]
            for name in sorted(os.listdir(art_dir))
            if name.startswith("environment_") and name.endswith(".json")]


def check_pair(art_dir, tag, frozen, required=True):
    """One (environment, toolchain) pair: what is on disk, and what it says about itself."""
    env_path = os.path.join(art_dir, "environment_%s.json" % tag)
    jsonl_path = os.path.join(art_dir, "portability_%s.jsonl" % tag)
    out = {"tag": tag, "required_by_matrix": required,
           "environment_record": env_path, "counters": jsonl_path,
           "problems": []}

    if not os.path.exists(env_path):
        out["problems"].append("environment record missing")
    else:
        try:
            with open(env_path, "r") as fh:
                env = json.load(fh)
        except ValueError as exc:
            out["problems"].append("environment record unparseable: %s" % exc)
            env = {}
        out["execution"] = env.get("execution")
        out["arch"] = env.get("arch")
        out["stdlib"] = env.get("stdlib")
        out["compiler"] = env.get("compiler")
        out["target_triple"] = env.get("target_triple")
        out["loader_version"] = env.get("loader_version")
        out["frozen_build_flags_verified"] = (
            env.get("build", {}).get("frozen_build_flags", {}).get("verified"))
        if out["execution"] != "native":
            out["problems"].append(
                "execution is %r; PROTOCOL.md 5.3 admits only native runs to the F7 "
                "decision" % out["execution"])
        if out["frozen_build_flags_verified"] is not True:
            out["problems"].append("frozen build flags not verified against the build log")

    if not os.path.exists(jsonl_path):
        out["problems"].append("counter file missing")
        return out

    if os.path.getsize(jsonl_path) == 0:
        out["problems"].append("counter file is empty")
        return out

    rows, bad = read_jsonl(jsonl_path)
    out["rows"] = len(rows)
    out["unparseable_rows"] = bad
    if bad:
        out["problems"].append("%d unparseable row(s)" % bad)
    if not rows:
        out["problems"].append("counter file holds no rows")
        return out

    n_key, _ = first_present(rows[0], ROW_N_KEYS)
    if n_key is None:
        out["frozen_size_check"] = "unverified: no attempt-count field in the row schema"
    else:
        sizes = sorted({r.get(n_key) for r in rows})
        out["attempt_counts"] = sizes
        want = frozen.get("n_scalar")
        if want is not None and sizes != [want]:
            out["problems"].append(
                "attempt counts %s do not match the frozen n_scalar %s" % (sizes, want))
        out["frozen_size_check"] = "verified" if sizes == [want] else "failed"

    s_key, _ = first_present(rows[0], ROW_SEED_KEYS)
    if s_key is None:
        out["seed_check"] = "unverified: no seed field in the row schema"
    else:
        seeds = sorted({int(r[s_key]) for r in rows if r.get(s_key) is not None})
        out["seeds"] = seeds
        want = frozen.get("production_seeds")
        if want is not None and seeds != sorted(want):
            out["problems"].append("seeds %s are not the frozen production block %s"
                                   % (seeds, sorted(want)))
        out["seed_check"] = "verified" if seeds == sorted(want or []) else "failed"

    d_key, _ = first_present(rows[0], ROW_DIGEST_KEYS)
    out["sample_digest_field"] = d_key
    if d_key is None:
        # Not fatal here.  PROTOCOL.md 5.3 predicts bitwise equality of the counters *and*
        # of the returned vectors; without a per-row digest the cross-stdlib test can only
        # reach the counters, and whether that is enough to close G4 is analyze.py's call.
        out["sample_digest_present"] = False
    else:
        out["sample_digest_present"] = True

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True,
                    help="directory holding one subdirectory per downloaded artifact")
    ap.add_argument("--expect", required=True,
                    help="JSON list of {env_id, arch, os, tags[]} the G4 matrix asked for")
    ap.add_argument("--protocol", required=True, help="config/protocol.json")
    ap.add_argument("--out", required=True, help="coverage.json to write")
    ap.add_argument("--allow-open", action="store_true",
                    help="report an open gate without failing (for a deliberately partial "
                         "run); the coverage file records it either way")
    args = ap.parse_args()

    with open(args.protocol, "r") as fh:
        cfg = json.load(fh)
    frozen = {"n_scalar": cfg["matrix"]["n_scalar"],
              "production_seeds": cfg["seeds"]["production"]}

    with open(args.expect, "r") as fh:
        expected = json.load(fh)

    results = []
    for env in expected:
        env_id = env["env_id"]
        art_dir = os.path.join(args.evidence, "portability-%s" % env_id)
        entry = {"env_id": env_id, "arch": env.get("arch"), "os": env.get("os"),
                 "runner_label": env.get("runner"), "artifact_dir": art_dir,
                 "artifact_present": os.path.isdir(art_dir), "pairs": []}
        if not entry["artifact_present"]:
            entry["reason_not_run"] = (
                "no artifact uploaded: the job did not run or did not complete")
        else:
            required = list(env["tags"])
            extra = [t for t in tags_present(art_dir) if t not in required]
            entry["extra_toolchains"] = extra
            for tag in required:
                entry["pairs"].append(check_pair(art_dir, tag, frozen, required=True))
            for tag in extra:
                entry["pairs"].append(check_pair(art_dir, tag, frozen, required=False))
        required_pairs = [p for p in entry["pairs"] if p["required_by_matrix"]]
        entry["complete"] = (entry["artifact_present"] and bool(required_pairs)
                             and all(not p["problems"] for p in required_pairs))
        results.append(entry)

    # Which comparisons the evidence on hand could support.  Counted per usable pair
    # rather than per complete environment: a toolchain that ran natively and cleanly is
    # available to analyze.py whether or not a sibling toolchain in the same environment
    # failed.  This is a statement about the evidence, not about the gate -- gate_state
    # below still turns on completeness.
    archs = {}
    for e in results:
        for p in e["pairs"]:
            if p["problems"] or p.get("execution") != "native" or not p.get("stdlib"):
                continue
            archs.setdefault(p.get("arch") or e.get("arch"), []).append(p["stdlib"])

    coverage = {
        "schema": SCHEMA,
        "gate": "G4",
        "family": "F7",
        "expected_environments": len(results),
        "complete_environments": sum(1 for e in results if e["complete"]),
        "environments": results,
        "stdlibs_per_arch": archs,
        "cross_stdlib_pairs_available": {a: sorted(set(v)) for a, v in archs.items()
                                         if len(set(v)) >= 2},
        "cross_arch_comparison_available": len([a for a in archs if archs[a]]) >= 2,
        "reproduction_commands": "experiments/exp7_confirmatory/results/portability_remote.md",
    }
    missing = [e["env_id"] for e in results if not e["complete"]]
    coverage["incomplete_environments"] = missing
    coverage["gate_state"] = "closable" if not missing else "open"

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(coverage, fh, indent=2)
        fh.write("\n")

    lines = ["## G4 portability evidence coverage", "",
             "| env_id | arch | os | artifact | toolchains | complete | problems |",
             "|---|---|---|---|---|---|---|"]
    for e in results:
        probs = "; ".join(p_ for pair in e["pairs"] for p_ in pair["problems"]) \
            or e.get("reason_not_run", "")
        lines.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            e["env_id"], e.get("arch"), e.get("os"),
            "yes" if e["artifact_present"] else "NO",
            ", ".join(p["tag"] + ("" if p["required_by_matrix"] else "*")
                      for p in e["pairs"]) or "-",
            "yes" if e["complete"] else "NO", probs or "-"))
    lines += ["", "`*` marks a toolchain the CI matrix does not require, contributed by a "
                  "hand run; it can supply a comparison but cannot complete an environment."]
    lines += ["", "Gate state: **%s**" % coverage["gate_state"], "",
              "An environment that did not run leaves G4 open with an exact command in "
              "`experiments/exp7_confirmatory/results/portability_remote.md`.  It is never "
              "recorded as a pass (PROTOCOL.md 5.3)."]
    summary = "\n".join(lines)
    print(summary)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a") as fh:
            fh.write(summary + "\n")

    if missing and not args.allow_open:
        sys.stderr.write(
            "\ncheck_evidence: G4 is OPEN. Missing or incomplete: %s\n"
            "This job fails deliberately: a green check would read as a passing gate.\n"
            % ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
