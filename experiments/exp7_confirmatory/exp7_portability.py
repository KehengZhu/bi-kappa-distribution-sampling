#!/usr/bin/env python3
"""Ingest portability evidence and decide family F7 / gate G4.

Consumes the artifact bundle the CI workflow produces: a directory holding one
``portability-<env_id>/`` subdirectory per environment, each with an
``environment_<tag>.json`` record and a ``portability_<tag>.jsonl`` counter file, plus a
top-level ``coverage.json``.

Two claims, tested with the rule appropriate to each (``PROTOCOL.md`` §5.3):

* **Same architecture, different standard library.** The candidate draws its Gamma, normal
  and uniform variates from primitives defined in its own header, so its output is a function
  of the engine alone and the prediction is *bitwise equality* of every counter and of the
  digest over the returned vectors. There is no tolerance, because the claim is exact. This
  is a sharper test than an interval overlap, and it is only available because the standard
  library was taken out of the sampling path.
* **Different architecture.** ``log``, ``exp`` and ``pow`` are not correctly rounded and the
  implementations differ, so bitwise agreement is impossible and the prediction is equality
  of *rates*, tested as equivalence against a pre-registered margin rather than as a failure
  to reject a null of exact equality.

Translated execution -- Rosetta 2, or Docker with the Rosetta binfmt handler -- is recorded
and reported but never enters the decision, and an environment that did not run leaves the
gate OPEN rather than passing. A missing platform can never be mistaken for a passing one.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict

import exp7_families as F

# Counter fields that must agree bit for bit between standard libraries on one architecture.
# Continuous diagnostics are compared separately and are allowed to differ, since they are
# summaries over draws rather than part of the sampled stream.
COUNTER_FIELDS = (
    "n_attempted", "n_returned", "n_finite", "n_avoidable", "n_honest",
    "n_finite_but_wrong", "n_nonfinite", "cat_finite", "cat_honest_overflow",
    "cat_avoidable", "cat_denominator_zero", "cat_subnormal_denominator",
    "gamma_variates", "uniform_variates", "engine_calls", "cap_rejects",
)
DIGEST_FIELDS = ("vector_digest", "digest", "sha256", "fnv1a")


def load_bundle(root: str) -> tuple[list[dict], list[dict]]:
    """Return (environment records, counter rows) found under ``root``."""
    envs, rows = [], []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in sorted(filenames):
            full = os.path.join(dirpath, fn)
            if fn.startswith("environment") and fn.endswith(".json"):
                with open(full) as fh:
                    rec = json.load(fh)
                rec.setdefault("_source", os.path.relpath(full, root))
                envs.append(rec)
            elif fn.startswith("portability") and fn.endswith(".jsonl"):
                tag = fn[len("portability_"):-len(".jsonl")] if "_" in fn else fn
                with open(fh_path := full) as fh:
                    for line in fh:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            r = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        r.setdefault("stdlib_tag", tag)
                        r.setdefault("_source", os.path.relpath(fh_path, root))
                        rows.append(r)
    return envs, rows


def env_key(rec: dict) -> tuple:
    return (rec.get("arch"), rec.get("stdlib"), rec.get("compiler"), rec.get("os"))


def row_key(r: dict) -> tuple:
    """Everything that identifies a configuration, excluding the environment."""
    return (r.get("method"), r.get("precision"), r.get("kappa"), r.get("seed"),
            r.get("cap"), r.get("phase"))


def digest_of(r: dict):
    for f in DIGEST_FIELDS:
        if f in r:
            return r[f]
    return None


def compare_bitwise(rows_a: list[dict], rows_b: list[dict]) -> dict:
    """Exact comparison of every counter and digest for matching configurations."""
    a = {row_key(r): r for r in rows_a}
    b = {row_key(r): r for r in rows_b}
    common = sorted(set(a) & set(b), key=lambda k: tuple(str(x) for x in k))
    only_a, only_b = sorted(set(a) - set(b), key=str), sorted(set(b) - set(a), key=str)
    diffs = []
    compared_digests = 0
    for k in common:
        ra, rb = a[k], b[k]
        for f in COUNTER_FIELDS:
            if f in ra and f in rb and ra[f] != rb[f]:
                diffs.append({"config": [str(x) for x in k], "field": f,
                              "a": ra[f], "b": rb[f]})
        da, db = digest_of(ra), digest_of(rb)
        if da is not None and db is not None:
            compared_digests += 1
            if da != db:
                diffs.append({"config": [str(x) for x in k], "field": "digest",
                              "a": da, "b": db})
    return {"n_common": len(common), "n_only_a": len(only_a), "n_only_b": len(only_b),
            "n_digests_compared": compared_digests, "n_differences": len(diffs),
            "differences": diffs[:50], "identical": not diffs and bool(common)}


def ingest(root: str, proto: F.Protocol) -> dict:
    envs, rows = load_bundle(root)
    if not envs:
        return {"ok": False, "reason": f"no environment records under {root}",
                "environments": [], "stdlib_pairs": [], "arch_pairs": []}

    by_env = defaultdict(list)
    env_by_key = {}
    for e in envs:
        env_by_key[env_key(e)] = e
    for r in rows:
        k = (r.get("arch"), r.get("stdlib"), r.get("compiler"), r.get("os"))
        if k not in env_by_key:
            # Fall back to matching on the stdlib tag when the row does not carry the
            # full environment identity.
            cands = [ek for ek in env_by_key
                     if ek[1] and r.get("stdlib_tag", "").startswith(str(ek[1])[:4].lower())]
            k = cands[0] if len(cands) == 1 else k
        by_env[k].append(r)

    native, translated, pending = [], [], []
    for k, e in env_by_key.items():
        rec = {"arch": e.get("arch"), "os": e.get("os"), "compiler": e.get("compiler"),
               "stdlib": e.get("stdlib"), "execution": e.get("execution"),
               "available": True, "completed": bool(by_env.get(k)),
               "n_rows": len(by_env.get(k, [])),
               "reason_not_run": None if by_env.get(k) else "no counter rows in the artifact"}
        if not rec["completed"]:
            pending.append(rec)
        elif e.get("execution") == "native":
            native.append(rec)
        else:
            rec["reason_not_run"] = (
                "execution is not native; PROTOCOL.md 5.3 records translated results as "
                "corroborating evidence and excludes them from the decision")
            translated.append(rec)

    # --- cross-standard-library, same architecture: bitwise -----------------------------
    stdlib_pairs = []
    by_arch = defaultdict(list)
    for k, e in env_by_key.items():
        if e.get("execution") == "native" and by_env.get(k):
            by_arch[e.get("arch")].append(k)
    for arch, keys in sorted(by_arch.items(), key=lambda kv: str(kv[0])):
        keys = sorted(keys, key=lambda k: str(k[1]))
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                ka, kb = keys[i], keys[j]
                if ka[1] == kb[1]:
                    continue
                cand_a = [r for r in by_env[ka] if r.get("method") == "CANDIDATE"]
                cand_b = [r for r in by_env[kb] if r.get("method") == "CANDIDATE"]
                cmp_ = compare_bitwise(cand_a, cand_b)
                cmp_["label"] = f"{arch}: {ka[1]} vs {kb[1]}"
                stdlib_pairs.append(cmp_)

    # --- cross-architecture: equivalence of rates ---------------------------------------
    arch_pairs = []
    archs = sorted({e.get("arch") for k, e in env_by_key.items()
                    if e.get("execution") == "native" and by_env.get(k)}, key=str)
    for i in range(len(archs)):
        for j in range(i + 1, len(archs)):
            a1, a2 = archs[i], archs[j]
            rows1 = [r for k, e in env_by_key.items() if e.get("arch") == a1
                     and e.get("execution") == "native"
                     for r in by_env.get(k, []) if r.get("method") == "CANDIDATE"]
            rows2 = [r for k, e in env_by_key.items() if e.get("arch") == a2
                     and e.get("execution") == "native"
                     for r in by_env.get(k, []) if r.get("method") == "CANDIDATE"]
            agg1, agg2 = defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
            for src, agg in ((rows1, agg1), (rows2, agg2)):
                for r in src:
                    key = (r.get("precision"), r.get("kappa"))
                    agg[key][0] += int(r.get("n_nonfinite", r.get("cat_honest_overflow", 0)) or 0)
                    agg[key][1] += int(r.get("n_attempted", 0) or 0)
            for key in sorted(set(agg1) & set(agg2), key=str):
                k1, n1 = agg1[key]
                k2, n2 = agg2[key]
                arch_pairs.append({"label": f"{a1} vs {a2} | {key[0]} kappa={key[1]}",
                                   "k1": k1, "n1": n1, "k2": k2, "n2": n2})

    f7 = F.family_F7(proto, stdlib_pairs, arch_pairs)
    return {"ok": f7.passed, "environments": native + translated + pending,
            "native": native, "translated": translated, "pending": pending,
            "stdlib_pairs": stdlib_pairs, "arch_pairs": arch_pairs,
            "F7": f7.as_dict(), "_family": f7}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle", help="directory of portability-<env_id>/ artifacts")
    ap.add_argument("--out", default=None, help="write the ingest result as JSON here")
    args = ap.parse_args()
    proto = F.Protocol()
    res = ingest(args.bundle, proto)
    res.pop("_family", None)
    text = json.dumps(res, indent=2, default=str)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    print(text)

    f7 = res.get("F7") or {}
    print(f"\nF7: {'no disagreement found' if f7.get('passed') else 'FAILS'} -- "
          f"{f7.get('detail', '')}", file=sys.stderr)

    # F7 passing is necessary but not sufficient for G4.  The gate also requires that both
    # architectures actually ran natively; otherwise the cross-architecture claim is
    # untested and the gate is OPEN, not closed.  Translated execution never counts.
    archs = sorted({e.get("arch") for e in res.get("native", [])}, key=str)
    reasons = []
    if not f7.get("passed"):
        reasons.append("F7 found a disagreement")
    if res.get("pending"):
        reasons.append(f"{len(res['pending'])} environment(s) produced no rows")
    if len(archs) < 2:
        reasons.append(f"only {archs} has native results, so the cross-architecture claim "
                       "is untested")
    if res.get("translated"):
        print(f"note: {len(res['translated'])} translated environment(s) recorded as "
              "corroborating evidence and excluded from the decision", file=sys.stderr)
    if reasons:
        print("G4 is NOT closed: " + "; ".join(reasons), file=sys.stderr)
        return 1
    print("G4 closes on this evidence.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
