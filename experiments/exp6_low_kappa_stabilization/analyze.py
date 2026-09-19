"""Experiment 6 analysis: raw counters and bulk draws in, frozen source data out.

Run with:  uv run --project ../../python python analyze.py

Nothing here recomputes a simulation.  Every number written under ``results/`` comes from
``raw/``, and every number plotted by ``make_figures.py`` comes from a committed CSV
written here, so each plotted point is recoverable from a human-readable row.

Reading order matches the runbook: E0 baseline, E1 pilot, E2 mechanism, E3 conditioning,
E4 complete loader, E5 portability, E6 performance, then the tables and the gate report.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy import special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exp6_io as IO       # noqa: E402
import exp6_stats as S     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(HERE, "raw")
RESULTS = os.path.join(HERE, "results")
CONFIG = os.path.join(HERE, "config")
EXP4 = os.path.join(ROOT, "experiments", "exp4_precision", "results", "exp4_results.json")

SMOKE = False


def out(name: str) -> str:
    os.makedirs(RESULTS, exist_ok=True)
    return os.path.join(RESULTS, name)


def write_csv(name: str, rows: list[dict], fieldnames: list[str] | None = None) -> str:
    path = out(name)
    if not rows:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        return path
    if fieldnames is None:
        # Union of every row's keys, in first-seen order.  Taking the first row's keys
        # would silently drop columns that only later rows carry, which is exactly how a
        # figure ends up reading a column that the CSV does not contain.
        fieldnames = []
        seen = set()
        for r in rows:
            for k in r:
                if k not in seen:
                    seen.add(k)
                    fieldnames.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in fieldnames})
    print(f"  wrote results/{name}  ({len(rows)} rows)")
    return path


def fmt(v):
    """Missing is ``NA``, never an empty cell.  A blank must never mean both zero and
    unavailable."""
    if v is None:
        return "NA"
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return repr(v)
    if isinstance(v, bool):
        return "true" if v else "false"
    return v


def shape_of(kappa: float, precision: str) -> float:
    """The shape the run actually used: ``kappa - 1/2`` formed in the working precision."""
    if precision == "float":
        return float(np.float32(kappa) - np.float32(0.5))
    return float(kappa - 0.5)


def group(rows, keys):
    g = collections.defaultdict(list)
    for r in rows:
        g[tuple(r.get(k) for k in keys)].append(r)
    return g


# ===========================================================================
# E0 - clean baseline reproduction
# ===========================================================================
def phase_baseline(report: dict) -> list[dict]:
    rows = IO.load_phase(RAW, "baseline", SMOKE)
    if not rows:
        print("  E0: no baseline rows")
        return []

    pooled = []
    for key, rs in group(rows, ["method", "precision", "stdlib", "tag", "kappa"]).items():
        method, precision, stdlib, tag, kappa = key
        n = sum(r["n_attempted"] for r in rs)
        nonfinite = sum(r["nonfinite_output"] for r in rs)
        est = S.rate_with_interval(nonfinite, n)
        pooled.append({
            "method": method, "precision": precision, "stdlib": stdlib, "tag": tag,
            "kappa": kappa,
            "shape_a": shape_of(kappa, precision), "n_seeds": len(rs), "n_attempted": n,
            "n_nonfinite": nonfinite, "rate": est["rate"], "ci_lo": est["lo"],
            "ci_hi": est["hi"], "is_upper_bound": est["is_upper_bound"],
            "x2_zero": sum(r.get("x2_zero", 0) for r in rs) or None,
            "max_finite_log_component": max(
                (r["max_finite_log_component"] for r in rs
                 if isinstance(r["max_finite_log_component"], float)), default=float("nan")),
        })
    pooled.sort(key=lambda r: (r["method"], r["precision"], r["stdlib"], r["tag"],
                               r["kappa"]))

    # Comparison against the historical Experiment 4 result.  exp4 is read, never written.
    comparisons = []
    if os.path.exists(EXP4):
        with open(EXP4, encoding="utf-8") as fh:
            exp4 = json.load(fh)
        e4 = {}
        for rec in exp4.get("released", []):
            e4[(rec["precision"], rec["stdlib"], round(float(rec["kappa"]), 6))] = rec
        # Rows produced by the provenance control build are not production rows; they are
        # the explanation set, looked up only where the primary build does not reproduce.
        control = {(r["precision"], round(float(r["kappa"]), 6)): r for r in pooled
                   if r["method"] == "SPLIT_released" and "contractfast" in str(r["tag"])}
        for p in pooled:
            if p["method"] != "SPLIT_released" or "contractfast" in str(p["tag"]):
                continue
            stdlib = p["stdlib"]
            k = (p["precision"], stdlib, round(float(p["kappa"]), 6))
            rec = e4.get(k)
            if rec is None:
                continue
            p4 = float(rec["nonfinite_frac"])
            n4 = int(rec["n_draws"])
            p6 = p["rate"]
            n6 = p["n_attempted"]
            pbar = (p4 * n4 + p6 * n6) / (n4 + n6)
            se = math.sqrt(max(pbar * (1 - pbar) * (1 / n4 + 1 / n6), 0.0))
            z = (p6 - p4) / se if se > 0 else 0.0
            identical = bool(int(round(p4 * n4)) == p["n_nonfinite"] and n4 == n6)
            ctrl = control.get((p["precision"], round(float(p["kappa"]), 6)))
            explained = bool(
                not identical and ctrl is not None and ctrl["stdlib"] == stdlib and
                int(round(p4 * n4)) == ctrl["n_nonfinite"] and n4 == ctrl["n_attempted"])
            comparisons.append({
                "precision": p["precision"], "stdlib": stdlib, "kappa": p["kappa"],
                "exp4_rate": p4, "exp4_n": n4, "exp6_rate": p6, "exp6_n": n6,
                "abs_difference": p6 - p4, "standardized_difference": z,
                "identical_counts": identical,
                "identical_under_contract_fast": explained,
            })

    exact = [c for c in comparisons if c["identical_counts"]]
    explained = [c for c in comparisons if c["identical_under_contract_fast"]]
    unexplained = [c for c in comparisons
                   if not c["identical_counts"] and not c["identical_under_contract_fast"]]
    qualitative = []
    for c in comparisons:
        # A qualitative change is a cell crossing between "failures observed" and "none
        # observed"; that is the statement the manuscript's operating range rests on.
        was = c["exp4_rate"] > 0
        now = c["exp6_rate"] > 0
        if was != now:
            qualitative.append(c)

    if not comparisons:
        verdict = "not comparable"
    elif len(exact) == len(comparisons):
        verdict = "reproduced"
    elif not unexplained and not qualitative:
        verdict = "reproduced, residual difference attributed to floating-point contraction"
    elif not qualitative:
        verdict = "revised"
    else:
        verdict = "unresolved"

    lines = [
        "# E0 - clean baseline reproduction",
        "",
        f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.",
        "",
        "Experiment 4's saved result records a dirty working tree, and the clean git object",
        "at the commit it names does not reproduce the sampler-header hash it recorded, so its",
        "numbers are a historical baseline rather than reproducible evidence.  This phase reruns",
        "the same protocol - same seeds 4001-4005, same 10^6 attempts per seed, same isotropic",
        "unrotated uncapped configuration, same released class - from a recorded source state,",
        "and states plainly whether the numbers come back.",
        "",
        f"**Verdict: {verdict}.**",
        "",
        f"Cells compared: {len(comparisons)}.  Bitwise-identical counts: {len(exact)}.",
        f"Cells that are bitwise identical once the historical build's floating-point "
        f"contraction setting is restored: {len(explained)}.  Cells still unexplained: "
        f"{len(unexplained)}.",
        "",
        "Experiment 4 was built without an explicit `-ffp-contract` flag and so inherited the",
        "compiler default, which enables fused multiply-add; this experiment disables it, "
        "because",
        "contraction changes the rounding of exactly the products whose overflow is under "
        "study.",
        "That setting was not recorded in the Experiment 4 manifest. `make "
        "baseline-contract-check`",
        "rebuilds the libstdc++ probe with the historical setting and reruns the "
        "single-precision",
        "baseline; the rows it produces carry the tag `libstdcxx-contractfast` and are a "
        "provenance",
        "control, not production data.",
        f"Cells that change a qualitative operating-range statement: {len(qualitative)}.",
        "",
        "Same seed and same class means the two runs should agree bit for bit unless the",
        "released header changed between them.  A standardized difference is reported as well,",
        "because agreement across standard libraries is statistical, not bitwise.",
        "",
        "| precision | stdlib | kappa | Exp4 rate | Exp6 rate | difference | z | identical |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for c in sorted(comparisons, key=lambda c: (c["precision"], c["stdlib"], c["kappa"])):
        lines.append(
            f"| {c['precision']} | {c['stdlib']} | {c['kappa']:g} | {c['exp4_rate']:.6g} | "
            f"{c['exp6_rate']:.6g} | {c['abs_difference']:+.3g} | {c['standardized_difference']:+.2f} | "
            f"{'yes' if c['identical_counts'] else ('contraction' if c['identical_under_contract_fast'] else 'no')} |")
    if qualitative:
        lines += ["", "## Cells that change a qualitative statement", ""]
        for c in qualitative:
            lines.append(f"- {c['precision']} / {c['stdlib']} / kappa={c['kappa']:g}: "
                         f"{c['exp4_rate']:.3g} -> {c['exp6_rate']:.3g}")
    with open(out("baseline_reproduction.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote results/baseline_reproduction.md  (verdict: {verdict})")

    write_csv("baseline_envelope.csv", pooled)
    report["E0"] = {"verdict": verdict, "cells_compared": len(comparisons),
                    "bitwise_identical": len(exact),
                    "explained_by_contraction": len(explained),
                    "unexplained": len(unexplained),
                    "qualitative_changes": len(qualitative),
                    "comparisons": comparisons}
    return pooled


# ===========================================================================
# E1 - pilot and primitive selection
# ===========================================================================
def phase_pilot(report: dict) -> None:
    rows = IO.load_phase(RAW, "pilot", SMOKE)
    if not rows:
        print("  E1: no pilot rows")
        return

    shapes = sorted({shape_of(r["kappa"], r["precision"]) for r in rows})
    oracle = S.validate_log_q_evaluator(shapes)
    write_csv("log_q_evaluator_validation.csv", oracle)
    evaluator_ok = all(r["pass"] for r in oracle)

    scalar_rows = []
    for r in rows:
        path = r.get("raw_file")
        if not path or not os.path.exists(path):
            continue
        a = shape_of(r["kappa"], r["precision"])
        hdr, rec = IO.read_records(path, 1)
        log_r = rec["log_r"].astype(float)
        log_w = rec["log_w"].astype(float)
        finite = np.isfinite(log_r) & np.isfinite(log_w)
        n_fin = int(finite.sum())

        z = S.z_from_log_w(log_w[finite], a)
        gof = S.gof_exponential(z, "Z")
        # W is used only where it is numerically resolved, and "resolved" has to mean every
        # draw, not most of them.  At the smallest shapes a fraction of W underflows to
        # exactly zero; testing the survivors is testing a truncated sample, and the Beta
        # test then rejects the *truncation* rather than the sampler.  The criterion is
        # numerical and is decided before any statistic is computed.
        w = np.exp(log_w[finite])
        resolved = np.isfinite(w) & (w > 0) & (w < 1)
        w_test_run = bool(n_fin and resolved.all())
        gof += S.gof_beta(w, a, label="W") if w_test_run else []

        pvals = [g.pvalue for g in gof]
        rejects = S.holm(pvals)
        entry = {
            "method": r["method"], "precision": r["precision"], "stdlib": r["stdlib"],
            "tag": r["tag"], "kappa": r["kappa"], "shape_a": a, "seed": r["seed"],
            "n_attempted": r["n_attempted"], "n_finite_log": n_fin,
            "log_primitive_failures": r["log_primitive_failures"],
            "out_of_domain": r["out_of_domain"],
            "log_attempts": r["log_attempts"],
            "measured_acceptance": r["measured_acceptance"],
            "published_acceptance_eq2": r.get("published_acceptance_eq2"),
            "published_acceptance_area_ratio": r.get("published_acceptance_area_ratio"),
            "uniform_variates": r["uniform_variates"],
            "gamma_variates": r["gamma_variates"],
            "endpoint_redraws": r["endpoint_redraws"],
            "seconds": r["seconds"],
            "w_resolved_fraction": float(resolved.mean()) if n_fin else float("nan"),
            "w_test_run": w_test_run,
            "familywise_alpha": S.FAMILYWISE_ALPHA,
            "n_tests": len(gof),
            "any_rejected": bool(any(rejects)),
        }
        for g, rej in zip(gof, rejects):
            key = g.name.replace(":", "_").lower()
            entry[f"{key}_stat"] = g.statistic
            entry[f"{key}_p"] = g.pvalue
            entry[f"{key}_rejected"] = rej

        for p in S.QUANTILE_LEVELS:
            qi = S.quantile_with_interval(log_r[finite], p, bonferroni=len(S.QUANTILE_LEVELS))
            target = S.exact_log_r_quantile(p, a)
            entry[f"log_r_q{p}_estimate"] = qi["estimate"]
            entry[f"log_r_q{p}_target"] = target
            entry[f"log_r_q{p}_error"] = qi["estimate"] - target
            entry[f"log_r_q{p}_lo"] = qi["lo"]
            entry[f"log_r_q{p}_hi"] = qi["hi"]
            entry[f"log_r_q{p}_resolved"] = qi["resolved"]
            entry[f"log_r_q{p}_covers_target"] = bool(
                qi["resolved"] and qi["lo"] <= target <= qi["hi"])
        scalar_rows.append(entry)

    write_csv("scalar_validation.csv", scalar_rows)
    write_ecdf_source_data(rows)

    # ---- selection rule, applied before any later phase is read -----------------
    summary = {}
    for name in sorted({r["method"] for r in scalar_rows}):
        sub = [r for r in scalar_rows if r["method"] == name]
        in_domain = [r for r in sub if r["out_of_domain"] == 0]
        fails = sum(r["log_primitive_failures"] for r in sub)
        rejected = sum(1 for r in sub if r["any_rejected"])
        qcov = [r[f"log_r_q{p}_covers_target"] for r in sub for p in S.QUANTILE_LEVELS
                if r[f"log_r_q{p}_resolved"]]
        per_variate = []
        for r in sub:
            if r["log_attempts"] and r["seconds"]:
                per_variate.append(r["seconds"] / r["log_attempts"])
        summary[name] = {
            "configurations": len(sub),
            "in_domain_configurations": len(in_domain),
            "log_primitive_failures": fails,
            "configurations_with_a_rejection": rejected,
            "quantile_intervals_covering_target": int(sum(qcov)),
            "quantile_intervals_resolved": len(qcov),
            "median_seconds_per_log_gamma_attempt":
                float(np.median(per_variate)) if per_variate else float("nan"),
            "defined_at_kappa_0.5001": bool(
                all(r["n_finite_log"] == r["n_attempted"] for r in sub
                    if abs(r["kappa"] - 0.5001) < 1e-9)),
        }
        summary[name]["passes"] = bool(
            fails == 0 and rejected == 0 and summary[name]["defined_at_kappa_0.5001"]
            and (not qcov or sum(qcov) >= 0.9 * len(qcov)))

    passing = [k for k, v in summary.items() if v["passes"]]
    if not passing:
        selected, reason = None, "neither candidate passed the selection rule: NO-GO"
    elif len(passing) == 1:
        selected = passing[0]
        reason = f"{selected} is the only candidate that passed every selection criterion"
    else:
        by_cost = sorted(passing,
                         key=lambda k: summary[k]["median_seconds_per_log_gamma_attempt"])
        selected = by_cost[0]
        reason = (
            f"both candidates passed; {selected} has the lower median time per accepted "
            f"log-Gamma attempt "
            f"({summary[selected]['median_seconds_per_log_gamma_attempt']:.3g} s vs "
            f"{summary[by_cost[1]]['median_seconds_per_log_gamma_attempt']:.3g} s). "
            "LOG-ID is additionally rejection-free and defined for every shape, while "
            "LOG-PUB carries a published acceptance-rejection envelope whose stated domain "
            "is 0 < a < 1; the cost tie-break is recorded rather than overridden."
        )

    lines = [
        "# E1 - pilot decision",
        "",
        f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.",
        "",
        "## Frozen evaluator check, run before any sampling result was inspected",
        "",
        f"The two-piece `log q` evaluator (switch at `log W = {S.LOG_W_SWITCH:g}`) was compared",
        "against an `mpmath` incomplete beta at 60 decimal digits on a 41-point grid spanning",
        "10 natural-log units on each side of the switch, for every pilot shape.",
        "",
        f"**{'PASS' if evaluator_ok else 'FAIL'}** - worst absolute `log q` disagreement "
        f"{max(r['max_abs_log_q_error'] for r in oracle):.3g}, tolerance {S.LOG_Q_ORACLE_TOL:g}.",
        "Source data: `log_q_evaluator_validation.csv`.",
        "",
        "## Candidates",
        "",
        "| candidate | configs | in-domain | log failures | configs with a rejection | "
        "quantile intervals covering target | median s per log-Gamma attempt | defined at 0.5001 | passes |",
        "|---|---:|---:|---:|---:|---|---:|---|---|",
    ]
    for name, v in summary.items():
        lines.append(
            f"| {name} | {v['configurations']} | {v['in_domain_configurations']} | "
            f"{v['log_primitive_failures']} | {v['configurations_with_a_rejection']} | "
            f"{v['quantile_intervals_covering_target']}/{v['quantile_intervals_resolved']} | "
            f"{v['median_seconds_per_log_gamma_attempt']:.3g} | "
            f"{'yes' if v['defined_at_kappa_0.5001'] else 'no'} | "
            f"{'yes' if v['passes'] else 'no'} |")

    # The acceptance-rate observation, reported because it is a check on the transcription.
    pub = [r for r in scalar_rows if r["method"] == "LOG-PUB" and r["out_of_domain"] == 0]
    if pub:
        z_eq2, z_area = [], []
        for r in pub:
            n = r["n_attempted"]
            m = r["measured_acceptance"]
            for pred, acc in ((r["published_acceptance_eq2"], z_eq2),
                              (r["published_acceptance_area_ratio"], z_area)):
                if not pred or not (0 < pred < 1):
                    continue
                sd = pred * math.sqrt((1 - pred) / n)
                acc.append((m - pred) / sd if sd > 0 else 0.0)
        lines += [
            "",
            "## Transcription check on LOG-PUB",
            "",
            "Liu, Martin and Syring (2017) state the acceptance rate of their scheme as",
            "`r(alpha) = {1 + w(alpha)}^-1` (Eq. 2, p. 1771).  Their target `h_alpha` is",
            "normalized (`c^-1 = Gamma(alpha+1)`) while the envelope `eta_alpha` integrates to",
            "`c(1+w)`, so the ratio of areas - which is what an acceptance-rejection scheme's",
            "measured rate converges to - is `Gamma(alpha+1)/{1+w(alpha)}`.  The two agree to",
            "leading order in the small-shape regime the paper addresses, and the paper's own",
            "prose identifies `1/(1+w)` as the probability of *proposing* from the Exp(1) branch.",
            "",
            f"Measured against the area ratio: worst standardized deviation "
            f"{max(z_area, key=abs, default=float('nan')):+.2f}.",
            f"Measured against Eq. (2) as published: worst standardized deviation "
            f"{max(z_eq2, key=abs, default=float('nan')):+.2f}.",
            "",
            "The sampled law is unaffected either way; the implementation is used as transcribed.",
            "",
        ]

    lines += [
        "",
        "## Decision",
        "",
        f"**Selected primitive: {selected or 'none (NO-GO)'}.**",
        "",
        reason,
        "",
        "The selection was made from the pilot alone.  No later phase was read before this",
        "file was written, and `config/frozen.json` records the choice together with the",
        "SHA-256 of the primitive's source, which every later phase checks.",
    ]
    with open(out("pilot_decision.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote results/pilot_decision.md  (selected: {selected})")

    if selected:
        import hashlib
        src = os.path.join(HERE, "src",
                           "log_gamma_identity.H" if selected == "LOG-ID"
                           else "log_gamma_published.H")
        with open(src, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        # The freeze happens once.  Re-running the analysis must not rewrite this file:
        # it is part of the dependency set that preflight checks against HEAD, so a new
        # timestamp on every run would make the tree dirty and reopen G6 for no reason.
        frozen_path = os.path.join(CONFIG, "frozen.json")
        existing = {}
        if os.path.exists(frozen_path):
            try:
                existing = json.load(open(frozen_path, encoding="utf-8"))
            except json.JSONDecodeError:
                existing = {}
        if (existing.get("frozen") and existing.get("log_primitive") == selected
                and existing.get("source_sha256") == digest):
            print(f"  config/frozen.json already frozen ({selected}, {digest[:12]}) - "
                  "left unchanged")
        else:
            frozen = {
                "log_primitive": selected,
                "frozen": True,
                "source_file": os.path.relpath(src, HERE),
                "source_sha256": digest,
                "frozen_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "selection_summary": summary,
                "evaluator_validation_passed": evaluator_ok,
            }
            with open(frozen_path, "w", encoding="utf-8") as fh:
                json.dump(frozen, fh, indent=2)
                fh.write("\n")
            print(f"  wrote config/frozen.json  ({selected}, {digest[:12]})")

    report["E1"] = {"selected": selected, "reason": reason, "summary": summary,
                    "evaluator_validation_passed": evaluator_ok}


SFP1_KAPPAS = (0.501, 0.505, 0.51, 2.0)     # panels a and b
SFP1_OVERLAP_KAPPAS = (0.55, 2.0)           # panel c, where W and Z are both resolved
ECDF_POINTS = 300


def write_ecdf_source_data(rows: list[dict]) -> None:
    """ECDF residual curves for the scalar figure, reduced to plottable source rows.

    The figure script reads only committed CSVs, so the residual curves and their
    simultaneous bands are computed here and written out at a fixed 300-point resolution.
    The residual is plotted against the null CDF value, which puts both transforms on the
    same bounded axis and makes the two routes directly comparable.
    """
    want = set(SFP1_KAPPAS) | set(SFP1_OVERLAP_KAPPAS)
    by_cfg = collections.defaultdict(list)
    for r in rows:
        if any(abs(r["kappa"] - k) < 1e-9 for k in want) and r.get("raw_file") \
                and os.path.exists(r["raw_file"]):
            by_cfg[(r["method"], r["precision"], r["tag"], r["kappa"])].append(r["raw_file"])

    out_rows = []
    for (method, precision, tag, kappa), paths in sorted(by_cfg.items()):
        a = shape_of(kappa, precision)
        log_w = np.concatenate([IO.read_records(p, 1)[1]["log_w"].astype(float)
                                for p in sorted(paths)])
        log_w = log_w[np.isfinite(log_w)]
        if log_w.size < 100:
            continue
        for transform in ("Z", "W"):
            if transform == "Z":
                x = S.z_from_log_w(log_w, a)
                u = -np.expm1(-x)                       # Exp(1) CDF, stable near zero
                resolved = np.isfinite(u)
            else:
                if not any(abs(kappa - k) < 1e-9 for k in SFP1_OVERLAP_KAPPAS):
                    continue
                w = np.exp(log_w)
                ok = np.isfinite(w) & (w > 0) & (w < 1)
                if ok.sum() < 100:
                    continue
                u = np.full(log_w.shape, np.nan)
                u[ok] = special.betainc(a, S.BETA_B, w[ok])
                resolved = ok
            uu = np.sort(u[resolved])
            n = uu.size
            ecdf = np.arange(1, n + 1) / n
            resid = ecdf - uu
            # Simultaneous Kolmogorov band at the frozen familywise level.
            c = math.sqrt(-0.5 * math.log(S.FAMILYWISE_ALPHA / 2.0))
            band = c / math.sqrt(n)
            idx = np.unique(np.linspace(0, n - 1, min(ECDF_POINTS, n)).astype(int))
            for i in idx:
                # Z = -log I_W, so F_Z(z) = 1 - I_W: the two routes describe the same
                # deviation mirrored about the middle of the axis and negated.  Both
                # orientations are stored so the overlap panel can overlay them without
                # computing anything the source data does not already contain.
                out_rows.append({
                    "method": method, "precision": precision, "tag": tag, "kappa": kappa,
                    "shape_a": a, "transform": transform, "null_cdf": float(uu[i]),
                    "ecdf_residual": float(resid[i]),
                    "null_cdf_z_orientation":
                        float(uu[i]) if transform == "Z" else float(1.0 - uu[i]),
                    "ecdf_residual_z_orientation":
                        float(resid[i]) if transform == "Z" else float(-resid[i]),
                    "band_lo": -band, "band_hi": band,
                    "n": n, "resolved_fraction": float(np.mean(resolved)),
                })
    write_csv("scalar_ecdf.csv", out_rows)


# ===========================================================================
# E2 / E5 - failure envelope
# ===========================================================================
def avoidable_within(row: dict, near: int) -> bool:
    """For the log path, avoidable loss must be zero or inside the near-threshold band.

    A log-domain attempt cannot lose a representable draw to arithmetic; if it appears to,
    the draw sat within rounding distance of the representability threshold and the
    working-precision and reference evaluations of ``log|V_j|`` fell on opposite sides of
    it.  Anything beyond that band would be a real defect.
    """
    if row["method"] != "LOG":
        return True
    return row["n_avoidable"] <= near


def envelope_rows(rows: list[dict], phase: str) -> list[dict]:
    """Counts, rates and intervals per method and configuration, at both scopes.

    Seed-level rows are kept alongside the pooled ones: a pooled count does not replace
    replicate information, and the figures carry seed-level source rows behind every
    pooled marker.
    """
    outrows = []
    for scope, keys in (("pooled", ["layer", "method", "precision", "stdlib", "kappa"]),
                        ("seed", ["layer", "method", "precision", "stdlib", "kappa", "seed"])):
        for key, rs in group(rows, keys).items():
            d = dict(zip(keys, key))
            n = sum(r["n_attempted"] for r in rs)
            cats = {c: sum(r.get(f"cat_{c}", 0) for r in rs) for c in IO.CATEGORY_NAMES}
            avoidable = (cats["denominator_zero"] + cats["quotient_first_loss"] +
                         cats["split_form_loss"] + cats["log_primitive_failure"])
            honest = cats["honest_overflow"]
            failed = avoidable + honest
            row = {
                "phase": phase, "scope": scope, **d,
                "shape_a": shape_of(d["kappa"], d["precision"]),
                "n_seeds": len({r["seed"] for r in rs}), "n_attempted": n,
                "n_finite": cats["finite"], "n_failed": failed,
                "n_avoidable": avoidable, "n_honest_overflow": honest,
            }
            for c in IO.CATEGORY_NAMES:
                row[f"cat_{c}"] = cats[c]
            for label, k in (("failure", failed), ("avoidable", avoidable),
                             ("honest", honest)):
                est = S.rate_with_interval(k, n)
                row[f"{label}_rate"] = est["rate"]
                row[f"{label}_ci_lo"] = est["lo"]
                row[f"{label}_ci_hi"] = est["hi"]
                row[f"{label}_is_upper_bound"] = est["is_upper_bound"]
            row["accounting_ok"] = bool(cats["finite"] + avoidable + honest == n)

            # The honest floor is a property of the draw, not of a method: it is the number
            # of attempts whose intended final vector is not representable at all.  A method
            # can still disagree with it by a few draws, because a working-precision
            # evaluation of log|V_j| rounds the threshold differently from the reference.
            # Those draws are exactly the ones inside the audited near-limit band, so the
            # discrepancy is reported against that band rather than swept up as a rate.
            ref_nr = sum(r.get("ref_nonrepresentable", 0) for r in rs)
            near = sum(r.get("near_limit", 0) for r in rs)
            floor = S.rate_with_interval(ref_nr, n)
            row["honest_floor_count"] = ref_nr
            row["honest_floor_rate"] = floor["rate"]
            row["honest_floor_ci_lo"] = floor["lo"]
            row["honest_floor_ci_hi"] = floor["hi"]
            row["honest_floor_is_upper_bound"] = floor["is_upper_bound"]
            row["near_limit_band_count"] = near
            row["honest_minus_floor"] = honest - ref_nr
            row["within_boundary_band"] = bool(abs(honest - ref_nr) <= near
                                               and avoidable_within(row, near))
            row["x2_zero"] = sum(r.get("x2_zero", 0) for r in rs)
            row["x2_subnormal"] = sum(r.get("x2_subnormal", 0) for r in rs)
            row["rotation_recoverable"] = sum(r.get("rotation_recoverable", 0) for r in rs)
            row["mean_abs_rel_err_log_vs_split"] = float(np.mean(
                [r["mean_abs_rel_err_log_vs_split"] for r in rs
                 if isinstance(r.get("mean_abs_rel_err_log_vs_split"), float)] or [np.nan]))
            row["seconds"] = sum(r.get("seconds", 0) or 0 for r in rs)
            for st in ("near_limit", "method_disagreement", "subnormal_denominator",
                       "bulk_failure", "uniform_sample"):
                row[f"audit_{st}_total"] = sum(r.get(f"audit_{st}_total", 0) for r in rs)
                row[f"audit_{st}_audited"] = sum(r.get(f"audit_{st}_audited", 0) for r in rs)
            outrows.append(row)
    outrows.sort(key=lambda r: (r["scope"], r["layer"], r["precision"], r["kappa"],
                                r["method"], r.get("seed") or 0))
    return outrows


def phase_mechanism(report: dict) -> list[dict]:
    rows = IO.load_phase(RAW, "mechanism", SMOKE)
    if not rows:
        print("  E2: no mechanism rows")
        return []
    env = envelope_rows(rows, "E2")
    write_csv("failure_envelope.csv", env)

    bad = [r for r in env if r["scope"] == "seed" and not r["accounting_ok"]]
    pooled = [r for r in env if r["scope"] == "pooled" and r["layer"] == "paired"]

    # The honest floor is a property of the draw, not of a method: it is the same number in
    # every method's row for a given configuration, and no method can beat it.
    floor = {}
    for r in pooled:
        floor.setdefault((r["precision"], r["stdlib"], r["kappa"]), []).append(
            r["honest_floor_rate"])
    floor_consistent = all(max(v) - min(v) < 1e-12 for v in floor.values())
    outside_band = [r for r in pooled if not r["within_boundary_band"]]

    beats_floor = [r for r in pooled if r["method"] == "LOG" and r["n_avoidable"] > 0]
    log_double = [r for r in pooled if r["method"] == "LOG" and r["precision"] == "double"]
    report["E2"] = {
        "accounting_failures": len(bad),
        "honest_floor_consistent_across_methods": floor_consistent,
        "configurations_outside_boundary_band": len(outside_band),
        "log_configurations_with_any_avoidable_loss": len(beats_floor),
        "log_avoidable_loss_double_total": sum(r["n_avoidable"] for r in log_double),
        "configurations": len(pooled),
    }
    print(f"  E2: accounting failures {len(bad)}, configurations outside the boundary band "
          f"{len(outside_band)}, LOG avoidable loss in double "
          f"{sum(r['n_avoidable'] for r in log_double)}")
    return env


def phase_portability(report: dict) -> list[dict]:
    rows = IO.load_phase(RAW, "portability", SMOKE)
    if not rows:
        print("  E5: no portability rows")
        report["E5"] = {"status": "not run"}
        return []
    env = envelope_rows(rows, "E5")
    for r in env:
        r["arch"] = "arm64"
        r["layer"] = "native_scalar"
        r["available"] = True
        r["completed"] = True
        r["reason_not_run"] = "NA"

    # The end-to-end boundary subset.  Portability is a release gate, not a reason to
    # repeat the full three-dimensional grid, so only the uncapped boundary cases C0, C1
    # and C3 are carried across environments; the loader phase already ran them on every
    # available toolchain, and the rates are lifted here so that an absent platform cannot
    # be mistaken for a passing one.
    loader_rows = IO.load_phase(RAW, "loader", SMOKE)
    for key, rs in group([r for r in loader_rows
                          if r.get("case") in ("C0", "C1", "C3") and r.get("cap") is None
                          and r["method"] in ("SPLIT", "LOG")],
                         ["tag", "case", "method", "precision", "stdlib"]).items():
        tag, case, method, precision, stdlib = key
        n = sum(r["n_returned"] for r in rs)
        bad = sum(r.get("nonfinite_returned", 0) or 0 for r in rs) + \
            sum(r.get("honest_overflow", 0) or 0 for r in rs) + \
            sum(r.get("log_primitive_failures", 0) or 0 for r in rs)
        est = S.rate_with_interval(bad, n)
        env.append({
            "phase": "E5", "scope": "pooled", "layer": "end_to_end", "method": method,
            "precision": precision, "stdlib": stdlib, "tag": tag, "case": case,
            "kappa": rs[0]["kappa"], "shape_a": shape_of(rs[0]["kappa"], precision),
            "arch": "arm64", "available": True, "completed": True, "reason_not_run": "NA",
            "n_seeds": len({r["seed"] for r in rs}), "n_attempted": n, "n_failed": bad,
            "failure_rate": est["rate"], "failure_ci_lo": est["lo"],
            "failure_ci_hi": est["hi"], "failure_is_upper_bound": est["is_upper_bound"],
            "out_of_domain": sum(r.get("out_of_domain", 0) or 0 for r in rs),
        })

    # Architectures that were not available locally are recorded as missing, with the exact
    # command to run them elsewhere, so an absent platform can never be read as a pass.
    for stdlib in ("libc++", "libstdc++"):
        env.append({
            "phase": "E5", "scope": "pooled", "layer": "native_scalar", "method": "SPLIT",
            "precision": "NA", "stdlib": stdlib, "kappa": float("nan"),
            "arch": "x86_64", "available": False, "completed": False,
            "reason_not_run":
                "no x86_64 host available on the machine this experiment ran on. Rosetta "
                "or any other emulation is not accepted as a substitute, so this cell is "
                "left open rather than filled. To close it, on an x86_64 host at the "
                "archived commit: `cd experiments/exp6_low_kappa_stabilization && make "
                "selftest && make preflight && make portability && make loader`, then copy "
                "back raw/portability/ and raw/loader/ and rerun `make analyze`. The "
                + ("clang/libc++" if stdlib == "libc++" else "GCC/libstdc++") +
                " toolchain must be present; the build flags are recorded in "
                "raw/environment.json.",
            "n_attempted": 0, "n_failed": 0,
        })
    write_csv("portability.csv", env)
    report["E5"] = {
        "environments_completed": sorted({(r["arch"], r["stdlib"]) for r in env
                                          if r.get("completed")}),
        "environments_missing": sorted({(r["arch"], r["stdlib"]) for r in env
                                        if not r.get("completed")}),
        "gate_open": True,
    }
    return env


# ===========================================================================
# E3 - conditioning and tail consequence
# ===========================================================================
Q_BIN_EDGES = [1.0, 0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 0.0]


def phase_conditioning(report: dict) -> None:
    rows = IO.load_phase(RAW, "conditioning", SMOKE)
    if not rows:
        print("  E3: no conditioning rows")
        return

    files = {}
    for r in rows:
        p = r.get("raw_file")
        if p and os.path.exists(p):
            files[(r["tag"], r["precision"], r["kappa"], r["seed"])] = p

    bin_rows, tail_rows = [], []
    by_config = collections.defaultdict(list)
    for (tag, precision, kappa, seed), path in files.items():
        by_config[(tag, precision, kappa)].append((seed, path))

    for (tag, precision, kappa), seeds in sorted(by_config.items()):
        a = shape_of(kappa, precision)
        log_r_targets = {p: S.exact_log_r_quantile(p, a) for p in S.TAIL_LEVELS}

        per_seed = []
        for seed, path in sorted(seeds):
            _, rec = IO.read_records(path, 2)
            log_r = rec["log_r_ref"].astype(float)
            log_w = rec["log_w_ref"].astype(float)
            log_q = S.log_q_from_log_w(log_w, a)
            q = np.exp(log_q)
            succ = {
                "SPLIT": rec["cat_split"] == 0,
                "LOG": rec["cat_log"] == 0,
                "QF": rec["cat_qf"] == 0,
            }
            honest = (rec["flags"] & IO.FLAG_BITS["honest_overflow_ref"]) != 0
            per_seed.append({"seed": seed, "log_r": log_r, "q": q, "succ": succ,
                             "honest": honest, "n": log_r.size})

        n_total = sum(s["n"] for s in per_seed)

        # (a) success by target-tail bin
        for i in range(len(Q_BIN_EDGES) - 1):
            hi, lo = Q_BIN_EDGES[i], Q_BIN_EDGES[i + 1]
            for method in ("SPLIT", "LOG", "QF"):
                k_tot = n_tot = 0
                for s in per_seed:
                    m = (s["q"] <= hi) & (s["q"] > lo)
                    n_tot += int(m.sum())
                    k_tot += int((m & s["succ"][method]).sum())
                est = S.rate_with_interval(k_tot, n_tot)
                bin_rows.append({
                    "tag": tag, "precision": precision, "kappa": kappa, "shape_a": a,
                    "method": method, "bin_kind": "target_upper_tail_q",
                    "bin_hi": hi, "bin_lo": lo, "n_in_bin": n_tot, "n_success": k_tot,
                    "success_rate": est["rate"], "ci_lo": est["lo"], "ci_hi": est["hi"],
                    "is_upper_bound": est["is_upper_bound"],
                    "merged": bool(n_tot < 20),
                    "n_seeds": len(per_seed), "n_total": n_total,
                })

        # (b) success by log-radius bin, on a fixed grid so panels share axes
        all_log_r = np.concatenate([s["log_r"] for s in per_seed])
        finite_lr = all_log_r[np.isfinite(all_log_r)]
        if finite_lr.size:
            edges = np.quantile(finite_lr, np.linspace(0, 1, 11))
            edges[0] -= 1e-9
            edges[-1] += 1e-9
            for i in range(len(edges) - 1):
                for method in ("SPLIT", "LOG"):
                    k_tot = n_tot = 0
                    for s in per_seed:
                        m = (s["log_r"] > edges[i]) & (s["log_r"] <= edges[i + 1])
                        n_tot += int(m.sum())
                        k_tot += int((m & s["succ"][method]).sum())
                    est = S.rate_with_interval(k_tot, n_tot)
                    bin_rows.append({
                        "tag": tag, "precision": precision, "kappa": kappa, "shape_a": a,
                        "method": method, "bin_kind": "log_radius",
                        "bin_hi": float(edges[i + 1]), "bin_lo": float(edges[i]),
                        "n_in_bin": n_tot, "n_success": k_tot, "success_rate": est["rate"],
                        "ci_lo": est["lo"], "ci_hi": est["hi"],
                        "is_upper_bound": est["is_upper_bound"], "merged": bool(n_tot < 20),
                        "n_seeds": len(per_seed), "n_total": n_total,
                    })

        # (c) tail ratios and quantile errors
        for method in ("SPLIT", "LOG"):
            n_succ = sum(int(s["succ"][method].sum()) for s in per_seed)
            fail = S.rate_with_interval(n_total - n_succ, n_total)
            for p in S.TAIL_LEVELS:
                thr = log_r_targets[p]
                # Counts per seed, not arrays: the bootstrap of a pooled proportion is done
                # on the counts, which is exact and finishes.
                a_counts = [int(((s["log_r"] > thr) & s["succ"][method]).sum())
                            for s in per_seed]
                a_sizes = [s["n"] for s in per_seed]
                c_counts = [int((s["log_r"][s["succ"][method]] > thr).sum())
                            for s in per_seed]
                c_sizes = [int(s["succ"][method].sum()) for s in per_seed]
                boot_a = S.seed_stratified_bootstrap_indicator(a_counts, a_sizes,
                                                               scale=1 - p)
                boot_c = S.seed_stratified_bootstrap_indicator(c_counts, c_sizes,
                                                               scale=1 - p)
                surv = np.concatenate([s["log_r"][s["succ"][method]] for s in per_seed]) \
                    if n_succ else np.array([])
                qi = S.quantile_with_interval(surv, p, bonferroni=len(S.TAIL_LEVELS))
                tail_rows.append({
                    "tag": tag, "precision": precision, "kappa": kappa, "shape_a": a,
                    "method": method, "p": p, "n_attempted": n_total, "n_success": n_succ,
                    "target_log_r_quantile": thr,
                    "rho_attempt": boot_a["estimate"], "rho_attempt_lo": boot_a["lo"],
                    "rho_attempt_hi": boot_a["hi"],
                    "rho_cond": boot_c["estimate"], "rho_cond_lo": boot_c["lo"],
                    "rho_cond_hi": boot_c["hi"],
                    "conditional_log_r_quantile": qi["estimate"],
                    "conditional_log_r_quantile_lo": qi["lo"],
                    "conditional_log_r_quantile_hi": qi["hi"],
                    "conditional_log_r_quantile_error":
                        qi["estimate"] - thr if qi["resolved"] else float("nan"),
                    "quantile_resolved": qi["resolved"],
                    "unconditional_failure_probability": fail["rate"],
                    "failure_ci_lo": fail["lo"], "failure_ci_hi": fail["hi"],
                    "failure_is_upper_bound": fail["is_upper_bound"],
                    "bootstrap_resamples": S.BOOTSTRAP_RESAMPLES,
                    "interval_conf": S.TAIL_INTERVAL_CONF,
                    "n_seeds": len(per_seed),
                })

    write_csv("conditioning_bins.csv", bin_rows)
    write_csv("tail_metrics.csv", tail_rows)

    monotone = []
    for key, rs in group([r for r in bin_rows if r["bin_kind"] == "target_upper_tail_q"
                          and r["method"] == "SPLIT"],
                         ["tag", "precision", "kappa"]).items():
        rs = [r for r in rs if r["n_in_bin"] >= 20]
        if len(rs) >= 3:
            rs.sort(key=lambda r: -r["bin_hi"])
            monotone.append({"config": key,
                             "success_first_bin": rs[0]["success_rate"],
                             "success_last_bin": rs[-1]["success_rate"]})
    report["E3"] = {"bins": len(bin_rows), "tail_rows": len(tail_rows),
                    "state_dependence": monotone}
    print(f"  E3: {len(bin_rows)} bin rows, {len(tail_rows)} tail rows")


# ===========================================================================
# E4 - complete loader validation
# ===========================================================================
def capped_uniform_transform(log_r, n_hat, kappa, cap, a):
    """Probability-integral transform of the radius under the declared bounded law.

    Conditional on the direction, the capped target is the uncapped radial law truncated at
    ``R_max(n) = cap / (sqrt(kappa) max_j |n_j|)``.  So
    ``U = F_R(R) / F_R(R_max(n))`` is Uniform(0,1) and independent of the direction - an
    exact test of the conditional law that does not assume the direction stays isotropic,
    which under a component-wise box it does not.

    ``F_R`` is evaluated as ``I_{1-w}(3/2, a)`` rather than as ``1 - I_w(a, 3/2)``: the
    accepted region sits where ``1 - F_R`` is close to one, and forming the complement by
    subtraction there loses every significant digit.
    """
    nmax = np.max(np.abs(n_hat), axis=1)
    ok = np.isfinite(log_r) & np.isfinite(nmax) & (nmax > 0)
    u = np.full(log_r.shape, np.nan)
    log_r_max = np.log(cap) - 0.5 * np.log(kappa) - np.log(nmax[ok])
    one_minus_w = special.expit(2.0 * log_r[ok])
    one_minus_w_max = special.expit(2.0 * log_r_max)
    f_r = special.betainc(S.BETA_B, a, one_minus_w)
    f_max = special.betainc(S.BETA_B, a, one_minus_w_max)
    with np.errstate(divide="ignore", invalid="ignore"):
        u[ok] = np.where(f_max > 0, f_r / f_max, np.nan)
    return u, nmax


def inject_radius_direction_coupling(log_r, cos_theta, fraction=0.5, seed=60002):
    """Couple radius and direction while leaving both one-dimensional marginals untouched.

    Within a random subset, the radii and the direction cosines are each sorted and then
    re-paired.  Every value that was in the sample is still in the sample exactly once, so
    the marginals are bit-identical; only the pairing changes.  A dependence diagnostic that
    cannot see this is not powerful enough to certify independence.
    """
    rng = np.random.default_rng(seed)
    lr = log_r.copy()
    ct = cos_theta.copy()
    n = lr.size
    idx = rng.choice(n, size=int(fraction * n), replace=False)
    lr[idx] = np.sort(lr[idx])
    ct[idx] = np.sort(ct[idx])
    return lr, ct


def loader_tests(hdr, rec, case, method, counters):
    """Every predeclared test for one (case, method), before multiplicity correction."""
    kappa = hdr["kappa"]
    tperp, tpar = hdr["theta_perp"], hdr["theta_par"]
    ub = np.array(hdr["ub"])
    cap = hdr["cap"]
    capped = math.isfinite(cap)
    a = shape_of(kappa, hdr["precision"])

    v = rec["v"]
    status = (rec["status"] & 0xFF).astype(int)
    finite = np.all(np.isfinite(v), axis=1)
    n_total = v.shape[0]
    n_finite = int(finite.sum())

    log_r, n_hat = IO.recover_radius_direction(v, kappa, tperp, tpar, ub)
    good = finite & np.isfinite(log_r)

    tests = []
    aux = {"n_total": n_total, "n_finite": n_finite,
           "nonfinite_fraction": 1.0 - n_finite / n_total if n_total else float("nan")}

    if good.sum() >= 100:
        lr = log_r[good]
        nh = n_hat[good]
        cos_theta = nh[:, 2]
        phi = np.arctan2(nh[:, 1], nh[:, 0])

        if not capped:
            log_w = S.log_w_from_log_r(lr)
            z = S.z_from_log_w(log_w, a)
            tests += [("radial", g) for g in S.gof_exponential(z, "Z")]
            w = np.exp(log_w)
            res = np.isfinite(w) & (w > 0) & (w < 1)
            if res.all():
                tests += [("radial", g) for g in S.gof_beta(w, a, label="W")]
            aux["w_resolved_fraction"] = float(res.mean())
            aux["w_test_run"] = float(bool(res.all()))
            tests += [("direction", g) for g in S.gof_uniform(cos_theta, -1, 1, "cosTheta")]
            tests += [("direction", g) for g in S.gof_uniform(phi, -math.pi, math.pi, "Phi")]
            tests.append(("independence", S.binned_independence_test(lr, cos_theta)))
        else:
            u, nmax = capped_uniform_transform(lr, nh, kappa, cap, a)
            uok = np.isfinite(u) & (u >= 0) & (u <= 1)
            tests += [("cap_law", g) for g in S.gof_uniform(u[uok], 0, 1, "U_cap")]
            tests.append(("cap_law", S.binned_independence_test(u[uok], nmax[uok])))
            # Every returned sample must satisfy the declared box.
            inside = np.all(np.abs(nh) * np.exp(lr)[:, None] * math.sqrt(kappa) <= cap * (1 + 1e-9),
                            axis=1)
            aux["cap_violations"] = int((~inside).sum())
            aux["u_resolved_fraction"] = float(uok.mean())

        # Anisotropy and frame, from the returned vectors and the field direction alone.
        q = IO.field_basis(ub)
        with np.errstate(over="ignore", invalid="ignore"):
            v_par = v[good] @ q[:, 2]
            v_perp = v[good] @ q[:, 0]
        fin = np.isfinite(v_par) & np.isfinite(v_perp) & (v_par != 0) & (v_perp != 0)
        aux["anisotropy_dropped_nonfinite"] = int((~fin).sum())
        if fin.sum() > 100:
            target = tpar / tperp
            mr = S.mad_ratio_interval(v_par[fin], v_perp[fin])
            aux["mad_ratio"] = mr["ratio"]
            aux["mad_ratio_lo"] = mr["lo"]
            aux["mad_ratio_hi"] = mr["hi"]
            aux["mad_ratio_bootstrap_subsample"] = mr["subsample"]
            aux["mad_ratio_target"] = target
            aux["mad_ratio_relative_error"] = (mr["ratio"] / target - 1.0 if target
                                               else float("nan"))
            aux["mad_ratio_covers_target"] = bool(
                np.isfinite(mr["lo"]) and mr["lo"] <= target <= mr["hi"])
            # Frame invariance: after dividing out the declared thermal speeds, the parallel
            # and one perpendicular projection must be identically distributed, because the
            # direction is isotropic before the anisotropic scaling is applied.
            lp = np.log(np.abs(v_par[fin])) - math.log(tpar)
            lq = np.log(np.abs(v_perp[fin])) - math.log(tperp)
            ks = stats.ks_2samp(lp, lq)
            tests.append(("frame", S.GofResult("frame:KS2", float(ks.statistic),
                                               float(ks.pvalue), int(fin.sum()))))

    # Finiteness and status consistency against the probe's own counters.
    if counters:
        # A log primitive used outside its published domain is not evidence about that
        # method.  C0 is the benign reference case at kappa = 2, where LOG-PUB's stated
        # domain 0 < a < 1 does not apply and every draw takes the documented fallback.
        aux["out_of_domain"] = counters.get("out_of_domain", 0)
        aux["out_of_domain_fraction"] = (
            counters.get("out_of_domain", 0) / n_total if n_total else float("nan"))
        declared = counters.get("nonfinite_returned")
        if declared is None:
            declared = counters.get("honest_overflow", 0) + \
                counters.get("log_primitive_failures", 0)
        aux["counter_nonfinite"] = declared
        aux["status_nonfinite"] = int((status != 0).sum())
        aux["finiteness_consistent"] = bool(
            abs(int((status != 0).sum()) - (n_total - n_finite)) <=
            counters.get("log_primitive_failures", 0))
    return tests, aux


def phase_loader(report: dict) -> None:
    rows = IO.load_phase(RAW, "loader", SMOKE)
    if not rows:
        print("  E4: no loader rows")
        return

    files = collections.defaultdict(list)
    counters = collections.defaultdict(dict)
    for r in rows:
        p = r.get("raw_file")
        if p and os.path.exists(p):
            files[(r["tag"], r["case"], r["method"])].append((r["seed"], p))
            for k, v in r.items():
                if isinstance(v, (int, float)) and k not in ("kappa", "shape_a", "seed"):
                    counters[(r["tag"], r["case"], r["method"])][k] = \
                        counters[(r["tag"], r["case"], r["method"])].get(k, 0) + v

    val_rows = []
    matrix = {}
    nc_store = {}

    for key in sorted(files):
        tag, case, method = key
        recs, hdr = [], None
        for seed, path in sorted(files[key]):
            h, a = IO.read_records(path, 3)
            hdr = h
            recs.append(a)
        rec = np.concatenate(recs)
        tests, aux = loader_tests(hdr, rec, case, method, counters[key])

        pvals = [g.pvalue for _, g in tests]
        rejects = S.holm(pvals)
        for (family, g), rej in zip(tests, rejects):
            val_rows.append({
                "tag": tag, "case": case, "method": method,
                "precision": hdr["precision"], "kappa": hdr["kappa"],
                "theta_par_over_perp": hdr["theta_par"] / hdr["theta_perp"],
                "cap": hdr["cap"], "family": family, "test": g.name,
                "statistic": g.statistic, "pvalue": g.pvalue, "n": g.n,
                "familywise_alpha": S.FAMILYWISE_ALPHA,
                "multiplicity_family_size": len(tests),
                "corrected_decision": "reject" if rej else "pass",
                "n_total": aux["n_total"], "n_finite": aux["n_finite"],
                "source_files": ";".join(os.path.relpath(p, HERE)
                                         for _, p in sorted(files[key])),
            })
        for k, v in aux.items():
            val_rows.append({
                "tag": tag, "case": case, "method": method,
                "precision": hdr["precision"], "kappa": hdr["kappa"],
                "theta_par_over_perp": hdr["theta_par"] / hdr["theta_perp"],
                "cap": hdr["cap"], "family": "diagnostic", "test": k,
                "statistic": v if isinstance(v, (int, float)) else float("nan"),
                "pvalue": float("nan"), "n": aux["n_total"],
                "familywise_alpha": S.FAMILYWISE_ALPHA,
                "multiplicity_family_size": len(tests),
                "corrected_decision": "NA",
                "n_total": aux["n_total"], "n_finite": aux["n_finite"],
                "source_files": ";".join(os.path.relpath(p, HERE)
                                         for _, p in sorted(files[key])),
            })

        cell = {}
        for family in ("radial", "direction", "independence", "frame", "cap_law"):
            fam = [(g, rej) for (f, g), rej in zip(tests, rejects) if f == family]
            if not fam:
                cell[family] = "na"
            elif any(rej for _, rej in fam):
                cell[family] = "fail"
            else:
                cell[family] = "pass"
        # The MAD ratio is an estimate, not a hypothesis test: the cell passes when its
        # bootstrap interval covers the declared thermal-speed ratio.
        cell["anisotropy"] = ("unresolved" if "mad_ratio_covers_target" not in aux
                              else "pass" if aux["mad_ratio_covers_target"] else "fail")
        cell["finiteness"] = "pass" if aux.get("finiteness_consistent", False) else (
            "unresolved" if "finiteness_consistent" not in aux else "fail")
        cell["n"] = aux["n_total"]
        cell["n_finite"] = aux["n_finite"]
        matrix[(tag, case, method)] = cell

        if case in ("C1", "C5") and method == "LOG":
            nc_store[(tag, case)] = (hdr, rec)

    # ---- negative controls -------------------------------------------------
    for (tag, case), (hdr, rec) in sorted(nc_store.items()):
        kappa, tperp, tpar, ub = hdr["kappa"], hdr["theta_perp"], hdr["theta_par"], hdr["ub"]
        a = shape_of(kappa, hdr["precision"])
        v = rec["v"]
        log_r, n_hat = IO.recover_radius_direction(v, kappa, tperp, tpar, ub)
        good = np.all(np.isfinite(v), axis=1) & np.isfinite(log_r)
        if good.sum() < 1000:
            continue

        if case == "C1":
            lr, ct = inject_radius_direction_coupling(log_r[good], n_hat[good][:, 2])
            g = S.binned_independence_test(lr, ct)
            detected = bool(np.isfinite(g.pvalue) and g.pvalue < S.FAMILYWISE_ALPHA)
            matrix[(tag, "NC1", "LOG")] = {
                "radial": "na", "direction": "na",
                "independence": "pass" if detected else "fail",
                "frame": "na", "cap_law": "na", "anisotropy": "na",
                "finiteness": "na", "n": int(good.sum()), "n_finite": int(good.sum())}
            val_rows.append({
                "tag": tag, "case": "NC1", "method": "LOG", "precision": hdr["precision"],
                "kappa": kappa, "theta_par_over_perp": tpar / tperp, "cap": hdr["cap"],
                "family": "negative_control", "test": "NC1:injected_radius_direction_coupling",
                "statistic": g.statistic, "pvalue": g.pvalue, "n": g.n,
                "familywise_alpha": S.FAMILYWISE_ALPHA, "multiplicity_family_size": 1,
                "corrected_decision": "detected" if detected else "MISSED",
                "n_total": int(good.sum()), "n_finite": int(good.sum()),
                "source_files": "derived from C1 LOG samples",
            })

        if case == "C5":
            # The capped sample tested against the *uncapped* target must be rejected.
            log_w = S.log_w_from_log_r(log_r[good])
            z = S.z_from_log_w(log_w, a)
            gs = S.gof_exponential(z, "Z")
            detected = bool(gs and min(g.pvalue for g in gs) < S.FAMILYWISE_ALPHA)
            matrix[(tag, "NC2", "LOG")] = {
                "radial": "pass" if detected else "fail", "direction": "na",
                "independence": "na", "frame": "na", "cap_law": "na", "anisotropy": "na",
                "finiteness": "na", "n": int(good.sum()), "n_finite": int(good.sum())}
            for g in gs:
                val_rows.append({
                    "tag": tag, "case": "NC2", "method": "LOG", "precision": hdr["precision"],
                    "kappa": kappa, "theta_par_over_perp": tpar / tperp, "cap": hdr["cap"],
                    "family": "negative_control",
                    "test": f"NC2:capped_sample_vs_uncapped_target:{g.name}",
                    "statistic": g.statistic, "pvalue": g.pvalue, "n": g.n,
                    "familywise_alpha": S.FAMILYWISE_ALPHA, "multiplicity_family_size": len(gs),
                    "corrected_decision": "detected" if detected else "MISSED",
                    "n_total": int(good.sum()), "n_finite": int(good.sum()),
                    "source_files": "derived from C5 LOG samples",
                })

    write_csv("loader_validation.csv", val_rows)

    cols = ["radial", "direction", "independence", "anisotropy", "frame", "finiteness",
            "cap_law"]
    matrix_rows = []
    for (tag, case, method), cell in sorted(matrix.items()):
        row = {"tag": tag, "case": case, "method": method, "n": cell["n"],
               "n_finite": cell["n_finite"]}
        for c in cols:
            row[c] = cell.get(c, "na")
        matrix_rows.append(row)
    write_csv("validation_matrix.csv", matrix_rows)
    lines = ["# E4 - complete-loader validation matrix", "",
             f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.", "",
             "Decisions are Holm-corrected within each configuration at a familywise "
             f"alpha of {S.FAMILYWISE_ALPHA}.  `pass` means the frozen test did not reject; "
             "it is consistency at the tested sample size, not proof of exactness.  NC1 and "
             "NC2 pass only when the injected defect *is* detected.", "",
             "| tag | case | method | n | n finite | " + " | ".join(cols) + " |",
             "|---|---|---|---:|---:|" + "|".join(["---"] * len(cols)) + "|"]
    for (tag, case, method), cell in sorted(matrix.items()):
        lines.append(f"| {tag} | {case} | {method} | {cell['n']} | {cell['n_finite']} | " +
                     " | ".join(cell.get(c, "na") for c in cols) + " |")
    failures = [k for k, c in matrix.items() if any(c.get(x) == "fail" for x in cols)]
    # A rejection of the *released* path's radial law on its finite survivors is the
    # measured consequence of discarding failures, not a defect in the loader under test.
    # It is separated here so that the fidelity gate is about the candidate, and reported
    # in its own right because it is the point of the experiment.
    candidate_failures = [k for k in failures if k[2] in ("LOG",)]
    conditional_rejections = [
        k for k in failures
        if k[2] in ("SPLIT", "SPLIT_released") and matrix[k].get("radial") == "fail"]
    if conditional_rejections:
        lines += [
            "", "## Released-path radial rejections on the finite survivors", "",
            "These are results, not defects. Where the released formation discards a "
            "non-negligible share of the intended draws, the draws it keeps are no longer "
            "distributed as the target: the radial test rejects on the survivor subset. The "
            "log path passes the same test on the same configuration. The size of the effect "
            "is in `tail_metrics.csv`.", ""]
        for k in conditional_rejections:
            c = matrix[k]
            lost = c["n"] - c["n_finite"]
            lines.append(f"- {k[0]} / {k[1]} / {k[2]}: {lost} of {c['n']} draws discarded "
                         f"({lost / c['n']:.2e} per attempt); radial law rejected on the "
                         "remaining survivors")
    other = [k for k in failures if k not in conditional_rejections]
    if other:
        lines += ["", "## Cells that failed", ""]
        for k in other:
            lines.append(f"- {k}: see `loader_validation.csv` rows with the matching "
                         "tag/case/method")
    with open(out("validation_matrix.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote results/validation_matrix.md  ({len(matrix)} cells, "
          f"{len(failures)} with a failure)")

    nc_ok = all(matrix.get((t, nc, "LOG"), {}).get(
        "independence" if nc == "NC1" else "radial") == "pass"
        for t, nc in [(k[0], k[1]) for k in matrix if k[1] in ("NC1", "NC2")])
    report["E4"] = {"cells": len(matrix), "cells_with_failure": len(failures),
                    "candidate_cells_with_failure": len(candidate_failures),
                    "released_path_conditional_rejections": [
                        {"tag": k[0], "case": k[1], "method": k[2],
                         "n": matrix[k]["n"], "n_finite": matrix[k]["n_finite"]}
                        for k in conditional_rejections],
                    "negative_controls_detected": bool(nc_ok),
                    "matrix": {"/".join(k): v for k, v in matrix.items()}}


# ===========================================================================
# E6 - performance
# ===========================================================================
def phase_performance(report: dict) -> list[dict]:
    rows = IO.load_phase(RAW, "performance", SMOKE)
    if not rows:
        print("  E6: no performance rows")
        return []

    blocks = [r for r in rows if r["layer"] != "reproducibility"]
    repro = [r for r in rows if r["layer"] == "reproducibility"]

    perf_rows = []
    by_case = group(blocks, ["tag", "case", "precision", "kappa", "cap", "method"])
    for key, rs in sorted(by_case.items(), key=lambda kv: tuple(str(x) for x in kv[0])):
        tag, case, precision, kappa, cap, method = key
        per_attempt = np.array([r["seconds"] / r["n_attempted"] for r in rs if r["n_attempted"]])
        per_return = np.array([r["seconds"] / r["n_returned"] for r in rs if r["n_returned"]])
        for r in rs:                      # block-level rows are preserved, not summarized away
            perf_rows.append({
                "scope": "block", "tag": tag, "case": case, "precision": precision,
                "kappa": kappa, "cap": cap, "method": method, "block": r["block"],
                "block_order": r["block_order"], "seconds": r["seconds"],
                "n_attempted": r["n_attempted"], "n_returned": r["n_returned"],
                "n_finite": r["n_finite"], "cap_rejects": r["cap_rejects"],
                "gamma_variates": r["gamma_variates"],
                "uniform_variates": r["uniform_variates"],
                "engine_calls": r["engine_calls"],
                "seconds_per_attempt": r["seconds"] / r["n_attempted"] if r["n_attempted"] else float("nan"),
                "seconds_per_return": r["seconds"] / r["n_returned"] if r["n_returned"] else float("nan"),
                "attempts_per_return": r["n_attempted"] / r["n_returned"] if r["n_returned"] else float("nan"),
                "log_primitive": r.get("log_primitive"),
            })
        perf_rows.append({
            "scope": "summary", "tag": tag, "case": case, "precision": precision,
            "kappa": kappa, "cap": cap, "method": method, "block": None,
            "n_blocks": len(rs),
            "median_seconds_per_attempt": float(np.median(per_attempt)) if per_attempt.size else float("nan"),
            "iqr_seconds_per_attempt": float(np.subtract(*np.percentile(per_attempt, [75, 25]))) if per_attempt.size else float("nan"),
            "median_seconds_per_return": float(np.median(per_return)) if per_return.size else float("nan"),
            "iqr_seconds_per_return": float(np.subtract(*np.percentile(per_return, [75, 25]))) if per_return.size else float("nan"),
            "total_attempts": sum(r["n_attempted"] for r in rs),
            "total_returned": sum(r["n_returned"] for r in rs),
            "gamma_variates_per_return": sum(r["gamma_variates"] for r in rs) / max(sum(r["n_returned"] for r in rs), 1),
            "uniform_variates_per_return": sum(r["uniform_variates"] for r in rs) / max(sum(r["n_returned"] for r in rs), 1),
            "engine_calls_per_return": sum(r["engine_calls"] for r in rs) / max(sum(r["n_returned"] for r in rs), 1),
            "cap_rejects_per_return": sum(r["cap_rejects"] for r in rs) / max(sum(r["n_returned"] for r in rs), 1),
            "log_primitive": rs[0].get("log_primitive"),
        })

    # Paired blockwise LOG/SPLIT ratio: the two methods ran the same block indices with the
    # same seeds, so the ratio is paired block by block rather than pooled.
    ratios = []
    for key, rs in group(blocks, ["tag", "case", "precision", "kappa", "cap"]).items():
        tag, case, precision, kappa, cap = key
        by_method = {}
        for r in rs:
            by_method.setdefault(r["method"], {})[r["block"]] = r
        if "LOG" not in by_method or "SPLIT" not in by_method:
            continue
        common = sorted(set(by_method["LOG"]) & set(by_method["SPLIT"]))
        if not common:
            continue
        rat = np.array([
            (by_method["LOG"][b]["seconds"] / by_method["LOG"][b]["n_returned"]) /
            (by_method["SPLIT"][b]["seconds"] / by_method["SPLIT"][b]["n_returned"])
            for b in common
            if by_method["LOG"][b]["n_returned"] and by_method["SPLIT"][b]["n_returned"]])
        if not rat.size:
            continue
        boot = S.seed_stratified_bootstrap([rat], lambda arrs: float(np.median(arrs[0])),
                                           conf=0.95)
        ratios.append({
            "scope": "ratio", "tag": tag, "case": case, "precision": precision,
            "kappa": kappa, "cap": cap, "method": "LOG/SPLIT", "n_blocks": rat.size,
            "median_time_per_return_ratio": float(np.median(rat)),
            "ratio_ci_lo": boot["lo"], "ratio_ci_hi": boot["hi"],
            "iqr_ratio": float(np.subtract(*np.percentile(rat, [75, 25]))),
        })
    perf_rows.extend(ratios)

    for r in repro:
        perf_rows.append({
            "scope": "reproducibility", "tag": r["tag"], "case": r["case"],
            "precision": r["precision"], "kappa": r["kappa"], "cap": r["cap"],
            "method": r["method"], "repeat_identical": r["repeat_identical"],
            "attempts_1": r["attempts_1"], "attempts_2": r["attempts_2"],
            "engine_calls_1": r["engine_calls_1"], "engine_calls_2": r["engine_calls_2"],
        })

    fieldnames = sorted({k for r in perf_rows for k in r})
    fieldnames = (["scope", "tag", "case", "method", "precision", "kappa", "cap"] +
                  [f for f in fieldnames if f not in
                   ("scope", "tag", "case", "method", "precision", "kappa", "cap")])
    write_csv("performance.csv", perf_rows, fieldnames)

    repro_ok = all(r.get("repeat_identical") is True for r in repro) if repro else None
    report["E6"] = {
        "ratios": ratios,
        "same_seed_reproducible": repro_ok,
        "worst_time_ratio": max((r["median_time_per_return_ratio"] for r in ratios),
                                default=float("nan")),
    }
    return perf_rows


# ===========================================================================
# Tables
# ===========================================================================
def latex_escape(s: str) -> str:
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def write_tex_table(name: str, caption: str, label: str, header: list[str],
                    body: list[list[str]], colspec: str | None = None) -> None:
    colspec = colspec or ("l" * len(header))
    lines = [r"% Generated by experiments/exp6_low_kappa_stabilization/analyze.py",
             r"% Do not edit by hand; regenerate with `make analyze`.",
             r"\begin{table}[t]", r"\centering",
             rf"\caption{{{caption}}}", rf"\label{{{label}}}",
             rf"\begin{{tabular}}{{{colspec}}}", r"\hline"]
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\hline")
    for row in body:
        lines.append(" & ".join(str(c) for c in row) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    with open(out(name), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote results/{name}")


def fmt_rate(rate, lo, hi, upper) -> str:
    if rate is None or (isinstance(rate, float) and math.isnan(rate)):
        return "--"
    if hi is None or (isinstance(hi, float) and math.isnan(hi)):
        return rf"${rate:.2e}$"
    if upper:
        return rf"$<{hi:.1e}$"
    return rf"${rate:.2e}$"


def build_tables(report: dict, env_rows, perf_rows, tail_rows_path: str) -> None:
    selected = report.get("E1", {}).get("selected") or "LOG"
    native = [r for r in env_rows if r["scope"] == "pooled" and r["layer"] == "native"]

    def no_failure_range(method, precision):
        """Smallest kappa on the ladder at and above which no failure was observed, with the
        upper bound that the sample size actually supports."""
        rs = sorted([r for r in native
                     if r["method"] == method and r["precision"] == precision],
                    key=lambda r: r["kappa"])
        best, bound = None, float("nan")
        for r in reversed(rs):
            if r["n_failed"] == 0:
                best, bound = r["kappa"], r["failure_ci_hi"]
            else:
                break
        return best, bound

    def rate_at(method, precision, kappa):
        rs = [r for r in native if r["method"] == method and r["precision"] == precision
              and abs(r["kappa"] - kappa) < 1e-9]
        if not rs:
            return None
        r = max(rs, key=lambda x: x["failure_rate"])   # worst across standard libraries
        return r

    tails = {}
    if os.path.exists(tail_rows_path):
        with open(tail_rows_path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if row["p"] == "0.999":
                    tails[(row["method"], row["precision"], round(float(row["kappa"]), 6))] = row

    ratio = {}
    for r in perf_rows:
        if r.get("scope") == "ratio":
            # A null cap is the uncapped benchmark; JSON has no infinity literal.
            ratio[(r["precision"], round(float(r["kappa"]), 6), r["cap"])] = r

    body = []
    for method, semantics, stream in (
        ("SPLIT", "returns non-finite components; the cap silently redraws them",
         "released stream"),
        ("LOG", "reports honest overflow; the cap rejects in the log domain",
         "different stream"),
    ):
        drange, dbound = no_failure_range(method, "double")
        frange, fbound = no_failure_range(method, "float")
        rd = rate_at(method, "double", 0.51)
        rf_ = rate_at(method, "float", 0.55)
        td = tails.get((method, "double", 0.51))
        tf = tails.get((method, "float", 0.55))
        rr = ratio.get(("double", 0.51, None))
        body.append([
            latex_escape(method if method == "SPLIT" else f"{method} ({selected})"),
            latex_escape(semantics),
            (rf"$\kappa\ge{drange:g}$ ($<{dbound:.1e}$)" if drange else "--") + r" / " +
            (rf"$\kappa\ge{frange:g}$ ($<{fbound:.1e}$)" if frange else "--"),
            (fmt_rate(rd["failure_rate"], rd["failure_ci_lo"], rd["failure_ci_hi"],
                      rd["failure_is_upper_bound"]) if rd else "--") + " / " +
            (fmt_rate(rf_["failure_rate"], rf_["failure_ci_lo"], rf_["failure_ci_hi"],
                      rf_["failure_is_upper_bound"]) if rf_ else "--"),
            (f"{float(td['rho_attempt']):.3f}/{float(td['rho_cond']):.3f}" if td else "--")
            + " / " +
            (f"{float(tf['rho_attempt']):.3f}/{float(tf['rho_cond']):.3f}" if tf else "--"),
            (f"{rr['median_time_per_return_ratio']:.2f}" if rr and method == "LOG"
             else "1.00" if method == "SPLIT" else "--"),
            latex_escape(stream),
        ])

    write_tex_table(
        "table_fp1.tex",
        "Operating range, overflow semantics and cost of the released Gamma-ratio loader "
        "and the log-domain path, in double / single precision. Ranges are empirical under "
        "the tested protocol: an entry $\\kappa\\ge k$ means no failure was observed at or "
        "above $k$, with the stated one-sided 95\\% upper bound on the per-draw rate. "
        "Tail-retention ratios are per-attempt / conditional on success at $p=0.999$.",
        "tab:fp1",
        ["method", "uncapped overflow semantics", "no observed failure (double / single)",
         "failure rate at $\\kappa{=}0.51$ / $0.55$",
         r"$\rho^{\rm attempt}/\rho^{\rm cond}$", "rel. time per return", "RNG stream"],
        body, colspec="lllllll")

    # Supplementary Table SFP1 - the complete failure grid.
    grid = [r for r in env_rows if r["scope"] == "pooled" and "failure_rate" in r]
    write_csv("table_sfp1_failure_grid.csv", grid)
    sfp1_body = []
    for r in sorted(grid, key=lambda r: (r["layer"], r["precision"], r["method"], r["kappa"]))[:400]:
        sfp1_body.append([
            latex_escape(r["layer"]), latex_escape(r["method"]), r["precision"],
            latex_escape(r["stdlib"]), f"{r['kappa']:g}", r["n_attempted"], r["n_failed"],
            fmt_rate(r.get("failure_rate"), r.get("failure_ci_lo"), r.get("failure_ci_hi"),
                     r.get("failure_is_upper_bound")),
            fmt_rate(r.get("avoidable_rate"), r.get("avoidable_ci_lo"),
                     r.get("avoidable_ci_hi"), r.get("avoidable_is_upper_bound")),
            fmt_rate(r.get("honest_floor_rate"), r.get("honest_floor_ci_lo"),
                     r.get("honest_floor_ci_hi"), r.get("honest_floor_is_upper_bound"))])
    write_tex_table(
        "table_sfp1_failure_grid.tex",
        "Complete failure grid: counts, rates and 95\\% Clopper--Pearson intervals for every "
        "layer, method, precision, standard library and $\\kappa$. Zero counts are reported "
        "as one-sided 95\\% upper bounds.",
        "tab:sfp1",
        ["layer", "method", "prec.", "stdlib", "$\\kappa$", "$N$", "failures",
         "failure rate", "avoidable", "honest"],
        sfp1_body, colspec="llllrrrlll")

    report["tables"] = {"table_fp1": "results/table_fp1.tex",
                        "table_sfp1": "results/table_sfp1_failure_grid.tex"}


def build_supplementary_tables(report: dict, perf_rows) -> None:
    val_path = out("loader_validation.csv")
    if os.path.exists(val_path):
        with open(val_path, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        stats_rows = [r for r in rows if r["family"] not in ("diagnostic",)]
        write_csv("table_sfp2_loader_validation.csv", stats_rows)
        body = [[latex_escape(r["case"]), latex_escape(r["method"]), r["precision"],
                 latex_escape(r["family"]), latex_escape(r["test"]),
                 f"{float(r['statistic']):.4g}" if r["statistic"] not in ("NA", "") else "--",
                 f"{float(r['pvalue']):.3g}" if r["pvalue"] not in ("NA", "") else "--",
                 r["n"], latex_escape(r["corrected_decision"])]
                for r in stats_rows[:400]]
        write_tex_table(
            "table_sfp2_loader_validation.tex",
            "Complete-loader validation: every predeclared statistic for cases C0--C6 and the "
            "two negative controls, with sample sizes and Holm-corrected decisions at a "
            "familywise $\\alpha$ of 0.01. NC1 and NC2 pass only when the injected defect is "
            "detected.",
            "tab:sfp2",
            ["case", "method", "prec.", "family", "test", "statistic", "$p$", "$n$",
             "decision"],
            body, colspec="lllllrrrl")

    summ = [r for r in perf_rows if r.get("scope") in ("summary", "ratio")]
    env_path = os.path.join(RAW, "environment.json")
    envj = json.load(open(env_path, encoding="utf-8")) if os.path.exists(env_path) else {}
    for r in summ:
        r["arch"] = envj.get("platform", {}).get("machine", "NA")
        r["os"] = envj.get("platform", {}).get("platform", "NA")
        r["cxxflags"] = envj.get("toolchains", {}).get("cxxflags", "NA")
        fp = (envj.get("floating_point") or {}).get(r.get("tag") or "", {}) or {}
        r["subnormals_flushed"] = fp.get("double_subnormals_flushed")
        r["rounding_mode"] = fp.get("rounding_mode")
    write_csv("table_sfp3_performance_environment.csv", summ)
    body = [[latex_escape(r["case"]), latex_escape(r["method"]), r["precision"],
             f"{r['kappa']:g}" if isinstance(r.get("kappa"), float) else r.get("kappa"),
             ("none" if r.get("cap") is None else f"{float(r['cap']):g}"),
             f"{r['median_seconds_per_attempt']:.3g}" if r.get("median_seconds_per_attempt") else "--",
             f"{r['median_seconds_per_return']:.3g}" if r.get("median_seconds_per_return") else "--",
             f"{r['attempts_per_return']:.3g}" if r.get("attempts_per_return") else
             (f"{r['total_attempts'] / max(r['total_returned'], 1):.3g}"
              if r.get("total_attempts") else "--"),
             f"{r['engine_calls_per_return']:.4g}" if r.get("engine_calls_per_return") else "--",
             latex_escape(r.get("arch", "NA")), latex_escape(str(r.get("rounding_mode")))]
            for r in summ if r.get("scope") == "summary"]
    write_tex_table(
        "table_sfp3_performance_environment.tex",
        "Time per attempt and per returned sample, random variates consumed, and the "
        "floating-point environment of each benchmark. Medians over ten timed blocks; the "
        "paired blockwise LOG/SPLIT ratios and their intervals are in the accompanying CSV.",
        "tab:sfp3",
        ["case", "method", "prec.", "$\\kappa$", "cap", "s/attempt", "s/return",
         "attempts/return", "engine calls/return", "arch", "rounding"],
        body, colspec="lllllrrrrll")


SOURCE_DATA_COLUMNS = {
    "failure_envelope.csv": [
        ("phase", "experiment phase that produced the row (E2 or E5)", "--"),
        ("scope", "`pooled` over seeds, or `seed` for one replicate", "--"),
        ("layer", "`paired` (common declared primitives) or `native` (each method's own primitive)", "--"),
        ("method", "QF, SPLIT or LOG", "--"),
        ("precision", "float or double; a real pipeline in that type, never drawn in double and cast", "--"),
        ("stdlib", "standard library whose Gamma generator was used", "--"),
        ("kappa", "spectral index", "dimensionless"),
        ("shape_a", "denominator Gamma shape kappa-1/2, formed in the working precision", "dimensionless"),
        ("n_attempted", "attempts; for uncapped runs this equals intended draws", "count"),
        ("n_finite", "attempts returning a finite three-vector", "count"),
        ("n_avoidable", "attempts lost to arithmetic whose intended vector was representable", "count"),
        ("n_honest_overflow", "attempts whose intended vector is not representable in the type", "count"),
        ("failure_rate / avoidable_rate / honest_rate", "per-attempt probabilities", "probability"),
        ("*_ci_lo / *_ci_hi", "two-sided 95% Clopper-Pearson interval", "probability"),
        ("*_is_upper_bound", "true when the count was zero and the interval is the one-sided 95% rule-of-three bound", "--"),
        ("accounting_ok", "true when finite + avoidable + honest equals attempts exactly", "--"),
        ("rotation_recoverable", "draws whose pre-rotation local vector overflows while every returned component is representable", "count"),
        ("audit_*_total / audit_*_audited", "size of each audit stratum and how much of it the 100-digit oracle adjudicated", "count"),
    ],
    "conditioning_bins.csv": [
        ("bin_kind", "`target_upper_tail_q` or `log_radius`", "--"),
        ("bin_hi / bin_lo", "bin edges; for q bins these are upper-tail probabilities", "probability / log units"),
        ("n_in_bin / n_success", "attempts in the bin and those returning a finite vector", "count"),
        ("success_rate", "Pr(success | bin)", "probability"),
        ("ci_lo / ci_hi", "two-sided 95% binomial interval", "probability"),
        ("merged", "true when the bin holds fewer than 20 attempts and is not interpreted", "--"),
    ],
    "tail_metrics.csv": [
        ("p", "target percentile", "probability"),
        ("target_log_r_quantile", "exact log R quantile of the uncapped target at this p", "log units"),
        ("rho_attempt", "Pr(R > r_p and success) / (1-p); failures stay outside the returned mass", "ratio"),
        ("rho_cond", "Pr(R > r_p | success) / (1-p); the law of finite survivors, which is also the law produced by silent redraw-to-success", "ratio"),
        ("rho_*_lo / rho_*_hi", "99% seed-stratified paired bootstrap percentile interval, 10000 resamples", "ratio"),
        ("conditional_log_r_quantile", "empirical quantile among survivors", "log units"),
        ("conditional_log_r_quantile_error", "survivor quantile minus target quantile", "log units"),
        ("quantile_resolved", "false when the sample cannot resolve the level; the estimate is then not interpreted", "--"),
        ("unconditional_failure_probability", "per-attempt failure probability, reported before any conditional statistic", "probability"),
    ],
    "loader_validation.csv": [
        ("family", "multiplicity family: radial, direction, independence, frame, cap_law, negative_control, or diagnostic", "--"),
        ("test", "statistic name", "--"),
        ("statistic / pvalue", "raw test statistic and p-value, preserved for source data", "--"),
        ("corrected_decision", "Holm-corrected decision at familywise alpha 0.01; `detected`/`MISSED` for negative controls", "--"),
        ("n_total / n_finite", "returned records and those with three finite components", "count"),
    ],
    "performance.csv": [
        ("scope", "`block` (one timed block), `summary` (median and IQR over blocks), `ratio` (paired blockwise LOG/SPLIT), `reproducibility`", "--"),
        ("seconds_per_attempt / seconds_per_return", "wall time per attempt and per returned sample", "s"),
        ("attempts_per_return", "cap rejection cost", "ratio"),
        ("engine_calls_per_return", "mt19937 outputs consumed per returned sample", "count"),
        ("median_time_per_return_ratio", "paired blockwise LOG/SPLIT ratio; 95% bootstrap interval in ratio_ci_lo/hi", "ratio"),
        ("repeat_identical", "true when the same build, method and seed reproduced identical counters", "--"),
    ],
    "portability.csv": [
        ("arch", "hardware architecture", "--"),
        ("available / completed", "whether the environment exists on this host and whether the job ran", "--"),
        ("reason_not_run", "exact reason and the command to run the job elsewhere; never blank", "--"),
    ],
    "scalar_validation.csv": [
        ("method", "LOG-ID or LOG-PUB", "--"),
        ("measured_acceptance", "accepted log-Gamma draws divided by proposals", "probability"),
        ("published_acceptance_eq2", "Liu, Martin and Syring (2017) Eq. (2) as printed", "probability"),
        ("published_acceptance_area_ratio", "Gamma(a+1)/(1+w), the ratio of areas implied by the same paper's target and envelope", "probability"),
        ("z_ks_stat / z_ks_p / z_cvm_stat / z_cvm_p", "KS and Cramer-von Mises against Exp(1) for Z = -log I_W(a,3/2)", "--"),
        ("w_ks_* / w_cvm_*", "the same against Beta(a,3/2) for W, computed only where W is numerically resolved", "--"),
        ("log_r_q*_estimate / _target / _error / _lo / _hi / _resolved / _covers_target",
         "log R quantiles with exact order-statistic intervals, Bonferroni-corrected over the five levels", "log units"),
    ],
}

MISSING_CODE = ("Missing values are the literal string `NA`. A blank cell never appears, so "
                "a blank can never be read as a zero.")


def write_source_data_readme(report: dict) -> None:
    lines = ["# Experiment 6 source data", "",
             f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.", "",
             "Every figure and table in this experiment is generated from the CSV files "
             "described here, and every plotted point is recoverable from a row in one of "
             "them. `make_figures.py` reads only these files; it never touches the bulk "
             "binaries under `raw/`.", "", MISSING_CODE, "",
             "Rates are per attempt unless a column says otherwise. Two-sided intervals are "
             "95% Clopper-Pearson; a zero count is reported as the one-sided 95% upper limit "
             "`1 - 0.05^(1/N)` and flagged with an `is_upper_bound` column, never as a zero "
             "probability. Tail ratios and method differences carry 99% seed-stratified "
             "bootstrap intervals from 10000 resamples. Quantile intervals are exact binomial "
             "order-statistic brackets, Bonferroni-corrected within each configuration.", ""]
    for name, cols in SOURCE_DATA_COLUMNS.items():
        if not os.path.exists(out(name)):
            continue
        lines += [f"## `{name}`", "", "| column | meaning | unit |", "|---|---|---|"]
        for col, meaning, unit in cols:
            lines.append(f"| `{col}` | {meaning} | {unit} |")
        lines.append("")
    lines += [
        "## Bulk files under `raw/`", "",
        "Not tracked in git, regenerable with the commands in `raw/environment.json`, and "
        "covered by `raw/raw_checksums.sha256`. Each begins with a 96-byte header naming its "
        "schema version, record size, configuration and record count; the reader refuses a "
        "file whose declared record count disagrees with its length.", "",
        "| file pattern | contents |", "|---|---|",
        "| `pilot/pilot_*.bin` | per-draw `log R` and `log W` from the E1 scalar pilot |",
        "| `conditioning/cond_*.bin` | per-attempt intended state and per-method terminal category |",
        "| `loader/loader_*.bin` | returned three-vectors, intended `log R`, status and attempt count |",
        "| `*/audit_*.bin` | primitives of audited attempts, handed to the 100-digit oracle |",
        "| `*/*.jsonl` | per-configuration counters, one JSON object per line |",
    ]
    with open(out("source_data_README.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  wrote results/source_data_README.md")


MANIFEST_FIELDS = [
    "phase", "algorithm", "variant", "precision", "kappa", "theta_perp", "theta_par",
    "ub_x", "ub_y", "ub_z", "cap", "seed", "n_attempted", "n_returned",
    "compiler", "compiler_version", "stdlib", "arch", "os", "cxxflags", "ftz_daz",
    "git_commit", "git_dirty", "sampler_header_sha256", "probe_sha256",
    "raw_file", "raw_sha256",
]


def write_manifest() -> None:
    """One row per produced configuration, naming the file it landed in and that file's hash.

    This is the index that makes the bulk data auditable without shipping it: the raw
    binaries are gitignored, so the manifest plus ``raw/raw_checksums.sha256`` is how a
    regenerated ``raw/`` is shown to be the one the committed results were computed from.
    """
    import hashlib

    env_path = os.path.join(RAW, "environment.json")
    envj = json.load(open(env_path, encoding="utf-8")) if os.path.exists(env_path) else {}
    deps = envj.get("dependencies", {})
    git = envj.get("git", {})
    fp = envj.get("floating_point") or {}

    cache: dict[str, str] = {}

    def digest(path: str) -> str:
        if path in cache:
            return cache[path]
        h = hashlib.sha256()
        try:
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(chunk)
            cache[path] = h.hexdigest()
        except OSError:
            cache[path] = "NA"
        return cache[path]

    phases = {"baseline": "E0", "pilot": "E1", "mechanism": "E2", "conditioning": "E3",
              "loader": "E4", "portability": "E5", "performance": "E6"}
    rows = []
    for phase_dir, phase_id in phases.items():
        for r in IO.load_phase(RAW, phase_dir, SMOKE):
            tag = r.get("tag", "")
            probe = os.path.join(HERE, f"exp6_probe_{tag}.exe")
            bulk = r.get("raw_file")
            src = bulk if bulk and os.path.exists(bulk) else None
            if src is None:
                sub = "smoke" if SMOKE else phase_dir
                cand = os.path.join(RAW, sub, f"{phase_dir}_{tag}.jsonl")
                src = cand if os.path.exists(cand) else None
            fpe = fp.get(tag) or {}
            rows.append({
                "phase": phase_id,
                "algorithm": r.get("method", "NA"),
                "variant": r.get("layer", "NA") + (
                    f"/{r['case']}" if r.get("case") else ""),
                "precision": r.get("precision", "NA"),
                "kappa": r.get("kappa"),
                "theta_perp": 1.0,
                "theta_par": r.get("theta_ratio", 1.0),
                "ub_x": "NA", "ub_y": "NA", "ub_z": "NA",
                "cap": r.get("cap"),
                "seed": r.get("seed"),
                "n_attempted": r.get("n_attempted"),
                "n_returned": r.get("n_returned", r.get("n_attempted")),
                "compiler": r.get("compiler", "NA").split()[0] if r.get("compiler") else "NA",
                "compiler_version": r.get("compiler", "NA"),
                "stdlib": r.get("stdlib", "NA"),
                "arch": r.get("arch", "NA"),
                "os": envj.get("platform", {}).get("platform", "NA"),
                "cxxflags": envj.get("toolchains", {}).get("cxxflags", "NA"),
                "ftz_daz": ("flushed" if fpe.get("double_subnormals_flushed")
                            else "ieee_subnormals" if fpe else "NA"),
                "git_commit": git.get("commit", "NA"),
                "git_dirty": git.get("repo_dirty"),
                "sampler_header_sha256":
                    (deps.get("../../cpp/bi_kappa_distribution.H") or {}).get("sha256", "NA"),
                "probe_sha256": digest(probe) if os.path.exists(probe) else "NA",
                "raw_file": os.path.relpath(src, HERE) if src else "NA",
                "raw_sha256": digest(src) if src else "NA",
            })
    path = os.path.join(RAW, "manifest.csv")
    os.makedirs(RAW, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in MANIFEST_FIELDS})
    print(f"  wrote raw/manifest.csv  ({len(rows)} rows, {len(cache)} files hashed)")


# ===========================================================================
# Gates and final report
# ===========================================================================
def oracle_summary() -> dict:
    path = out("oracle_audit.jsonl")
    if not os.path.exists(path):
        return {"ran": False}
    rows = IO.read_jsonl(path)
    if not rows:
        return {"ran": False}
    return {
        "ran": True,
        "files": len({r["file"] for r in rows}),
        "records": sum(r["n_records"] for r in rows),
        "disagreements": sum(r["disagreements"] for r in rows),
        "conversion_failures": sum(r["conversion_failures"] for r in rows),
        "max_rel_err_split_subnormal": max(
            (r.get("max_rel_err_split_subnormal_denominator", 0.0) for r in rows), default=0.0),
        "max_rel_err_log_subnormal": max(
            (r.get("max_rel_err_log_subnormal_denominator", 0.0) for r in rows), default=0.0),
        "n_subnormal_denominator_audited": sum(
            r.get("n_subnormal_denominator", 0) for r in rows),
        "n_radius_compared": sum(r.get("n_radius_compared", 0) for r in rows),
        "by_precision": {
            p: {
                "max_rel_err_split": max((r["max_rel_err_split"] for r in rows
                                          if r.get("precision") == p), default=float("nan")),
                "max_rel_err_log": max((r["max_rel_err_log"] for r in rows
                                        if r.get("precision") == p), default=float("nan")),
                "max_abs_log_x2_error": max((r["max_abs_log_x2_error"] for r in rows
                                             if r.get("precision") == p), default=float("nan")),
            } for p in ("double", "float")
        },
    }


def evaluate_gates(report: dict) -> dict:
    env_path = os.path.join(RAW, "environment.json")
    envj = json.load(open(env_path, encoding="utf-8")) if os.path.exists(env_path) else {}
    oracle = oracle_summary()
    report["oracle"] = oracle

    e1 = report.get("E1", {})
    e2 = report.get("E2", {})
    e3 = report.get("E3", {})
    e4 = report.get("E4", {})
    e5 = report.get("E5", {})
    e6 = report.get("E6", {})

    gates = {}
    gates["G0"] = {
        "name": "primary-source audit",
        "pass": bool(e1.get("evaluator_validation_passed") and e1.get("selected")),
        "evidence": "config/pilot.json records the exact equations, stated parameter domain, "
                    "output scale and limitations of both candidates with their DOIs; the "
                    "frozen log-q evaluator agrees with an arbitrary-precision incomplete "
                    "beta to better than 1e-10; the LOG-PUB transcription reproduces the "
                    "acceptance rate implied by the paper's own target and envelope and the "
                    "sampled law satisfies E[log X] = psi(alpha).",
    }
    sel = e1.get("summary", {}).get(e1.get("selected"), {})
    gates["G1"] = {
        "name": "scalar correctness",
        "pass": bool(sel.get("passes")),
        "evidence": f"selected primitive {e1.get('selected')}: "
                    f"{sel.get('log_primitive_failures', 'NA')} log-domain failures, "
                    f"{sel.get('configurations_with_a_rejection', 'NA')} of "
                    f"{sel.get('configurations', 'NA')} pilot configurations with a "
                    "Holm-corrected rejection.",
    }
    gates["G2"] = {
        "name": "mechanism closure",
        "pass": bool(e2.get("accounting_failures") == 0 and
                     e2.get("honest_floor_consistent_across_methods") and
                     e2.get("configurations_outside_boundary_band") == 0 and
                     oracle.get("ran") and oracle.get("disagreements") == 0 and
                     oracle.get("conversion_failures") == 0),
        "evidence": f"{e2.get('accounting_failures')} accounting-identity failures across "
                    f"{e2.get('configurations')} configurations; every method's honest count "
                    f"agrees with the reference floor to within the audited near-threshold "
                    f"band ({e2.get('configurations_outside_boundary_band')} exceptions); "
                    f"the log path loses "
                    f"{e2.get('log_avoidable_loss_double_total')} representable draws in "
                    f"double across the whole ladder; the independent cpp_dec_float_100 "
                    f"oracle adjudicated {oracle.get('records', 0)} audited attempts with "
                    f"{oracle.get('disagreements', 'NA')} disagreements.",
    }
    gates["G3"] = {
        "name": "complete-loader fidelity",
        "pass": bool(e4.get("candidate_cells_with_failure") == 0 and
                     e4.get("negative_controls_detected") and e3.get("tail_rows")),
        "evidence": f"{e4.get('cells')} validation cells; "
                    f"{e4.get('candidate_cells_with_failure')} of the candidate's cells "
                    f"contain a failure; negative controls "
                    f"{'detected' if e4.get('negative_controls_detected') else 'NOT detected'}; "
                    f"{e3.get('tail_rows')} tail-metric rows quantify the conditioning. "
                    f"{len(e4.get('released_path_conditional_rejections') or [])} released-path "
                    "cells reject the radial law on their finite survivors, which is the "
                    "conditioning result rather than a defect.",
    }
    missing = e5.get("environments_missing") or []
    gates["G4"] = {
        "name": "portability",
        "pass": False if missing else bool(e5.get("environments_completed")),
        "evidence": f"completed: {e5.get('environments_completed')}; "
                    f"not available on this host: {missing}. An unavailable architecture is "
                    "recorded with the exact command to run it elsewhere; emulation is not "
                    "accepted as a substitute, so this gate stays open.",
    }
    ratio = e6.get("worst_time_ratio")
    gates["G5"] = {
        "name": "operational viability",
        "pass": bool(e6.get("same_seed_reproducible") and
                     isinstance(ratio, float) and math.isfinite(ratio)),
        "evidence": f"worst paired LOG/SPLIT time per returned sample: "
                    f"{ratio:.2f}x; same-seed reproducibility "
                    f"{e6.get('same_seed_reproducible')}. Whether that cost is acceptable "
                    "for a header-only C++11 reference implementation is an author decision, "
                    "not a statistical one." if isinstance(ratio, float) else "not measured",
    }
    clean = bool(envj.get("git", {}).get("dependencies_clean"))
    gates["G6"] = {
        "name": "reproducible artifact",
        "pass": bool(clean and os.path.exists(os.path.join(HERE, "checksums.sha256"))),
        "evidence": f"dependency set clean: {clean}; commit "
                    f"{envj.get('git', {}).get('commit', 'NA')[:12]}; "
                    f"exploratory flag: {envj.get('exploratory')}.",
    }

    mandatory_for_claim = ["G0", "G1", "G2", "G3", "G6"]
    if all(gates[g]["pass"] for g in gates):
        verdict = "GO"
    elif all(gates[g]["pass"] for g in mandatory_for_claim):
        verdict = "PARTIAL"
    else:
        verdict = "NO-GO"
    report["gates"] = gates
    report["verdict"] = verdict
    return gates


def write_analysis_report(report: dict) -> None:
    gates = report["gates"]
    verdict = report["verdict"]
    e1, e2, e3, e4 = (report.get(k, {}) for k in ("E1", "E2", "E3", "E4"))
    oracle = report.get("oracle", {})

    lines = [f"{verdict}", "",
             "# Experiment 6 - analysis report", "",
             f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}.", "",
             "## Gates", "",
             "| gate | name | verdict | evidence |", "|---|---|---|---|"]
    for g in sorted(gates):
        lines.append(f"| {g} | {gates[g]['name']} | "
                     f"{'PASS' if gates[g]['pass'] else 'OPEN'} | {gates[g]['evidence']} |")

    lines += ["", "## What the verdict means", ""]
    if verdict == "GO":
        lines.append("Every gate passed. The attributed log-domain integration may be reported "
                     "as a mitigation and considered for production adoption.")
    elif verdict == "PARTIAL":
        lines += [
            "The gates that a manuscript mitigation claim requires (G0-G3, G6) passed; at "
            "least one gate that production adoption requires did not.",
            "",
            "The log path stays an unreleased prototype. The measured operating envelope, the "
            "mechanism decomposition and the conditioning result stand on their own and are "
            "reportable within the boundary those gates establish. The released default in "
            "`cpp/bi_kappa_distribution.H` is not changed by this experiment, and no public "
            "API is altered.",
        ]
    else:
        lines += ["A gate that a manuscript mitigation claim requires did not pass. Retain the "
                  "current operating-envelope result as the finite-precision layer, keep the "
                  "log path explicitly unreleased, and do not substitute a different "
                  "high-level sampler without a new scope decision."]

    lines += ["", "## Headline numbers", "",
              f"- Selected log primitive: **{e1.get('selected')}** - {e1.get('reason')}",
              f"- Accounting-identity failures: {e2.get('accounting_failures')} across "
              f"{e2.get('configurations')} paired configurations.",
              f"- Independent 100-digit oracle: {oracle.get('records', 0)} audited attempts, "
              f"{oracle.get('disagreements', 'NA')} disagreements, "
              f"{oracle.get('conversion_failures', 'NA')} conversion failures.",
              f"- Complete-loader cells: {e4.get('cells')}, of which "
              f"{e4.get('candidate_cells_with_failure')} of the candidate's contain a "
              f"failure; negative controls "
              f"{'detected' if e4.get('negative_controls_detected') else 'NOT detected'}.",
              f"- Released-path cells whose radial law is rejected on the finite survivors: "
              f"{len(e4.get('released_path_conditional_rejections') or [])}. That is the "
              "conditioning result: discarding the failures changes the law of what is "
              "kept, and the log path passes the same test on the same configuration."]
    if oracle.get("ran"):
        lines += [
            "",
            "## Accuracy of the surviving draws",
            "",
            "The oracle also re-derives each formation's radius exactly, which is how this "
            "experiment can say something the counters cannot. Surviving a draw is not the "
            "same as getting it right. On the "
            f"{oracle.get('n_subnormal_denominator_audited', 0):,} audited draws whose "
            "denominator Gamma landed in the subnormal range - the draws the split form is "
            "closest to losing altogether - the split form's returned radius carries up to "
            f"{oracle.get('max_rel_err_split_subnormal', float('nan')):.3g} relative error, "
            f"against {oracle.get('max_rel_err_log_subnormal', float('nan')):.3g} for the log "
            "path. A subnormal denominator has already lost most of its significand, and the "
            "square root halves the exponent without restoring it; the log path never "
            "materializes the denominator at all. This is conditional on the denominator "
            "being subnormal, a stratum the audit covers at a fixed sampling rate, and it "
            "says nothing about draws where the denominator is normal - there the two agree "
            "to rounding.",
        ]
    lines += [
        "",
        "## Claim boundary",
        "",
        "Set by the independent novelty audit, not by this result. The Gamma-ratio Kappa "
        "construction, log-domain Gamma generation, the small-shape Gamma underflow and the "
        "low-parameter denominator hazard are all established prior art and are attributed "
        "as such.",
        "",
        "Supportable on the evidence here:",
        "",
        "- Forming the intermediate quotient `X1/X2` can overflow even where its square root "
        "is representable; evaluating the algebraically identical `sqrt(X1)/sqrt(X2)` removes "
        "that mode, and in paired tests the two formulations receive identical variates.",
        "- That change does not remedy zeros produced by the small-shape Gamma primitive and "
        "does not extend the mathematical support of the target law.",
        "- Working in the log domain is established for small-shape Gamma variates; the "
        "contribution here is its integration and validation in this loader.",
        "- Failure rates are reported as observed counts with intervals, and a configuration "
        "with no observed failure is reported as `no failures observed in N draws under the "
        "tested configuration`, with its one-sided upper bound.",
        "",
        "Not supportable, and not claimed:",
        "",
        "- discovery of the small-shape Gamma underflow problem, or of division by zero in a "
        "low-parameter Kappa loader;",
        "- a first or novel log-domain Gamma generator, or a novel Gamma-ratio, Beta-prime, "
        "Student-t or rejection Kappa sampler;",
        "- that the split expression solves the low-kappa finite-precision problem;",
        "- that the log representation is exact - it is a stable representation carrying "
        "about one ulp of `(log U)/a`;",
        "- that the log path removes all avoidable loss, beyond what the final-vector oracle "
        "has closed;",
        "- any universal mathematical lower bound on kappa, or `reliable down to`, without a "
        "sample size, platform, failure criterion and confidence bound;",
        "- that both standard libraries use an identical `Y U^(1/a)` implementation;",
        "- complete-loader correctness inferred from tests run after deleting non-finite "
        "outputs. Conditional results are labelled conditional throughout.",
        "",
        "## Unresolved exceptions", ""]
    unresolved = [f"- {g} ({gates[g]['name']}): {gates[g]['evidence']}"
                  for g in sorted(gates) if not gates[g]["pass"]]
    lines += unresolved or ["- none"]
    lines += ["", "## Decisive source data", "",
              "- `baseline_reproduction.md`, `baseline_envelope.csv` (E0)",
              "- `pilot_decision.md`, `scalar_validation.csv`, "
              "`log_q_evaluator_validation.csv` (E1)",
              "- `failure_envelope.csv` (E2), `oracle_audit.jsonl`, "
              "`oracle_disagreements.jsonl`",
              "- `conditioning_bins.csv`, `tail_metrics.csv` (E3)",
              "- `loader_validation.csv`, `validation_matrix.md` (E4)",
              "- `portability.csv` (E5), `performance.csv` (E6)",
              "- `source_data_README.md` defines every column."]
    with open(out("analysis_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote results/analysis_report.md  ({verdict})")


def main() -> None:
    global SMOKE
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="analyze raw/smoke/ instead of the production phases; smoke output "
                         "is diagnostic only and must never feed a figure or table")
    args = ap.parse_args()
    SMOKE = args.smoke
    if SMOKE:
        print("  SMOKE MODE - diagnostic only")

    os.makedirs(RESULTS, exist_ok=True)
    report = {"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "smoke": SMOKE}

    env_path = os.path.join(RAW, "environment.json")
    if os.path.exists(env_path):
        report["environment"] = json.load(open(env_path, encoding="utf-8"))

    phase_baseline(report)
    phase_pilot(report)
    env_rows = phase_mechanism(report)
    phase_conditioning(report)
    phase_loader(report)
    port_rows = phase_portability(report)
    perf_rows = phase_performance(report)

    all_env = env_rows + [r for r in port_rows if r.get("scope")]
    if all_env:
        build_tables(report, all_env, perf_rows, out("tail_metrics.csv"))
    build_supplementary_tables(report, perf_rows)
    write_source_data_readme(report)
    write_manifest()
    evaluate_gates(report)
    write_analysis_report(report)

    with open(out("exp6_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
        fh.write("\n")
    print("  wrote results/exp6_results.json")
    print(f"\n  VERDICT: {report['verdict']}")


if __name__ == "__main__":
    main()
