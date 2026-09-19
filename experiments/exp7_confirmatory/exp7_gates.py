#!/usr/bin/env python3
"""Gate evaluation and the final verdict, per ``PROTOCOL.md`` §6.

Every gate below is a computation over evidence the run produced. None is satisfied by a
prose string, by a file existing, or by a number being finite -- which is how Experiment 6's
G0 passed on an unrelated computation, its G3 on ``all([])`` over an empty list, its G4 on a
hard-coded constant, its G5 on ``isfinite(ratio)``, and its G6 on ``os.path.exists``.

Two rules are enforced structurally rather than by care:

* **Absent evidence fails.** Every gate states the evidence it requires; if a required input
  is missing the gate FAILS. A gate can never pass because nothing was measured.
* **The protocol is the only source of thresholds.** ``exp7_families.Protocol`` raises on a
  missing key, so a criterion that was not pre-registered cannot be applied.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import exp7_families as F

GATE_NAMES = {
    "G0": "primary-source audit",
    "G1": "scalar correctness",
    "G2": "mechanism closure",
    "G3": "complete-loader fidelity",
    "G4": "portability",
    "G5": "operational viability",
    "G6": "reproducible artifact",
}

# G0-G3 and G6 are what a manuscript mitigation claim requires; all seven are what production
# adoption requires.  Frozen in the plan, restated in PROTOCOL.md §6.
MANUSCRIPT_GATES = ("G0", "G1", "G2", "G3", "G6")


class Gate:
    def __init__(self, name: str, status: str, evidence: str, inputs: dict | None = None):
        assert status in ("PASS", "FAIL", "OPEN"), status
        self.name = name
        self.status = status
        self.evidence = evidence
        self.inputs = inputs or {}

    @property
    def passed(self) -> bool:
        return self.status == "PASS"

    def as_dict(self) -> dict:
        return {"gate": self.name, "title": GATE_NAMES[self.name], "status": self.status,
                "evidence": self.evidence, "inputs": self.inputs}


def _missing(name: str, what: str) -> Gate:
    return Gate(name, "FAIL",
                f"required evidence is absent: {what}. An absent measurement fails this gate; "
                "it never passes it.")


# ---------------------------------------------------------------------------
def gate_G0(proto: F.Protocol, ev: dict) -> Gate:
    """Evaluator accuracy, a machine-checked attribution checklist, and a predictor named
    in advance."""
    tol = float(proto.get("statistics", "log_q_oracle_tolerance"))
    worst = ev.get("log_q_worst_abs_error")
    if worst is None:
        return _missing("G0", "the log-q evaluator was never validated against the "
                              "arbitrary-precision oracle")
    checklist = proto.get("attribution_checklist")
    required_fields = ("method", "doi", "source", "equation", "domain", "output_scale",
                       "limitation")
    incomplete = [c.get("method", "?") for c in checklist
                  if any(not c.get(f) for f in required_fields)]
    predictor = proto.get("acceptance_rate_predictor", "rule")
    ok = bool(worst < tol) and not incomplete and bool(predictor)
    return Gate("G0", "PASS" if ok else "FAIL",
                f"log-q evaluator agrees with the 60-digit oracle to {worst:.3g} "
                f"(tolerance {tol:g}); attribution checklist complete for "
                f"{len(checklist) - len(incomplete)}/{len(checklist)} cited methods"
                + (f"; INCOMPLETE: {incomplete}" if incomplete else "")
                + f"; acceptance-rate predictor pre-registered as '{predictor}'",
                {"log_q_worst_abs_error": worst, "tolerance": tol,
                 "incomplete_attributions": incomplete})


def gate_G1(proto: F.Protocol, fams: dict, ev: dict) -> Gate:
    """F1-F4 pass on in-domain configurations, the benign control exercises the candidate,
    and avoidable loss is exactly zero with no elastic band."""
    needed = ("F1_radial_law", "F2_quantile_coverage", "F3_quantile_direction",
              "F4_upper_tail_mass")
    absent = [f for f in needed if f not in fams]
    if absent:
        return _missing("G1", f"families {absent} produced no result")
    failed = [f for f in needed if not fams[f].passed]
    avoidable = ev.get("candidate_avoidable_total")
    if avoidable is None:
        return _missing("G1", "the candidate's avoidable-loss count was never computed")
    benign = ev.get("benign_control_exercises_candidate")
    if benign is None:
        return _missing("G1", "no record of whether the benign control exercised the "
                              "candidate rather than a fallback")
    ok = not failed and int(avoidable) == 0 and bool(benign)
    return Gate("G1", "PASS" if ok else "FAIL",
                f"{len(needed) - len(failed)}/{len(needed)} scalar families pass"
                + (f" (failed: {failed})" if failed else "")
                + f"; candidate avoidable loss {avoidable} (required exactly 0)"
                + f"; benign control exercises the candidate: {bool(benign)}",
                {"failed_families": failed, "avoidable": avoidable,
                 "benign_control_ok": bool(benign)})


def gate_G2(proto: F.Protocol, ev: dict) -> Gate:
    """Accounting exact, oracle clean, honest count equal to the floor, and every audit
    stratum covered at its declared rate."""
    for key, what in (("accounting_failures", "the accounting identity was never checked"),
                      ("oracle_disagreements", "the oracle never adjudicated anything"),
                      ("oracle_conversion_failures", "no oracle conversion record"),
                      ("honest_minus_floor_unadjudicated",
                       "no comparison of the honest count against the oracle floor")):
        if ev.get(key) is None:
            return _missing("G2", what)
    declared = proto.get("audit", "strata_rates")
    coverage = ev.get("audit_coverage") or {}
    short = {k: coverage.get(k) for k, r in declared.items()
             if coverage.get(k) is not None and coverage[k] + 1e-12 < r}
    never = [k for k in declared if k not in coverage]
    ok = (int(ev["accounting_failures"]) == 0 and int(ev["oracle_disagreements"]) == 0
          and int(ev["oracle_conversion_failures"]) == 0
          and int(ev["honest_minus_floor_unadjudicated"]) == 0
          and not short and not never)
    return Gate("G2", "PASS" if ok else "FAIL",
                f"{ev['accounting_failures']} accounting-identity failures; "
                f"{ev['oracle_disagreements']} oracle disagreements in "
                f"{ev.get('oracle_audited_records', '?')} adjudicated attempts; "
                f"{ev['honest_minus_floor_unadjudicated']} honest-count discrepancies left "
                "unadjudicated"
                + (f"; strata below their declared rate: {short}" if short else "")
                + (f"; strata never reported: {never}" if never else ""),
                {"audit_coverage": coverage, "declared_rates": declared})


def gate_G3(proto: F.Protocol, fams: dict, ev: dict) -> Gate:
    """F5 jointly, F6 by measured power with every control present, and conditional cells
    labelled."""
    if "F5_loader_battery" not in fams:
        return _missing("G3", "the loader battery produced no result")
    if "F6_negative_controls" not in fams:
        return _missing("G3", "no negative controls ran")
    f5, f6 = fams["F5_loader_battery"], fams["F6_negative_controls"]
    unlabelled = ev.get("conditional_cells_without_loss_fraction")
    if unlabelled is None:
        return _missing("G3", "no record of which loader cells are conditional on success")
    ok = f5.passed and f6.passed and int(unlabelled) == 0
    return Gate("G3", "PASS" if ok else "FAIL",
                f"F5 {'passes' if f5.passed else 'FAILS'} jointly over its cells "
                f"({f5.detail}); F6 {'passes' if f6.passed else 'FAILS'} ({f6.detail}); "
                f"{unlabelled} conditional cells lack a loss fraction",
                {"F5": f5.as_dict(), "F6": f6.as_dict()})


def gate_G4(proto: F.Protocol, fams: dict, ev: dict) -> Gate:
    """Bitwise across standard libraries, equivalence across architectures, and an
    unavailable environment leaves the gate OPEN rather than passing."""
    env = ev.get("environments") or []
    if not env:
        return _missing("G4", "no environment records were produced")
    native = [e for e in env if e.get("execution") == "native" and e.get("completed")]
    translated = [e for e in env if e.get("execution") == "translated"]
    pending = [e for e in env if not e.get("completed")]
    if "F7_portability" not in fams:
        return _missing("G4", "no portability comparison was produced")
    f7 = fams["F7_portability"]
    archs = {e.get("arch") for e in native}
    detail = (f"{len(native)} native environments completed "
              f"({sorted({(e.get('arch'), e.get('stdlib')) for e in native})}); "
              f"{len(translated)} translated results recorded as corroborating only and "
              f"excluded from the decision; {f7.detail}")
    if not f7.passed:
        return Gate("G4", "FAIL", detail + f"; offenders: {f7.offenders}",
                    {"F7": f7.as_dict()})
    if pending or len(archs) < 2:
        return Gate("G4", "OPEN",
                    detail + "; "
                    + (f"{len(pending)} environments not yet run: "
                       f"{[(e.get('arch'), e.get('stdlib')) for e in pending]}. "
                       if pending else "")
                    + ("only one architecture has native results, so the "
                       "cross-architecture claim is untested. " if len(archs) < 2 else "")
                    + "The exact commands are in results/portability_remote.md. Emulation "
                      "is not accepted as a substitute, so the gate stays open rather than "
                      "being quietly passed.",
                    {"F7": f7.as_dict(), "pending": pending})
    return Gate("G4", "PASS", detail, {"F7": f7.as_dict()})


def gate_G5(proto: F.Protocol, ev: dict) -> Gate:
    """A computed rule: the paired ratio's interval against a pre-registered bound."""
    bound = float(proto.get("gates", "G5", "time_per_return_ratio_bound"))
    ratio = ev.get("time_per_return_ratio")
    hi = ev.get("time_per_return_ratio_hi")
    repro = ev.get("same_seed_reproducible")
    documented = ev.get("rng_stream_break_documented")
    if ratio is None or hi is None:
        return _missing("G5", "no paired performance ratio with an interval")
    if repro is None:
        return _missing("G5", "same-seed reproducibility was never checked")
    if documented is None:
        return _missing("G5", "no record that the RNG-stream break is documented")
    ok = float(hi) < bound and bool(repro) and bool(documented)
    return Gate("G5", "PASS" if ok else "FAIL",
                f"candidate/LEGACY time per returned sample {ratio:.3g} "
                f"(cluster-bootstrap upper limit {hi:.3g}; pre-registered bound {bound:g}); "
                f"same-seed reproducible: {bool(repro)}; RNG-stream break documented: "
                f"{bool(documented)}",
                {"ratio": ratio, "hi": hi, "bound": bound})


def gate_G6(proto: F.Protocol, ev: dict) -> Gate:
    """Verification that actually verifies, byte-identical regeneration, and an exact-version
    archive identifier."""
    for key, what in (("make_verify_exit_code", "`make verify` was never run"),
                      ("make_reverify_identical", "derived artifacts were never regenerated "
                                                  "and compared"),
                      ("dependencies_clean", "no provenance record for the production run"),
                      ("baseline_comparison_resolved",
                       "the baseline comparison against the previous release was never "
                       "resolved")):
        if ev.get(key) is None:
            return _missing("G6", what)
    archive = ev.get("archive_identifier")
    protocol_sha = ev.get("protocol_sha256")
    sha_ok = protocol_sha == proto.sha256
    ok = (int(ev["make_verify_exit_code"]) == 0 and bool(ev["make_reverify_identical"])
          and bool(ev["dependencies_clean"]) and bool(ev["baseline_comparison_resolved"])
          and sha_ok and bool(archive))
    status = "PASS" if ok else ("OPEN" if (not archive and int(ev["make_verify_exit_code"]) == 0
                                           and ev["make_reverify_identical"]) else "FAIL")
    return Gate("G6", status,
                f"make verify exit {ev['make_verify_exit_code']}; derived artifacts "
                f"regenerate byte-identically: {bool(ev['make_reverify_identical'])}; "
                f"dependency set clean: {bool(ev['dependencies_clean'])}; baseline "
                f"comparison resolved: {bool(ev['baseline_comparison_resolved'])}; "
                f"protocol hash matches the frozen document: {sha_ok}; "
                + (f"archive identifier {archive}" if archive
                   else "NO exact-version archive identifier yet"),
                {"archive_identifier": archive, "protocol_sha256_matches": sha_ok})


# ---------------------------------------------------------------------------
def evaluate(proto: F.Protocol, families: dict, evidence: dict) -> dict:
    gates = {
        "G0": gate_G0(proto, evidence),
        "G1": gate_G1(proto, families, evidence),
        "G2": gate_G2(proto, evidence),
        "G3": gate_G3(proto, families, evidence),
        "G4": gate_G4(proto, families, evidence),
        "G5": gate_G5(proto, evidence),
        "G6": gate_G6(proto, evidence),
    }
    if all(g.passed for g in gates.values()):
        verdict = "GO"
    elif all(gates[k].passed for k in MANUSCRIPT_GATES):
        verdict = "PARTIAL"
    else:
        verdict = "NO-GO"
    return {"verdict": verdict, "gates": gates}


def write_report(path: str, proto: F.Protocol, result: dict, families: dict,
                 evidence: dict, generated_utc: str) -> None:
    """``results/analysis_report.md``: the verdict on the first line, then the gate table."""
    gates = result["gates"]
    unresolved = [g for g in gates.values() if not g.passed]
    L = [result["verdict"], "", "# Experiment 7 - confirmatory analysis report", "",
         f"Generated {generated_utc}.", "",
         f"Protocol `config/protocol.json` version "
         f"{proto.get('protocol_version')}, SHA-256 `{proto.sha256}`.", ""]

    L += ["## Gates", "", "| gate | name | verdict | evidence |", "|---|---|---|---|"]
    for k in sorted(gates):
        g = gates[k]
        L.append(f"| {k} | {GATE_NAMES[k]} | {g.status} | {g.evidence} |")
    L += ["", "## Families", "",
          "| family | alpha | passed | statistic | p | detail |", "|---|---:|---|---:|---:|---|"]
    for name in sorted(families):
        f = families[name]
        a = "-" if f.alpha != f.alpha else f"{f.alpha:g}"
        st = "-" if f.statistic != f.statistic else f"{f.statistic:.4g}"
        pv = "-" if f.p_value != f.p_value else f"{f.p_value:.4g}"
        L.append(f"| {name} | {a} | {'yes' if f.passed else 'NO'} | {st} | {pv} | {f.detail} |")

    L += ["", "## What the verdict means", ""]
    if result["verdict"] == "GO":
        L += ["Every gate passed on the pre-registered rules, on a seed block disjoint from "
              "every other experiment in this repository, against a protocol committed "
              "before the first draw existed."]
    elif result["verdict"] == "PARTIAL":
        L += ["The gates a manuscript mitigation claim requires (G0-G3, G6) passed; at least "
              "one gate that production adoption requires did not. The precise unresolved "
              "item is listed below and is not argued away."]
    else:
        L += ["At least one gate a manuscript claim requires did not pass. The failure is "
              "preserved. Per PROTOCOL.md §8 the only legitimate next step is to identify a "
              "concrete implementation defect, fix it, freeze a new implementation hash and "
              "a new protocol, draw a further disjoint seed block, and rerun. Repeating the "
              "run to obtain a more favourable result is not permitted."]
    L += ["", f"Expected false-failure rate of this design, computed from the frozen family "
              f"structure before the run: at most "
              f"{proto.get('statistics', 'expected_false_failure_bound'):g}.", ""]

    L += ["## Unresolved", ""]
    if unresolved:
        for g in unresolved:
            L.append(f"- **{g.name} ({GATE_NAMES[g.name]}) - {g.status}.** {g.evidence}")
    else:
        L.append("None.")
    L += ["", "## Decisive source data", "",
          "- `PROTOCOL.md`, `config/protocol.json`, `config/power_study.json`, "
          "`config/honest_floor.md` - the frozen rules and what they can detect",
          "- `results/scalar_validation.csv`, `results/tail_metrics.csv` (F1-F4)",
          "- `results/failure_envelope.csv`, `results/oracle_audit.jsonl` (G2)",
          "- `results/loader_validation.csv`, `results/negative_controls.csv` (F5, F6)",
          "- `results/portability.csv`, `results/portability_remote.md` (F7, G4)",
          "- `results/performance.csv` (G5)",
          "- `results/source_data_README.md` defines every column.", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
