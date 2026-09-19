"""Emit ``src/exp7_protocol.H`` from the frozen ``config/protocol.json``.

Run with:  uv run --project ../../python python src/gen_protocol_header.py

The probe has to record order statistics at exactly the index brackets the protocol
froze, and has to know the quantile levels, the upper-tail thresholds and the incomplete-
beta switch point that the frozen statistics use.  Writing those numbers into the C++ by
hand would let the source and ``config/protocol.json`` drift apart silently, and the
protocol is explicit that every decision rule is read from the JSON rather than from a
constant in the code.  This script is the one place the two are joined: the header is
generated, never edited, and it carries the SHA-256 of the JSON it came from so that every
raw file the probe writes can name the protocol revision it was produced under.

``src/exp7_common.H`` declares the kappa ladder and both seed blocks independently, in
ordinary C++.  The probe compares its own declarations against this generated header at
start-up and aborts on any difference, so the two statements of the matrix check each
other instead of one being trusted.

Deterministic: same input, byte-identical output.  ``make protocol-check`` regenerates
``config/protocol.json`` and diffs it, and the header is rebuilt whenever the JSON changes.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROTOCOL = os.path.join(ROOT, "config", "protocol.json")
OUT = os.path.join(HERE, "exp7_protocol.H")


def c_double(x: float) -> str:
    """A double literal that round-trips exactly."""
    return repr(float(x))


def main() -> int:
    with open(PROTOCOL, "rb") as fh:
        raw = fh.read()
    sha = hashlib.sha256(raw).hexdigest()
    p = json.loads(raw.decode("utf-8"))

    matrix = p["matrix"]
    stats = p["statistics"]
    levels = list(matrix["quantile_levels"])
    n_scalar = int(matrix["n_scalar"])

    # The F2 brackets depend only on (n, p); the protocol stores one cell per
    # (precision, kappa, seed, level), so collapse them and refuse any inconsistency.
    brackets: dict[float, tuple[int, int, float]] = {}
    for cell in p["F2_cells"]:
        if not cell.get("resolved"):
            continue
        key = float(cell["p"])
        val = (int(cell["lo_index"]), int(cell["hi_index"]), float(cell["achieved_coverage"]))
        if key in brackets and brackets[key] != val:
            sys.stderr.write(
                f"exp7: protocol.json gives two different brackets for p={key}\n")
            return 1
        brackets[key] = val
    missing = [lv for lv in levels if lv not in brackets]
    if missing:
        sys.stderr.write(f"exp7: protocol.json resolves no bracket for levels {missing}\n")
        return 1

    lines = []
    w = lines.append
    w("// GENERATED FILE -- DO NOT EDIT.")
    w("//")
    w("// Written by src/gen_protocol_header.py from config/protocol.json.  It carries the")
    w("// parts of the frozen protocol that the probe itself has to act on: the quantile")
    w("// levels and their order-statistic index brackets, the upper-tail thresholds, the")
    w("// incomplete-beta switch point, the frozen sample sizes and the matrix.  Everything")
    w("// else the protocol fixes is a decision rule, and decision rules belong to the")
    w("// analysis, which reads the JSON directly.")
    w("//")
    w(f"// config/protocol.json SHA-256: {sha}")
    w("")
    w("#ifndef EXP7_PROTOCOL_H")
    w("#define EXP7_PROTOCOL_H")
    w("")
    w("namespace exp7")
    w("{")
    w("namespace protocol")
    w("{")
    w("")
    w(f'static const char *const kProtocolVersion= "{p["protocol_version"]}";')
    w(f'static const char *const kProtocolSha256= "{sha}";')
    w("")
    w("/// Intended attempts per seed, by phase, exactly as PROTOCOL.md Sec. 4 freezes them.")
    w(f"static const long long kNScalar= {int(matrix['n_scalar'])};")
    w(f"static const long long kNMechanism= {int(matrix['n_mechanism'])};")
    w(f"static const long long kNConditioning= {int(matrix['n_conditioning'])};")
    w(f"static const long long kNLoader= {int(matrix['n_loader'])};")
    w("")
    w("/// Quantile levels of F2/F3, and the exact 0-based order-statistic bracket the")
    w("/// protocol froze for each of them at n = kNScalar.  The probe emits the sample")
    w("/// values at these indices; the analysis owns the coverage rule.")
    w(f"static const int kNumQuantileLevels= {len(levels)};")
    w("static const double kQuantileLevels[kNumQuantileLevels]= {")
    w("    " + ", ".join(c_double(lv) for lv in levels) + "};")
    w("static const long long kQuantileLoIndex[kNumQuantileLevels]= {")
    w("    " + ", ".join(str(brackets[lv][0]) for lv in levels) + "};")
    w("static const long long kQuantileHiIndex[kNumQuantileLevels]= {")
    w("    " + ", ".join(str(brackets[lv][1]) for lv in levels) + "};")
    w("static const double kQuantileAchievedCoverage[kNumQuantileLevels]= {")
    w("    " + ", ".join(c_double(brackets[lv][2]) for lv in levels) + "};")
    w("")
    w("/// F4 upper-tail thresholds: exceedance counts of Z above z0 = -log(q0).")
    tail = list(matrix["tail_q0"])
    w(f"static const int kNumTailQ0= {len(tail)};")
    w("static const double kTailQ0[kNumTailQ0]= {")
    w("    " + ", ".join(c_double(q) for q in tail) + "};")
    w("")
    w("/// Switch point of the two-piece log I_W(a, 3/2) evaluator, in log W.")
    w(f"static const double kLogQSwitch= {c_double(stats['log_q_switch'])};")
    w("")
    w("/// Amendment 1.1.0: a returned draw whose radius differs from the exact value by")
    w("/// more than this fraction is FINITE_BUT_WRONG -- a loss, not a success.  The")
    w("/// threshold is 2^-(digits/2), i.e. more than half the significand of its type gone.")
    acc = p["accuracy"]["max_relative_error"]
    w(f"static const double kMaxRelErrorDouble= {c_double(acc['double'])};")
    w(f"static const double kMaxRelErrorFloat= {c_double(acc['float'])};")
    w("")
    w("/// Amendment 1.1.0: the audit stratum sampling rates are a protocol datum, not an")
    w("/// implementation detail.  A rate of 1 means the stratum is adjudicated in full.")
    rates = p["audit"]["strata_rates"]
    order = ["near_limit_2_log_units", "method_disagreement", "finite_but_wrong_candidate",
             "public_path_failure", "subnormal_denominator", "uniform_sample"]
    missing = [k for k in order if k not in rates]
    if missing:
        sys.stderr.write(f"exp7: protocol.json declares no rate for {missing}\n")
        return 1
    extra = [k for k in rates if k not in order]
    if extra:
        sys.stderr.write(f"exp7: protocol.json declares audit strata this build does not "
                         f"implement: {extra}\n")
        return 1
    w(f"static const int kNumAuditRates= {len(order)};")
    w("static const char *const kAuditStratumNames[kNumAuditRates]= {")
    w("    " + ", ".join('"%s"' % k for k in order) + "};")
    w("static const double kAuditStratumRates[kNumAuditRates]= {")
    w("    " + ", ".join(c_double(rates[k]) for k in order) + "};")
    w("")
    w("/// The matrix, restated here so that the hand-written declarations in")
    w("/// exp7_common.H can be checked against the frozen JSON at start-up.")
    ladder = list(matrix["kappa_ladder"])
    w(f"static const int kNumKappa= {len(ladder)};")
    w("static const double kKappaLadder[kNumKappa]= {")
    w("    " + ", ".join(c_double(k) for k in ladder) + "};")
    prod = list(p["seeds"]["production"])
    perf = list(p["seeds"]["performance"])
    w(f"static const int kNumProductionSeeds= {len(prod)};")
    w("static const unsigned kProductionSeeds[kNumProductionSeeds]= {")
    w("    " + ", ".join(str(int(s)) + "u" for s in prod) + "};")
    w(f"static const int kNumPerformanceSeeds= {len(perf)};")
    w("static const unsigned kPerformanceSeeds[kNumPerformanceSeeds]= {")
    w("    " + ", ".join(str(int(s)) + "u" for s in perf) + "};")
    w("")
    w("} // namespace protocol")
    w("} // namespace exp7")
    w("")
    w("#endif // EXP7_PROTOCOL_H")
    w("")

    text = "\n".join(lines)
    old = None
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as fh:
            old = fh.read()
    if old != text:
        with open(OUT, "w", encoding="utf-8") as fh:
            fh.write(text)
    print(f"  src/exp7_protocol.H from config/protocol.json ({sha[:12]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
