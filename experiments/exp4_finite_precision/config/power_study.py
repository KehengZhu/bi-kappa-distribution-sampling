#!/usr/bin/env python3
"""Pre-registration power study for the frozen families of ``PROTOCOL.md`` §5.

Run and committed **before any Experiment 7 data existed**, from simulated draws only.  It
changes no decision rule; it states what the frozen rules can and cannot detect, which is a
thing a reader is entitled to know before being shown a result that says "no defect found".

Two questions:

1. **Level.** Does each family reject a correct sampler at its own alpha_F?
2. **Power against survivor conditioning.**  The draws a loader loses near kappa = 1/2 are
   the ones in the far upper tail, so the finite survivors are the target conditioned on
   ``Z < -log q`` where ``q`` is the loss fraction.  That is precisely the effect the battery
   is used to certify *absent* in the candidate, so its power against it at the frozen sample
   size is what decides whether "the candidate passes" means anything.

The answer, in one line: the bulk-weighted statistics are blind below a 1e-3 loss fraction,
and the F4 exceedance tests are not.  Experiment 6 had no F4, which is why its radial battery
detected conditioning only above roughly 1e-3 -- and the cells it certified sat below that.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.stats as st

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import exp7_families as F  # noqa: E402

N = 1_000_000                 # the frozen P1/P3 sample size, per configuration per seed
LOSS_FRACTIONS = [1e-2, 1e-3, 1e-4, 1e-5]
REPLICATES_POWER = 200
REPLICATES_LEVEL = 2000
RNG_SEED = 20260919           # declared, not tuned


def survivors(rng, n: int, q: float) -> np.ndarray:
    """``n`` draws from Exp(1) conditioned on ``Z < -log q`` -- the finite-survivor law."""
    if q <= 0.0:
        return rng.exponential(1.0, n)
    out, need = [], n
    cut = -np.log(q)
    while need > 0:
        z = rng.exponential(1.0, int(need / (1 - q)) + 1000)
        z = z[z < cut]
        out.append(z[:need])
        need -= z[:need].size
    return np.concatenate(out)


def one_cell(proto, z):
    a2, p_ad = F.anderson_darling_exp1(z)
    p_ks = st.kstest(z, "expon").pvalue
    p_cvm = st.cramervonmises(z, "expon").pvalue
    f1 = F.family_F1(proto, [{"label": "c", "pvalues": [p_ad, p_ks, p_cvm]}])
    tests = [F.exceedance_test(z, q0) for q0 in proto.get("matrix", "tail_q0")]
    f4 = F.family_F4(proto, [{"label": "c", "tests": tests}])
    return {"A2": a2, "p_ad": p_ad, "p_ks": p_ks, "p_cvm": p_cvm,
            "F1_rejects": not f1.passed, "F4_rejects": not f4.passed}


def main() -> None:
    proto = F.Protocol()
    rng = np.random.default_rng(RNG_SEED)
    report = {"n_per_configuration": N, "rng_seed": RNG_SEED,
              "protocol_sha256": proto.sha256,
              "alpha_F1": proto.alpha("F1_radial_law"),
              "alpha_F4": proto.alpha("F4_upper_tail_mass"),
              "level": {}, "power": {}}

    # ---- level, at a smaller n where many replicates are affordable ---------------------
    n_level = 50_000
    hits = {"ad": 0, "ks": 0, "cvm": 0, "F1": 0, "F4": 0}
    ps = []
    for _ in range(REPLICATES_LEVEL):
        r = one_cell(proto, rng.exponential(1.0, n_level))
        ps.append(r["p_ad"])
        hits["ad"] += r["p_ad"] < 0.01
        hits["ks"] += r["p_ks"] < 0.01
        hits["cvm"] += r["p_cvm"] < 0.01
        hits["F1"] += r["F1_rejects"]
        hits["F4"] += r["F4_rejects"]
    report["level"] = {k: v / REPLICATES_LEVEL for k, v in hits.items()}
    report["level"]["n"] = n_level
    report["level"]["replicates"] = REPLICATES_LEVEL
    report["level"]["ad_pvalue_uniformity_p"] = float(st.kstest(ps, "uniform").pvalue)

    # ---- power against survivor conditioning at the frozen n ---------------------------
    for q in LOSS_FRACTIONS:
        hits = {"ad": 0, "ks": 0, "cvm": 0, "F1": 0, "F4": 0}
        for _ in range(REPLICATES_POWER):
            r = one_cell(proto, survivors(rng, N, q))
            hits["ad"] += r["p_ad"] < 0.01
            hits["ks"] += r["p_ks"] < 0.01
            hits["cvm"] += r["p_cvm"] < 0.01
            hits["F1"] += r["F1_rejects"]
            hits["F4"] += r["F4_rejects"]
        report["power"][f"{q:.0e}"] = {k: v / REPLICATES_POWER for k, v in hits.items()}
        print(f"  loss {q:.0e}: " + "  ".join(
            f"{k}={v / REPLICATES_POWER:.2f}" for k, v in hits.items()), flush=True)

    report["replicates_power"] = REPLICATES_POWER
    required = float(proto.get("F6", "required_power"))
    report["required_power"] = required
    report["minimum_detectable_loss_fraction"] = {
        fam: next((q for q in sorted(LOSS_FRACTIONS)
                   if report["power"][f"{q:.0e}"][fam] >= required), None)
        for fam in ("F1", "F4")}

    out = os.path.join(HERE, "power_study.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    print(f"wrote {out}")
    print("  minimum detectable loss fraction at power >= "
          f"{required:g}: {report['minimum_detectable_loss_fraction']}")


if __name__ == "__main__":
    main()
