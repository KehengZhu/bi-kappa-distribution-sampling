#!/usr/bin/env python3
"""Compare regenerated files with the committed ones, for ``make reverify`` and
``make protocol-check``.

Run with:  uv run --project ../../python python compare_regenerated.py outputs|protocol [--strict]

By default two text files agree when everything but their floating-point numbers is
identical and every pair of floating-point numbers agrees to a relative tolerance
(``--rtol``, default 1e-9).  Integers, words, keys, column names and row counts must match
exactly, so every count, verdict and decision is still compared exactly.  A different
numpy or scipy changes the last digits of some p-values and quantiles, which this accepts;
a changed result does not pass.  ``--strict`` asks for byte-for-byte equality instead, which
is what the committed files satisfy under the package versions recorded in
``raw/environment.json``.

``outputs`` compares ``.reverify/results_before`` with ``.reverify/results_after`` and
``.reverify/figures_before/figure_manifest.json`` with its regenerated copy.  Rendered
figures are compared byte for byte and only reported, as before.

``protocol`` checks ``config/protocol.json`` against ``config/make_protocol.py`` without
overwriting it.  The only part of the protocol that depends on a random stream is the Monte
Carlo estimate of the F2 miss-count distribution, and numpy does not keep the stream of
``Generator.multinomial`` fixed across versions.  So the check has two parts:

  1. the protocol is regenerated with the committed Monte Carlo counts substituted, and
     must agree with the committed file under the comparison above.  Every value derived
     from the counts -- the critical miss count, the attained level, the stability rows --
     is therefore recomputed and checked;
  2. the Monte Carlo is re-run, and its counts must be consistent with the committed ones
     (chi-square test of homogeneity, rejected below ``--mc-alpha``).

Under ``--strict`` the protocol is regenerated with a fresh Monte Carlo and compared byte
for byte, as ``make protocol-check`` did before.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# A SHA-256 digest, and a number that is not part of a longer word.  Digests are taken out
# first so that the digits inside them are never read as numbers.
DIGEST = re.compile(r"(?<![0-9A-Za-z])[0-9a-f]{64}(?![0-9A-Za-z])")
NUMBER = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![\w.])")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def split_text(text: str) -> tuple[str, list[str], list[str]]:
    """The text with digests and numbers replaced by placeholders, and the two lists."""
    digests = DIGEST.findall(text)
    text = DIGEST.sub("\x00D", text)
    numbers = NUMBER.findall(text)
    return NUMBER.sub("\x00N", text), digests, numbers


def is_integer_token(tok: str) -> bool:
    return re.fullmatch(r"[-+]?\d+", tok) is not None


class Result:
    """What one comparison found."""

    def __init__(self, name: str) -> None:
        self.name = name
        self.problems: list[str] = []
        self.n_float_diff = 0
        self.n_digest_diff = 0
        self.max_rel = 0.0
        self.identical = False

    @property
    def ok(self) -> bool:
        return not self.problems

    def line(self) -> str:
        if self.identical:
            return f"  {self.name}: identical"
        if self.ok:
            parts = []
            if self.n_float_diff:
                parts.append(f"{self.n_float_diff} numbers differ, max relative difference "
                             f"{self.max_rel:.1e}")
            if self.n_digest_diff:
                parts.append(f"{self.n_digest_diff} digests of regenerated files differ")
            return f"  {self.name}: agrees within tolerance ({'; '.join(parts)})"
        shown = "\n".join("      " + p for p in self.problems[:8])
        more = len(self.problems) - 8
        return (f"  {self.name}: DIFFERS\n{shown}"
                + (f"\n      ... and {more} more" if more > 0 else ""))


def compare_text(name: str, a: str, b: str, rtol: float, atol: float,
                 digest_pairs: set[tuple[str, str]]) -> Result:
    res = Result(name)
    ska, da, na = split_text(a)
    skb, db, nb = split_text(b)
    if ska != skb:
        la, lb = ska.splitlines(), skb.splitlines()
        for i, (x, y) in enumerate(zip(la, lb)):
            if x != y:
                res.problems.append(f"text differs at line {i + 1}: "
                                    f"{x.replace(chr(0), '#')[:100]!r} vs "
                                    f"{y.replace(chr(0), '#')[:100]!r}")
                break
        else:
            res.problems.append(f"line count differs: {len(la)} vs {len(lb)}")
        return res
    for x, y in zip(da, db):
        # A digest may differ only when it is the digest of a regenerated file, and that
        # file is compared in its own right.
        if x == y:
            continue
        if (x, y) in digest_pairs:
            res.n_digest_diff += 1
        else:
            res.problems.append(f"digest differs: {x[:16]}... vs {y[:16]}...")
    for x, y in zip(na, nb):
        if x == y:
            continue
        if is_integer_token(x) and is_integer_token(y):
            res.problems.append(f"integer differs: {x} vs {y}")
            continue
        fx, fy = float(x), float(y)
        if math.isclose(fx, fy, rel_tol=rtol, abs_tol=atol):
            res.n_float_diff += 1
            scale = max(abs(fx), abs(fy))
            if scale > 0:
                res.max_rel = max(res.max_rel, abs(fx - fy) / scale)
        else:
            res.problems.append(f"number differs beyond tolerance: {x} vs {y}")
    return res


def compare_file(name: str, pa: str, pb: str, strict: bool, rtol: float, atol: float,
                 digest_pairs: set[tuple[str, str]]) -> Result:
    with open(pa, "rb") as fh:
        ba = fh.read()
    with open(pb, "rb") as fh:
        bb = fh.read()
    if ba == bb:
        res = Result(name)
        res.identical = True
        return res
    if strict:
        res = Result(name)
        res.problems.append("bytes differ (--strict)")
        return res
    try:
        ta, tb = ba.decode("utf-8"), bb.decode("utf-8")
    except UnicodeDecodeError:
        res = Result(name)
        res.problems.append("binary file, bytes differ")
        return res
    return compare_text(name, ta, tb, rtol, atol, digest_pairs)


def list_files(root: str) -> set[str]:
    out = set()
    for d, _, files in os.walk(root):
        for f in files:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def digest_pairs_of(pairs: list[tuple[str, str]]) -> set[tuple[str, str]]:
    """(old digest, new digest) of every regenerated file whose bytes changed."""
    out = set()
    for pa, pb in pairs:
        if os.path.isfile(pa) and os.path.isfile(pb):
            with open(pa, "rb") as fa, open(pb, "rb") as fb:
                ha, hb = sha256_bytes(fa.read()), sha256_bytes(fb.read())
            if ha != hb:
                out.add((ha, hb))
    return out


def run_outputs(args: argparse.Namespace) -> int:
    base = os.path.join(HERE, ".reverify")
    rb, ra = os.path.join(base, "results_before"), os.path.join(base, "results_after")
    fb, fa = os.path.join(base, "figures_before"), os.path.join(base, "figures_after")
    ok = True

    before, after = list_files(rb), list_files(ra)
    for f in sorted(before - after):
        print(f"  results/{f}: not regenerated")
        ok = False
    for f in sorted(after - before):
        print(f"  results/{f}: new file, not in the committed results")
        ok = False
    common = sorted(before & after)
    pairs = [(os.path.join(rb, f), os.path.join(ra, f)) for f in common]
    pairs += [(os.path.join(fb, f), os.path.join(fa, f)) for f in list_files(fb) & list_files(fa)]
    digests = digest_pairs_of(pairs)

    n_identical = 0
    for f in common:
        res = compare_file(f"results/{f}", os.path.join(rb, f), os.path.join(ra, f),
                           args.strict, args.rtol, args.atol, digests)
        if res.identical:
            n_identical += 1
        else:
            print(res.line())
        ok &= res.ok
    print(f"  reverify: {n_identical} of {len(common)} files in results/ are byte-identical")

    res = compare_file("figures/figure_manifest.json",
                       os.path.join(fb, "figure_manifest.json"),
                       os.path.join(fa, "figure_manifest.json"),
                       args.strict, args.rtol, args.atol, digests)
    print(res.line())
    ok &= res.ok

    # Rendered figures: reported, never a failure.
    rendered = sorted(f for f in list_files(fb) | list_files(fa) if f != "figure_manifest.json")
    changed = [f for f in rendered
               if not (os.path.isfile(os.path.join(fb, f)) and os.path.isfile(os.path.join(fa, f)))
               or open(os.path.join(fb, f), "rb").read() != open(os.path.join(fa, f), "rb").read()]
    if changed:
        print("  reverify: NOTE -- the bytes of these figure files differ; not a failure:")
        for f in changed[:8]:
            print(f"    figures/{f}")
    else:
        print("  reverify: rendered figures also reproduced byte for byte")

    mode = "byte for byte" if args.strict else f"within relative tolerance {args.rtol:g}"
    if ok:
        print(f"  reverify: every derived artifact regenerated {mode}")
        return 0
    print(f"  reverify: a derived artifact did not regenerate {mode}")
    return 1


def load_make_protocol():
    path = os.path.join(HERE, "config", "make_protocol.py")
    spec = importlib.util.spec_from_file_location("make_protocol", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def regenerate_protocol(mod, outdir: str, counts=None) -> str:
    """Run make_protocol.main() with its output directed to ``outdir``.  With ``counts``,
    the Monte Carlo is replaced by those counts."""
    import numpy as np

    mod.HERE = outdir
    if counts is not None:
        c = np.asarray(counts, dtype=np.int64)
        mod.group_miss_pmf = lambda: (c / float(mod.F2_GROUP_MC_REPLICATES), c)
    with contextlib.redirect_stdout(io.StringIO()):
        mod.main()
    return os.path.join(outdir, "protocol.json")


def homogeneity_pvalue(c1: list[int], c2: list[int]) -> float:
    """Chi-square test that two count vectors come from one distribution.  Trailing cells
    are pooled until every expected count is at least 5."""
    import numpy as np
    from scipy import stats

    t = np.array([c1, c2], dtype=float)
    while t.shape[1] > 2:
        expected = t.sum(axis=0) * t.sum(axis=1)[:, None] / t.sum()
        if expected[:, -1].min() >= 5:
            break
        t[:, -2] += t[:, -1]
        t = t[:, :-1]
    return float(stats.chi2_contingency(t, correction=False)[1])


def run_protocol(args: argparse.Namespace) -> int:
    committed_path = os.path.join(HERE, "config", "protocol.json")
    with open(committed_path) as fh:
        committed = json.load(fh)
    mod = load_make_protocol()

    with tempfile.TemporaryDirectory() as tmp:
        if args.strict:
            fresh = regenerate_protocol(mod, tmp)
            res = compare_file("config/protocol.json", committed_path, fresh, True,
                               args.rtol, args.atol, set())
            print(res.line())
            if res.ok:
                print("  protocol-check: config/protocol.json matches config/make_protocol.py "
                      "byte for byte")
                return 0
            print("  protocol-check: config/protocol.json is NOT what config/make_protocol.py "
                  "produces")
            return 1

        counts = committed["F2"]["group_miss_counts"]
        if committed["F2"]["monte_carlo"]["replicates"] != mod.F2_GROUP_MC_REPLICATES:
            print("  protocol-check: the committed Monte Carlo size is not the script's")
            return 1
        regen = regenerate_protocol(mod, tmp, counts)
        res = compare_file("config/protocol.json", committed_path, regen, False,
                           args.rtol, args.atol, set())
        print(res.line())

        mod = load_make_protocol()
        _, fresh_counts = mod.group_miss_pmf()
        fresh_counts = [int(x) for x in fresh_counts]
        if fresh_counts == counts:
            print("  F2 Monte Carlo: the re-run reproduces the committed counts exactly")
            mc_ok = True
        else:
            p = homogeneity_pvalue(counts, fresh_counts)
            mc_ok = p >= args.mc_alpha
            print(f"  F2 Monte Carlo: the re-run gives other counts ({fresh_counts}); "
                  f"homogeneity with the committed ones p = {p:.4g} "
                  f"({'consistent' if mc_ok else 'INCONSISTENT'} at {args.mc_alpha:g})")

    if res.ok and mc_ok:
        print("  protocol-check: config/protocol.json matches config/make_protocol.py")
        return 0
    print("  protocol-check: config/protocol.json is NOT what config/make_protocol.py produces")
    return 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", choices=["outputs", "protocol"])
    ap.add_argument("--strict", action="store_true", help="require byte-for-byte equality")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--atol", type=float, default=1e-15)
    ap.add_argument("--mc-alpha", type=float, default=1e-4)
    args = ap.parse_args()
    return run_outputs(args) if args.what == "outputs" else run_protocol(args)


if __name__ == "__main__":
    sys.exit(main())
