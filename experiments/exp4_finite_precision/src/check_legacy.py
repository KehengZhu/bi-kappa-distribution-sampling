"""Verify that the vendored comparator is still the released 1.0.0 header.

Run with:  uv run --project ../../python python src/check_legacy.py

``src/legacy/bi_kappa_distribution_v1.H`` is ``cpp/bi_kappa_distribution.H`` exactly as
released at version 1.0.0, with two changes that are forced by having both headers in one
translation unit: the include guard is renamed and the whole file is wrapped in
``namespace bikappa_v1``.  A banner recording that fact was added above the wrap.

This script undoes precisely those three edits and compares the result, byte for byte,
against ``git show <commit>:cpp/bi_kappa_distribution.H``.  Anything else that has been
changed in the comparator -- a "harmless" reformatting, a compiler warning silenced, a
constant nudged -- shows up as a diff and fails the check.  A comparator that has drifted
is not a comparator.
"""

from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
VENDORED = os.path.join(HERE, "legacy", "bi_kappa_distribution_v1.H")

# The commit whose cpp/bi_kappa_distribution.H is release 1.0.0.
RELEASE_COMMIT = "0fc2c95"
RELEASE_PATH = "cpp/bi_kappa_distribution.H"

GUARD_VENDORED = "_BI_KAPPA_DISTRIBUTION_V1_H_"
GUARD_ORIGINAL = "_BI_KAPPA_DISTRIBUTION_H_"
NAMESPACE_OPEN = "namespace bikappa_v1"
NAMESPACE_CLOSE = "} // namespace bikappa_v1"
BANNER_MARK = "FROZEN COMPARATOR -- DO NOT EDIT."
RULE = "// " + "-" * 75


def unwrap(text: str) -> tuple[str, list[str]]:
    """Undo the guard rename and the namespace wrap.  Returns (text, notes)."""
    notes = []
    lines = text.split("\n")

    # 1. The banner: a comment block whose rule lines bracket the FROZEN COMPARATOR mark.
    mark = [i for i, ln in enumerate(lines) if BANNER_MARK in ln]
    if len(mark) != 1:
        notes.append(f"expected exactly one {BANNER_MARK!r} line, found {len(mark)}")
    else:
        i = mark[0]
        start = i
        while start > 0 and lines[start - 1].startswith("//"):
            start -= 1
        end = i
        while end + 1 < len(lines) and lines[end + 1].startswith("//"):
            end += 1
        if not (lines[start] == RULE and lines[end] == RULE):
            notes.append("the banner is not delimited by the expected rule lines")
        del lines[start:end + 1]
        # The banner is followed by the namespace opening; drop the blank line it left.
        if start < len(lines) and lines[start] == "":
            pass

    # 2. The namespace: its opening line, the brace on the next line, the blank line after
    #    it, and the closing line with the blank line before it.
    try:
        o = lines.index(NAMESPACE_OPEN)
    except ValueError:
        notes.append(f"{NAMESPACE_OPEN!r} not found")
        o = None
    if o is not None:
        span = 1
        if o + span < len(lines) and lines[o + span] == "{":
            span += 1
        else:
            notes.append("the namespace opening is not followed by a bare '{'")
        if o + span < len(lines) and lines[o + span] == "":
            span += 1
        del lines[o:o + span]

    try:
        c = lines.index(NAMESPACE_CLOSE)
    except ValueError:
        notes.append(f"{NAMESPACE_CLOSE!r} not found")
        c = None
    if c is not None:
        start = c
        if start > 0 and lines[start - 1] == "":
            start -= 1
        del lines[start:c + 1]

    out = "\n".join(lines)

    # 3. The include guard.
    if GUARD_VENDORED not in out:
        notes.append(f"{GUARD_VENDORED} not found; the guard was not renamed as expected")
    out = out.replace(GUARD_VENDORED, GUARD_ORIGINAL)
    return out, notes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", default=RELEASE_COMMIT)
    ap.add_argument("--show-diff", action="store_true", default=True)
    args = ap.parse_args()

    proc = subprocess.run(
        ["git", "show", f"{args.commit}:{RELEASE_PATH}"],
        cwd=ROOT, capture_output=True, check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(
            f"check-legacy: cannot read {args.commit}:{RELEASE_PATH}\n"
            f"{proc.stderr.decode('utf-8', 'replace')}")
        return 1
    original = proc.stdout.decode("utf-8")

    with open(VENDORED, encoding="utf-8") as fh:
        vendored = fh.read()

    recovered, notes = unwrap(vendored)
    for note in notes:
        sys.stderr.write(f"check-legacy: {note}\n")

    if recovered == original:
        print(f"  check-legacy: src/legacy/bi_kappa_distribution_v1.H reproduces "
              f"{args.commit}:{RELEASE_PATH} exactly")
        return 1 if notes else 0

    sys.stderr.write(
        "check-legacy: the vendored comparator no longer reproduces the released 1.0.0\n"
        "              header after undoing the guard rename and the namespace wrap.\n")
    if args.show_diff:
        diff = difflib.unified_diff(
            original.split("\n"), recovered.split("\n"),
            fromfile=f"{args.commit}:{RELEASE_PATH}",
            tofile="src/legacy/bi_kappa_distribution_v1.H (unwrapped)", lineterm="")
        for i, line in enumerate(diff):
            if i > 200:
                sys.stderr.write("              ... diff truncated\n")
                break
            sys.stderr.write("              " + line + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
