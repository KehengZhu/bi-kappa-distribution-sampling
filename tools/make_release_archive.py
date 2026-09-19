#!/usr/bin/env python3
"""Build a bit-reproducible release archive, and a manifest that can verify it.

Two people running this on the same source revision must get byte-identical archives, or the
checksum in a paper is a checksum of one machine's filesystem rather than of the release.
`tar` alone does not give that: it records mtimes, uids, gids, usernames and directory order,
and GNU tar's `--sort=name` is not available in the BSD tar that ships with macOS. So the
archive is assembled here instead, with every one of those fields pinned.

What goes in:

  * every file tracked by git at the chosen revision -- that is the source of record, and
    using `git ls-tree` rather than the working tree means an uncommitted edit cannot leak in;
  * optionally, the experiment's untracked bulk results (`--with-results`), which are too
    large to track but are covered by their own checksum manifest.

Usage:
    tools/make_release_archive.py --ref v2.0.0 --out dist/
    tools/make_release_archive.py --ref HEAD --with-results experiments/exp7_confirmatory
"""
from __future__ import annotations

import argparse
import hashlib
import io
import os
import subprocess
import sys
import tarfile

# Fixed epoch for every entry.  Any constant works; what matters is that it is a constant.
FIXED_MTIME = 1_000_000_000  # 2001-09-09T01:46:40Z


def run(*args: str, cwd: str | None = None) -> str:
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True,
                          text=True).stdout.rstrip("\n")


def tracked_files(repo: str, ref: str) -> list[str]:
    out = run("git", "-C", repo, "ls-tree", "-r", "--name-only", ref)
    return sorted(p for p in out.split("\n") if p)


def blob(repo: str, ref: str, path: str) -> bytes:
    return subprocess.run(["git", "-C", repo, "show", f"{ref}:{path}"],
                          check=True, capture_output=True).stdout


def add(tar: tarfile.TarFile, name: str, data: bytes, mode: int = 0o644) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mtime = FIXED_MTIME
    info.mode = mode
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.type = tarfile.REGTYPE
    tar.addfile(info, io.BytesIO(data))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD", help="git revision to archive")
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    ap.add_argument("--out", default="dist")
    ap.add_argument("--name", default=None, help="archive basename (default bi-kappa-<ref>)")
    ap.add_argument("--with-results", action="append", default=[],
                    metavar="DIR",
                    help="also include this directory's untracked result files "
                         "(repeatable); they are read from the working tree, so the tree "
                         "must be the one that produced them")
    args = ap.parse_args()

    repo = os.path.abspath(args.repo)
    commit = run("git", "-C", repo, "rev-parse", args.ref)
    short = commit[:12]
    name = args.name or f"bi-kappa-{args.ref.replace('/', '-')}-{short}"
    outdir = os.path.abspath(args.out)
    os.makedirs(outdir, exist_ok=True)
    tarpath = os.path.join(outdir, f"{name}.tar")

    files = tracked_files(repo, args.ref)
    entries: list[tuple[str, bytes]] = [(p, blob(repo, args.ref, p)) for p in files]

    # Untracked bulk results, read from the working tree and listed explicitly so that what
    # is included is visible in the manifest rather than implied by a glob.
    for d in args.with_results:
        root = os.path.join(repo, d)
        if not os.path.isdir(root):
            print(f"error: {d} is not a directory", file=sys.stderr)
            return 2
        extra = []
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = sorted(n for n in dirnames if n not in ("__pycache__", ".git"))
            for fn in sorted(filenames):
                full = os.path.join(dirpath, fn)
                rel = os.path.relpath(full, repo)
                if rel in files or fn.endswith((".part", ".exe", ".pyc")):
                    continue
                extra.append(rel)
        for rel in sorted(extra):
            with open(os.path.join(repo, rel), "rb") as fh:
                entries.append((rel, fh.read()))

    entries.sort(key=lambda e: e[0])

    # A manifest of what is inside, hashed, and the provenance of the archive itself.
    lines = [f"# bi-kappa release archive manifest",
             f"# source-revision {commit}",
             f"# ref {args.ref}",
             f"# entries {len(entries)}",
             f"# fixed-mtime {FIXED_MTIME}"]
    for rel, data in entries:
        lines.append(f"{hashlib.sha256(data).hexdigest()}  {rel}")
    manifest = ("\n".join(lines) + "\n").encode()

    with tarfile.open(tarpath, "w", format=tarfile.PAX_FORMAT) as tar:
        add(tar, f"{name}/ARCHIVE_MANIFEST.sha256", manifest)
        for rel, data in entries:
            add(tar, f"{name}/{rel}", data)

    # gzip separately with mtime=0, since GzipFile otherwise stamps the current time.
    import gzip
    gzpath = tarpath + ".gz"
    with open(tarpath, "rb") as src, open(gzpath, "wb") as dst:
        with gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0) as gz:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                gz.write(chunk)
    os.remove(tarpath)

    digest = hashlib.sha256(open(gzpath, "rb").read()).hexdigest()
    with open(gzpath + ".sha256", "w") as fh:
        fh.write(f"{digest}  {os.path.basename(gzpath)}\n")

    print(f"archive      {gzpath}")
    print(f"sha256       {digest}")
    print(f"revision     {commit}")
    print(f"entries      {len(entries)}")
    print(f"size         {os.path.getsize(gzpath):,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
