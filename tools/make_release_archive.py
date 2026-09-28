#!/usr/bin/env python3
"""Build a bit-reproducible release archive, and a manifest that can verify it.

Two people running this on the same source revision must get byte-identical archives, or a
published checksum is a checksum of one machine's filesystem rather than of the release.
`tar` alone does not give that: it records mtimes, uids, gids, usernames and directory order,
and GNU tar's `--sort=name` is not available in the BSD tar that ships with macOS. So the
archive is assembled here instead, with every one of those fields pinned.

What goes in: exactly the files `git archive <ref>` exports, which is also what the GitHub
source archive and the Zenodo deposit of a tag contain. That is every file tracked at the
chosen revision except those marked `export-ignore` in `.gitattributes` (the release tools,
the CI machinery of the finite-precision study, and the parts of that study that stay in the
repository). Reading from the revision rather than the working tree means an uncommitted
edit cannot leak in. Each file keeps git's mode: 0755 if it is executable, 0644 otherwise.

The archive holds one extra file, `ARCHIVE_MANIFEST.sha256`, listing the SHA-256 of every
entry; `tools/verify_release.py` checks an archive against it.

Usage:
    tools/make_release_archive.py --ref v3.0.0 --out dist/
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import subprocess
import tarfile

# Fixed epoch for every entry.  Any constant works; what matters is that it is a constant.
FIXED_MTIME = 1_000_000_000  # 2001-09-09T01:46:40Z


def run(*args: str, cwd: str | None = None) -> str:
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True,
                          text=True).stdout.rstrip("\n")


def exported_files(repo: str, ref: str) -> list[tuple[str, bytes, int]]:
    """The regular files `git archive` exports at `ref`, as (path, bytes, mode)."""
    raw = subprocess.run(["git", "-C", repo, "archive", "--format=tar", ref],
                         check=True, capture_output=True).stdout
    entries = []
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as tar:
        for info in tar:
            if info.isdir():
                continue
            if not info.isreg():
                raise SystemExit(f"error: {info.name} is not a regular file; "
                                 "this script handles regular files only")
            data = tar.extractfile(info).read()
            mode = 0o755 if info.mode & 0o100 else 0o644
            entries.append((info.name, data, mode))
    return sorted(entries)


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
    ap = argparse.ArgumentParser(
        description="Build a byte-reproducible .tar.gz of what `git archive` exports.")
    ap.add_argument("--ref", default="HEAD", help="git revision to archive (default HEAD)")
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), help="repository root (default: this script's repo)")
    ap.add_argument("--out", default="dist", help="output directory (default dist/)")
    ap.add_argument("--name", default=None,
                    help="archive basename (default bi-kappa-<ref>-<first 12 hex of commit>)")
    args = ap.parse_args()

    repo = os.path.abspath(args.repo)
    commit = run("git", "-C", repo, "rev-parse", f"{args.ref}^{{commit}}")
    short = commit[:12]
    name = args.name or f"bi-kappa-{args.ref.replace('/', '-')}-{short}"
    outdir = os.path.abspath(args.out)
    os.makedirs(outdir, exist_ok=True)
    tarpath = os.path.join(outdir, f"{name}.tar")

    entries = exported_files(repo, commit)

    # A manifest of what is inside, hashed, and the provenance of the archive itself.
    lines = ["# bi-kappa release archive manifest",
             f"# source-revision {commit}",
             f"# ref {args.ref}",
             f"# entries {len(entries)}",
             f"# fixed-mtime {FIXED_MTIME}"]
    for rel, data, _ in entries:
        lines.append(f"{hashlib.sha256(data).hexdigest()}  {rel}")
    manifest = ("\n".join(lines) + "\n").encode()

    with tarfile.open(tarpath, "w", format=tarfile.PAX_FORMAT) as tar:
        add(tar, f"{name}/ARCHIVE_MANIFEST.sha256", manifest)
        for rel, data, mode in entries:
            add(tar, f"{name}/{rel}", data, mode)

    # gzip separately with mtime=0, since GzipFile otherwise stamps the current time.
    gzpath = tarpath + ".gz"
    with open(tarpath, "rb") as src, open(gzpath, "wb") as dst:
        with gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0) as gz:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                gz.write(chunk)
    os.remove(tarpath)

    with open(gzpath, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
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
