# P5 portability: the environments, and the exact command for each

Gate G4 asks two different questions of the same phase, and they need different hardware.

On one architecture the candidate's stream is a function of the engine alone — the Gamma,
normal and uniform primitives are in the header rather than taken from `<random>` — so
libc++ and libstdc++ are predicted to agree **bitwise**, with no tolerance. Across
architectures `log`, `exp` and `pow` are not correctly rounded and the arm64 and x86_64
libm implementations differ, so only the **rates** are predicted to agree, inside the
pre-registered equivalence margin of ±0.15 on the log rate ratio (PROTOCOL.md §5.3).

The second question needs an x86_64 machine. This project is developed on arm64 macOS.
Rosetta 2 and Docker/Rosetta will run an x86_64 binary there, and those runs are worth
having, but PROTOCOL.md §5.3 records them as `execution = translated` and excludes them
from the F7 decision. They corroborate; they do not close the gate. Native x86_64 is
reachable here only through CI runners, which is what `.github/workflows/portability.yml`
is for.

Everything below runs the **frozen** sizes: 10⁶ attempts × 5 seeds, production block
7001–7005, over the full κ ladder. The phase costs about a minute of single-core time, so
there is no reason to reduce it and no provision for doing so. A run at any other size or
seed block is not P5 and must not be filed as P5 evidence.

An environment that has not run leaves G4 **open**. It is never recorded as a pass and
never replaced by a zero. That is what the `available` / `completed` / `reason_not_run`
columns in the table below exist to make visible.

---

## The recipe

Run from the repository root. Three values change per environment; everything else is
identical, which is the point — the comparison is only meaningful if the environments
differ in exactly the way being measured.

```bash
ENV_ID="macos-arm64"                      # see the per-environment table below
CXX_LIBCXX="clang++ -stdlib=libc++"       # empty string if this environment has no libc++
CXX_LIBSTDCXX="g++-15"                    # empty string if it has no libstdc++

mkdir -p evidence

# Build the probe and run P5.  The sizes come from the makefile's frozen defaults;
# do not pass N_PORTABILITY or SEEDS_PROD.  The build log is kept because the
# environment record quotes the compile line make actually used and fails if the
# frozen flags (-std=c++11 -O2 -ffp-contract=off) are not on it.
make -C experiments/exp7_confirmatory portability \
  CXX_LIBCXX="$CXX_LIBCXX" \
  CXX_LIBSTDCXX="$CXX_LIBSTDCXX" \
  2>&1 | tee evidence/run.log

# One environment record per standard library.
for tag in libcxx libstdcxx; do
  case "$tag" in
    libcxx)    cxx="$CXX_LIBCXX" ;;
    libstdcxx) cxx="$CXX_LIBSTDCXX" ;;
  esac
  [ -n "$cxx" ] || continue
  $cxx -std=c++11 -O2 -ffp-contract=off -Wall -Wextra \
       -I cpp -o "/tmp/env_probe_$tag" .github/ci/env_probe.cpp
  "/tmp/env_probe_$tag" > "/tmp/env_probe_$tag.json"
  python3 .github/ci/make_env_record.py \
    --probe-json "/tmp/env_probe_$tag.json" \
    --tag "$tag" \
    --env-id "$ENV_ID" \
    --cxx "$cxx" \
    --repo-root "$PWD" \
    --build-log evidence/run.log \
    --phase P5 \
    --out "evidence/environment_$tag.json"
  cp "experiments/exp7_confirmatory/raw/portability/portability_$tag.jsonl" evidence/
done

# Package it the way the workflow's artifact is packaged, so that a hand run and a CI
# run are ingested by the same code path.
mkdir -p evidence-bundle/"portability-$ENV_ID"
cp evidence/* evidence-bundle/"portability-$ENV_ID"/
```

`execution` is **computed**, not declared: `make_env_record.py` compares the architecture
the binary was compiled for with the architecture the host reports, and consults
`sysctl sysctl.proc_translated` on macOS. A run under Rosetta is recorded as
`translated` whether or not whoever started it noticed.

Ship `evidence-bundle/` to whichever machine runs the decision, merge it with the other
environments' bundles into one directory, and run:

```bash
python3 .github/ci/check_evidence.py \
  --evidence evidence-bundle \
  --expect .github/ci/expected_environments.json \
  --protocol experiments/exp7_confirmatory/config/protocol.json \
  --out evidence-bundle/coverage.json

cd experiments/exp7_confirmatory
uv run --project ../../python python analyze.py \
  --portability-ingest ../../evidence-bundle
```

The first command is bookkeeping — which environments are present, whether each ran
native, whether the frozen sizes and seeds were used. It exits non-zero if an environment
is missing, because a silent pass on absent evidence is the specific failure this file
exists to prevent. The second decides F7 and G4 against `config/protocol.json`.

---

## Per-environment settings

### arm64 · macOS · Apple clang/libc++ and Homebrew GCC 15/libstdc++ · **native**

Both toolchains are present on the development host (`brew install gcc` supplies the
second). One run of the recipe produces both records, side by side in one bundle, which
is what the bitwise test needs — and it is the only cross-standard-library pair
obtainable without CI. The `macos-latest` runner has no libstdc++ toolchain, so it
contributes `environment_libcxx.json` alone; a hand run on this host adds
`environment_libstdcxx.json` to the same `portability-macos-arm64` bundle, and
`analyze.py` then has an arm64 pair to test bitwise even if the `ubuntu-24.04-arm` leg
never runs.

```bash
ENV_ID="macos-arm64"
CXX_LIBCXX="clang++ -stdlib=libc++"
CXX_LIBSTDCXX="g++-15"
```

### x86_64 · Ubuntu 24.04 · clang/libc++ and g++/libstdc++ · **native**

Not available locally. Reached through the `ubuntu-latest` runner in
`.github/workflows/portability.yml`, or on any native x86_64 Linux machine after

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends build-essential clang libc++-dev libc++abi-dev
```

```bash
ENV_ID="linux-x86_64"
CXX_LIBCXX="clang++ -stdlib=libc++"
CXX_LIBSTDCXX="g++"
```

This environment carries the cross-architecture half of G4. Without it the gate cannot
close, however many arm64 environments pass.

### arm64 · Ubuntu 24.04 · clang/libc++ and g++/libstdc++ · **native**

Reached through the `ubuntu-24.04-arm` runner. That label's availability depends on the
account's plan and on repository visibility, so the job is marked `continue-on-error` and
the environment is listed here as pending rather than assumed. The same commands as the
x86_64 Linux row, with

```bash
ENV_ID="linux-arm64"
```

It is the useful control against macOS/arm64: same architecture, different operating
system and libm, so a difference here separates an architecture effect from an OS one.

### x86_64 · macOS via Rosetta 2 · Apple clang · libc++ · **translated**

Available locally; corroborating only.

```bash
ENV_ID="macos-x86_64-rosetta"
CXX_LIBCXX="arch -x86_64 clang++ -stdlib=libc++ -arch x86_64"
CXX_LIBSTDCXX=""
```

`make_env_record.py` will record `execution = translated`, `check_evidence.py` will refuse
to count it toward coverage, and `analyze.py` drops it before the F7 test. Run it for
information; do not file it as closure.

### x86_64 · Linux in Docker via Rosetta · g++ · libstdc++ · **translated**

Available locally; corroborating only.

```bash
docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w ubuntu:24.04 bash -lc '
  apt-get update &&
  apt-get install -y --no-install-recommends build-essential clang libc++-dev libc++abi-dev python3 &&
  make -C experiments/exp7_confirmatory portability CXX_LIBCXX="clang++ -stdlib=libc++" CXX_LIBSTDCXX="g++"'
```

Same standing as the row above.

---

## Status

`execution` is the mode the run was or would be performed in. `available` says whether the
environment can be reached from this project as it stands. `completed` says whether P5 has
actually been run there and its evidence filed.

| arch | os | compiler | stdlib | execution | available | completed | reason_not_run |
|---|---|---|---|---|---|---|---|
| arm64 | macOS 26.7 (dev host) | Apple clang 21.0.0 | libc++ | native | yes | no | P5 not yet run; environment present and the command is above |
| arm64 | macOS 26.7 (dev host) | Homebrew GCC 15.2.0 | libstdc++ | native | yes | no | P5 not yet run; environment present and the command is above |
| x86_64 | Ubuntu 24.04 (`ubuntu-latest`) | clang 18 | libc++ | native | pending | no | CI never executed: the branch has not been pushed |
| x86_64 | Ubuntu 24.04 (`ubuntu-latest`) | g++ 13 | libstdc++ | native | pending | no | CI never executed: the branch has not been pushed |
| arm64 | Ubuntu 24.04 (`ubuntu-24.04-arm`) | clang 18 | libc++ | native | pending | no | CI never executed; runner label availability not yet confirmed for this account |
| arm64 | Ubuntu 24.04 (`ubuntu-24.04-arm`) | g++ 13 | libstdc++ | native | pending | no | CI never executed; runner label availability not yet confirmed for this account |
| arm64 | macOS (`macos-latest`) | Apple clang | libc++ | native | pending | no | CI never executed: the branch has not been pushed |
| x86_64 | macOS 26.7 via Rosetta 2 | Apple clang 21.0.0 | libc++ | translated | yes | no | excluded from the F7 decision by PROTOCOL.md §5.3; corroborating only |
| x86_64 | Ubuntu 24.04 in Docker/Rosetta | g++ 13 | libstdc++ | translated | yes | no | excluded from the F7 decision by PROTOCOL.md §5.3; corroborating only |

**G4 is open.** No environment has completed P5. The compiler versions given for the CI
rows are the current defaults of those runner images and are recorded as expectations; the
environment record each job uploads carries the versions actually used, and that record,
not this table, is the evidence.
