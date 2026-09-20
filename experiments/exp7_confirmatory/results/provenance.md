# Experiment 7 - provenance of this analysis

This is the only file under results/ that carries wall-clock times and build identity. Every value below is copied from raw/manifest.csv and the environment record written at run time; nothing is read from the clock or the host while the analysis runs, so this file regenerates byte for byte and `make reverify` remains a meaningful check.

- protocol: `config/protocol.json` version 3.0.0, SHA-256 `c3b4a0a07877c0f8e2370ff7fdae98d3f3d96c15465aadb9eae00a280391c6a6`
- mode: production
- manifest: `raw/manifest.csv`, 846 distinct output files recorded in 846 manifest rows

- probe run times (UTC, from the writing processes): 2026-09-20T20:47:05Z .. 2026-09-20T20:58:07Z

## Builds that produced the data

| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |
|---|---|---|---|---|---|---|
| libcxx | clang 21.0.0 (clang-2100.3.34.2) | libc++ | arm64 | native | `7874afeb97cb5ba7d7a5b1c5771c15ff2d40228ba92d0c52efe8dcd87431ac67` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |
| libstdcxx | gcc 15.2.0 | libstdc++ | arm64 | native | `c999a548375281e6afa30d545d4d1acaaba8a77de5a571b76493780e81c2fc01` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |

## Adjudication

- 4 audit stream(s) adjudicated by `boost::multiprecision::cpp_dec_float_100`:
  - `raw/p2/audit_p2_libcxx.bin`
  - `raw/p2/audit_p2_libstdcxx.bin`
  - `raw/p3/audit_p3_libcxx.bin`
  - `raw/p3/audit_p3_libstdcxx.bin`

## Analysis inputs

| file | SHA-256 |
|---|---|
| `PROTOCOL.md` | `45992ff0ce8f877b58ec7e5e6af4574a45f964c6b53117864dc30bb42756dd7c` |
| `config/protocol.json` | `c3b4a0a07877c0f8e2370ff7fdae98d3f3d96c15465aadb9eae00a280391c6a6` |
| `config/power_study.json` | `bac1196402a4dbe345d9f35a8e2e445979e0896f7957860764942924f60a080b` |
| `config/honest_floor.md` | `a86bbe17b579e92f5596bbc7a1110b23bd9dfabf567de56c9792668404dee8cd` |
| `results/schema.md` | `e2543ad027cd9234596fe5c2820c038f7ca7e56ebb953542c0a6c7c70512a41d` |
| `analyze.py` | `a50aa97d0d9db4302d4847194eed716150b38a49a418235229fe44f8dfede101` |
| `exp7_families.py` | `c8950b3e2bbc08d7958ab0529b2f5745d3861a88a3b4b73aa7c1bc4c19d113a9` |
| `exp7_gates.py` | `be9dba2be45e621c087696c77b4c8a8aa3b91ecfba6ce6091d35e145e813b8a1` |
| `exp7_stats.py` | `df4d0df9dfd8ce72bc3971bdc37e0b2b8afca596654a7db9c5f726a7c976e21a` |
| `exp7_io.py` | `4fed3570ba6ebae8a92c0f412aeaded2fb5d14ab958e501eb82023f007ebad15` |
| `exp7_portability.py` | `3df7ff069a37f435014bb60fb0086b982418a880c5782e4858d6df80d62d06fc` |
