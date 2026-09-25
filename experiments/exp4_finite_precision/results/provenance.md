# Experiment 7 - provenance of this analysis

This is the only file under results/ that carries wall-clock times and build identity. Every value below is copied from raw/manifest.csv and the environment record written at run time; nothing is read from the clock or the host while the analysis runs, so this file regenerates byte for byte and `make reverify` remains a meaningful check.

- protocol: `config/protocol.json` version 4.0.0, SHA-256 `664d906899b1726ae6d0bd3856483c4590af875a9dfdc3fee14fdcb79dbc4e34`
- mode: production
- manifest: `raw/manifest.csv`, 846 distinct output files recorded in 846 manifest rows

- probe run times (UTC, from the writing processes): 2026-09-25T05:25:13Z .. 2026-09-25T05:55:38Z

## Builds that produced the data

| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |
|---|---|---|---|---|---|---|
| libcxx | clang 21.0.0 (clang-2100.3.34.2) | libc++ | arm64 | native | `2f81b486545542fc7b2159aaeec4d7169505046154a97bb1469c18bc24d28fac` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |
| libstdcxx | gcc 15.2.0 | libstdc++ | arm64 | native | `62a584d372abbb1724b224dc18b88c277d6ee16b754484eecb896f7fbdabeca5` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |

## Adjudication

- 4 audit stream(s) adjudicated by `boost::multiprecision::mpfr_float_100`:
  - `raw/p2/audit_p2_libcxx.bin`
  - `raw/p2/audit_p2_libstdcxx.bin`
  - `raw/p3/audit_p3_libcxx.bin`
  - `raw/p3/audit_p3_libstdcxx.bin`

## Analysis inputs

| file | SHA-256 |
|---|---|
| `PROTOCOL.md` | `018ce02e01e33939ae6b14eed1f7844a392bb3602183ddcf1f7ae852a19b9768` |
| `config/protocol.json` | `664d906899b1726ae6d0bd3856483c4590af875a9dfdc3fee14fdcb79dbc4e34` |
| `config/power_study.json` | `bac1196402a4dbe345d9f35a8e2e445979e0896f7957860764942924f60a080b` |
| `config/honest_floor.md` | `a86bbe17b579e92f5596bbc7a1110b23bd9dfabf567de56c9792668404dee8cd` |
| `results/schema.md` | `0617e2703b8f5e0ebe390b3ebfce17d9fc7f70ea9636a27342f000c47b8a6c8a` |
| `analyze.py` | `1faac741f55a6160f293b92eab32d6d71cfe6af76c83fc2923715092d3455fa7` |
| `exp7_families.py` | `c8950b3e2bbc08d7958ab0529b2f5745d3861a88a3b4b73aa7c1bc4c19d113a9` |
| `exp7_gates.py` | `7a7ab5f51aed01008d451816ba6c40cc348e275f63b73ea3fb99cbbe7530426d` |
| `exp7_stats.py` | `df4d0df9dfd8ce72bc3971bdc37e0b2b8afca596654a7db9c5f726a7c976e21a` |
| `exp7_io.py` | `4fed3570ba6ebae8a92c0f412aeaded2fb5d14ab958e501eb82023f007ebad15` |
| `exp7_portability.py` | `3df7ff069a37f435014bb60fb0086b982418a880c5782e4858d6df80d62d06fc` |
