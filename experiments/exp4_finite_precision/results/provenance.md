# Experiment 7 - provenance of this analysis

This is the only file under results/ that carries wall-clock times and build identity. Every value below is copied from raw/manifest.csv and the environment record written at run time; nothing is read from the clock or the host while the analysis runs, so this file regenerates byte for byte and `make reverify` remains a meaningful check.

- protocol: `config/protocol.json` version 5.0.0, SHA-256 `8fa713c11bebc4a76f9d1a412633a7fdb72f7be1ddce397112801c1b400f3596`
- mode: production
- manifest: `raw/manifest.csv`, 846 distinct output files recorded in 846 manifest rows

- probe run times (UTC, from the writing processes): 2026-09-25T14:30:50Z .. 2026-09-25T15:01:41Z

## Builds that produced the data

| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |
|---|---|---|---|---|---|---|
| libcxx | clang 21.0.0 (clang-2100.3.34.2) | libc++ | arm64 | native | `d737ba2ebff2de789a6eee738c7f82ff8adbccef9ecbebdd68b41dbdc540003c` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |
| libstdcxx | gcc 15.2.0 | libstdc++ | arm64 | native | `daa0ffec476d4a076f68b0099908b47f90093065dfcf60ea999fa32a5c5cfad3` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |

## Adjudication

- 4 audit stream(s) adjudicated by `boost::multiprecision::mpfr_float_100`:
  - `raw/p2/audit_p2_libcxx.bin`
  - `raw/p2/audit_p2_libstdcxx.bin`
  - `raw/p3/audit_p3_libcxx.bin`
  - `raw/p3/audit_p3_libstdcxx.bin`

## Analysis inputs

| file | SHA-256 |
|---|---|
| `PROTOCOL.md` | `fdaf227dec99e069f5a166388ea5836b83a258dbb9117831d721578a77901dc3` |
| `config/protocol.json` | `8fa713c11bebc4a76f9d1a412633a7fdb72f7be1ddce397112801c1b400f3596` |
| `config/power_study.json` | `bac1196402a4dbe345d9f35a8e2e445979e0896f7957860764942924f60a080b` |
| `config/honest_floor.md` | `a86bbe17b579e92f5596bbc7a1110b23bd9dfabf567de56c9792668404dee8cd` |
| `results/schema.md` | `a4f9e12ab539296dd442b7cd821c6ccb6f83500e1ed580336a4d557f15c9d25f` |
| `analyze.py` | `dce9182090d91e00cc188ed139c53396aa7010634731d1d28877fce9ca045e14` |
| `exp7_families.py` | `c8950b3e2bbc08d7958ab0529b2f5745d3861a88a3b4b73aa7c1bc4c19d113a9` |
| `exp7_gates.py` | `f10bd5bc8f088a88282eed3f570c9e5074e6bb7693a96c3ebb4926d4b9a4b7e0` |
| `exp7_stats.py` | `df4d0df9dfd8ce72bc3971bdc37e0b2b8afca596654a7db9c5f726a7c976e21a` |
| `exp7_io.py` | `4fed3570ba6ebae8a92c0f412aeaded2fb5d14ab958e501eb82023f007ebad15` |
| `exp7_portability.py` | `3df7ff069a37f435014bb60fb0086b982418a880c5782e4858d6df80d62d06fc` |
