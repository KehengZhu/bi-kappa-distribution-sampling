# Experiment 7 - provenance of this analysis

This is the only file under results/ that carries wall-clock times and build identity. Every value below is copied from raw/manifest.csv and the environment record written at run time; nothing is read from the clock or the host while the analysis runs, so this file regenerates byte for byte and `make reverify` remains a meaningful check.

- protocol: `config/protocol.json` version 2.0.0, SHA-256 `914c5d1d2b6900d872402a545f46a99819038987416d1b882909eb4ca2bafa48`
- mode: production
- manifest: `raw/manifest.csv`, 846 distinct output files recorded in 846 manifest rows

- probe run times (UTC, from the writing processes): 2026-09-20T16:08:30Z .. 2026-09-20T16:19:39Z

## Builds that produced the data

| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |
|---|---|---|---|---|---|---|
| libcxx | clang 21.0.0 (clang-2100.3.34.2) | libc++ | arm64 | native | `e481c63e0c9236e265479a1a5a3a13e29c2751ac5322aff3184fb12359063b1b` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |
| libstdcxx | gcc 15.2.0 | libstdc++ | arm64 | native | `a84df0e922eda88008fd78f22358f105c3696f0471ac0919c21eadd3f035d27f` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |

## Adjudication

- 4 audit stream(s) adjudicated by `boost::multiprecision::cpp_dec_float_100`:
  - `raw/p2/audit_p2_libcxx.bin`
  - `raw/p2/audit_p2_libstdcxx.bin`
  - `raw/p3/audit_p3_libcxx.bin`
  - `raw/p3/audit_p3_libstdcxx.bin`

## Analysis inputs

| file | SHA-256 |
|---|---|
| `PROTOCOL.md` | `021c0b91ae4ee912c26259d11f23e82e5c81d76fade12e75ca6013c89de34857` |
| `config/protocol.json` | `914c5d1d2b6900d872402a545f46a99819038987416d1b882909eb4ca2bafa48` |
| `config/power_study.json` | `bac1196402a4dbe345d9f35a8e2e445979e0896f7957860764942924f60a080b` |
| `config/honest_floor.md` | `a86bbe17b579e92f5596bbc7a1110b23bd9dfabf567de56c9792668404dee8cd` |
| `results/schema.md` | `71a683e382bc89648a30a2a6ea7434b9f8d38ec7cb728992e9c5fcb6f205ce09` |
| `analyze.py` | `7388b7c12100e437310e03b1d2cb9f4e7320d5160341fb72f11281b973fd1d0d` |
| `exp7_families.py` | `b04996449bd4ed6057b27adf129f8ae8e51ae9087ceaf47f83c7d5ca86e491ce` |
| `exp7_gates.py` | `c149b2cdb0b626ff8ffca873d7e070d35b87795e694a5a7b8dcb15c1f224645a` |
| `exp7_stats.py` | `df4d0df9dfd8ce72bc3971bdc37e0b2b8afca596654a7db9c5f726a7c976e21a` |
| `exp7_io.py` | `4fed3570ba6ebae8a92c0f412aeaded2fb5d14ab958e501eb82023f007ebad15` |
| `exp7_portability.py` | `3df7ff069a37f435014bb60fb0086b982418a880c5782e4858d6df80d62d06fc` |
