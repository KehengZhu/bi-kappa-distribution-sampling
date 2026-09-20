# Experiment 7 - provenance of this analysis

This is the only file under results/ that carries wall-clock times and build identity. Every value below is copied from raw/manifest.csv and the environment record written at run time; nothing is read from the clock or the host while the analysis runs, so this file regenerates byte for byte and `make reverify` remains a meaningful check.

- protocol: `config/protocol.json` version 1.3.0, SHA-256 `7d6e05fbbf3cd9f864b265ffa6122b0ea4fd88757752dbc00ba634762d91b76b`
- mode: production
- manifest: `raw/manifest.csv`, 846 distinct output files recorded in 846 manifest rows

- probe run times (UTC, from the writing processes): 2026-09-20T04:27:45Z .. 2026-09-20T04:38:24Z

## Builds that produced the data

| tag | compiler | stdlib | arch | execution | probe SHA-256 | cxxflags |
|---|---|---|---|---|---|---|
| libcxx | clang 21.0.0 (clang-2100.3.34.2) | libc++ | arm64 | native | `c100ec2bd199e42a0d72c99320d3818dc291fb254537ffb976dbac8dbd4e257a` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |
| libstdcxx | gcc 15.2.0 | libstdc++ | arm64 | native | `ae547cf3e084cebfc7c7188d56d9db247a61be16e59bebcdf5168bc86375464e` | `-Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src` |

## Adjudication

- 4 audit stream(s) adjudicated by `boost::multiprecision::cpp_dec_float_100`:
  - `raw/p2/audit_p2_libcxx.bin`
  - `raw/p2/audit_p2_libstdcxx.bin`
  - `raw/p3/audit_p3_libcxx.bin`
  - `raw/p3/audit_p3_libstdcxx.bin`

## Analysis inputs

| file | SHA-256 |
|---|---|
| `PROTOCOL.md` | `a476c12898f1697744780e94d6dcf5ae62b44af14d7fd7ffa7e3ca8547b55aa1` |
| `config/protocol.json` | `7d6e05fbbf3cd9f864b265ffa6122b0ea4fd88757752dbc00ba634762d91b76b` |
| `config/power_study.json` | `bac1196402a4dbe345d9f35a8e2e445979e0896f7957860764942924f60a080b` |
| `config/honest_floor.md` | `a86bbe17b579e92f5596bbc7a1110b23bd9dfabf567de56c9792668404dee8cd` |
| `results/schema.md` | `71a683e382bc89648a30a2a6ea7434b9f8d38ec7cb728992e9c5fcb6f205ce09` |
| `analyze.py` | `717c28e0f280469a158dcbdba2398ff563c14c06be2acc1ad1ec5bf1953e7813` |
| `exp7_families.py` | `b04996449bd4ed6057b27adf129f8ae8e51ae9087ceaf47f83c7d5ca86e491ce` |
| `exp7_gates.py` | `de2e5591560c8d521ae14e2b85bdbb25ef66cb6d497ff728bac83ca450d00dac` |
| `exp7_stats.py` | `df4d0df9dfd8ce72bc3971bdc37e0b2b8afca596654a7db9c5f726a7c976e21a` |
| `exp7_io.py` | `4fed3570ba6ebae8a92c0f412aeaded2fb5d14ab958e501eb82023f007ebad15` |
| `exp7_portability.py` | `3df7ff069a37f435014bb60fb0086b982418a880c5782e4858d6df80d62d06fc` |
