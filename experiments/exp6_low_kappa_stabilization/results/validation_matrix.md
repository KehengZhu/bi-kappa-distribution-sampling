# E4 - complete-loader validation matrix

Generated 2026-09-19T17:46:07+00:00.

Decisions are Holm-corrected within each configuration at a familywise alpha of 0.01.  `pass` means the frozen test did not reject; it is consistency at the tested sample size, not proof of exactness.  NC1 and NC2 pass only when the injected defect *is* detected.

| tag | case | method | n | n finite | radial | direction | independence | anisotropy | frame | finiteness | cap_law |
|---|---|---|---:|---:|---|---|---|---|---|---|---|
| libcxx | C0 | LOG | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C0 | SPLIT | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C0 | SPLIT_released | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C1 | LOG | 500000 | 499998 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C1 | SPLIT | 500000 | 499722 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C1 | SPLIT_released | 500000 | 499722 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C2 | LOG | 500000 | 499998 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C2 | SPLIT | 500000 | 499722 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C2 | SPLIT_released | 500000 | 499722 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C3 | LOG | 500000 | 499922 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C3 | SPLIT | 500000 | 497137 | fail | pass | pass | pass | pass | pass | na |
| libcxx | C3 | SPLIT_released | 500000 | 497137 | fail | pass | pass | pass | pass | pass | na |
| libcxx | C4 | LOG | 500000 | 499622 | pass | pass | pass | pass | pass | pass | na |
| libcxx | C4 | SPLIT | 500000 | 487946 | fail | pass | pass | pass | pass | pass | na |
| libcxx | C4 | SPLIT_released | 500000 | 487946 | fail | pass | pass | pass | pass | pass | na |
| libcxx | C5 | LOG | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libcxx | C5 | SPLIT | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libcxx | C6 | LOG | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libcxx | C6 | SPLIT | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libcxx | NC1 | LOG | 499998 | 499998 | na | na | pass | na | na | na | na |
| libcxx | NC2 | LOG | 500000 | 500000 | pass | na | na | na | na | na | na |
| libstdcxx | C0 | LOG | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C0 | SPLIT | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C0 | SPLIT_released | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C1 | LOG | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C1 | SPLIT | 500000 | 499702 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C1 | SPLIT_released | 500000 | 499702 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C2 | LOG | 500000 | 500000 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C2 | SPLIT | 500000 | 499702 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C2 | SPLIT_released | 500000 | 499702 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C3 | LOG | 500000 | 499926 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C3 | SPLIT | 500000 | 497105 | fail | pass | pass | pass | pass | pass | na |
| libstdcxx | C3 | SPLIT_released | 500000 | 497105 | fail | pass | pass | pass | pass | pass | na |
| libstdcxx | C4 | LOG | 500000 | 499611 | pass | pass | pass | pass | pass | pass | na |
| libstdcxx | C4 | SPLIT | 500000 | 487936 | fail | pass | pass | pass | pass | pass | na |
| libstdcxx | C4 | SPLIT_released | 500000 | 487936 | fail | pass | pass | pass | pass | pass | na |
| libstdcxx | C5 | LOG | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libstdcxx | C5 | SPLIT | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libstdcxx | C6 | LOG | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libstdcxx | C6 | SPLIT | 500000 | 500000 | na | na | na | pass | pass | pass | pass |
| libstdcxx | NC1 | LOG | 500000 | 500000 | na | na | pass | na | na | na | na |
| libstdcxx | NC2 | LOG | 500000 | 500000 | pass | na | na | na | na | na | na |

## Released-path radial rejections on the finite survivors

These are results, not defects. Where the released formation discards a non-negligible share of the intended draws, the draws it keeps are no longer distributed as the target: the radial test rejects on the survivor subset. The log path passes the same test on the same configuration. The size of the effect is in `tail_metrics.csv`.

- libcxx / C3 / SPLIT: 2863 of 500000 draws discarded (5.73e-03 per attempt); radial law rejected on the remaining survivors
- libcxx / C3 / SPLIT_released: 2863 of 500000 draws discarded (5.73e-03 per attempt); radial law rejected on the remaining survivors
- libcxx / C4 / SPLIT: 12054 of 500000 draws discarded (2.41e-02 per attempt); radial law rejected on the remaining survivors
- libcxx / C4 / SPLIT_released: 12054 of 500000 draws discarded (2.41e-02 per attempt); radial law rejected on the remaining survivors
- libstdcxx / C3 / SPLIT: 2895 of 500000 draws discarded (5.79e-03 per attempt); radial law rejected on the remaining survivors
- libstdcxx / C3 / SPLIT_released: 2895 of 500000 draws discarded (5.79e-03 per attempt); radial law rejected on the remaining survivors
- libstdcxx / C4 / SPLIT: 12064 of 500000 draws discarded (2.41e-02 per attempt); radial law rejected on the remaining survivors
- libstdcxx / C4 / SPLIT_released: 12064 of 500000 draws discarded (2.41e-02 per attempt); radial law rejected on the remaining survivors
