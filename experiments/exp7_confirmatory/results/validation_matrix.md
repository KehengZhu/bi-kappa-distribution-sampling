# Experiment 7 - complete-loader validation matrix

Protocol `config/protocol.json` SHA-256 `7d6e05fbbf3cd9f864b265ffa6122b0ea4fd88757752dbc00ba634762d91b76b`.

Rows are the curated configurations of PROTOCOL.md 4.1, pooled over the frozen seeds; columns are the tests of family F5 -- the five frozen with the protocol, then the upper-tail members amendment 1.3.0 added at each threshold in `F5_tail.q0`. Decisions are Holm-corrected jointly over **all** cells and tests at a familywise alpha of 0.01, not within each cell.

`n.a.` marks a test that is not a property of the cell's law rather than one that passed: under a cap the accepted set couples the radius to the direction and is not rotationally symmetric in the azimuth, so direction uniformity, independence and frame invariance do not apply there and the cap law is tested by its own conditional transform instead.

**A conditional cell cannot support the fidelity claim on its own.** A cell is conditional when the statistics ran on a subsample that something correlated with the radius had already filtered -- a non-finite component is exactly what a large radius produces. Every such cell carries its loss fraction in the last column and in `loader_validation.csv`.

A tail count is the number of draws above `z0 = -log q0`, the unresolved draws included, against `Binomial(attempts, q0)`; the excess test is a Kolmogorov-Smirnov of the resolved excesses against `Exp(1)`. Because an unresolved draw counts toward every threshold, the count member also reads on any cell whose loss fraction is itself above `q0`, and `loader_validation.csv` publishes the resolved and unresolved parts of every count separately.

| case | method | tag | direction uniformity | independence | frame invariance | anisotropy | cap law | tail 1e-2 count | tail 1e-2 excess | tail 1e-3 count | tail 1e-3 excess | tail 1e-4 count | tail 1e-4 excess | n analysed | conditional |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C0 | CANDIDATE | libcxx | pass (p=0.844) | pass (p=0.673) | pass (p=0.113) | pass (p=0.309) | n.a. | pass (p=0.121) | pass (p=0.951) | pass (p=0.283) | pass (p=0.92) | pass (p=0.887) | pass (p=0.917) | 500000 | unconditional |
| C0 | CANDIDATE | libstdcxx | p=0.844 (not in F5) | p=0.673 (not in F5) | p=0.113 (not in F5) | p=0.309 (not in F5) | n.a. | p=0.121 (not in F5) | p=0.951 (not in F5) | p=0.283 (not in F5) | p=0.92 (not in F5) | p=0.887 (not in F5) | p=0.917 (not in F5) | 500000 | unconditional |
| C0 | LEGACY | libcxx | p=0.853 (not in F5) | p=0.66 (not in F5) | p=0.209 (not in F5) | p=0.0421 (not in F5) | n.a. | p=0.527 (not in F5) | p=0.493 (not in F5) | p=0.134 (not in F5) | p=0.648 (not in F5) | p=0.0232 (not in F5) | p=0.916 (not in F5) | 500000 | unconditional |
| C0 | LEGACY | libstdcxx | p=0.353 (not in F5) | p=0.644 (not in F5) | p=0.911 (not in F5) | p=0.469 (not in F5) | n.a. | p=0.274 (not in F5) | p=0.993 (not in F5) | p=0.516 (not in F5) | p=0.538 (not in F5) | p=0.288 (not in F5) | p=0.115 (not in F5) | 500000 | unconditional |
| C1 | CANDIDATE | libcxx | pass (p=0.802) | pass (p=0.662) | pass (p=0.296) | pass (p=0.924) | n.a. | pass (p=0.744) | pass (p=0.882) | pass (p=0.858) | pass (p=0.234) | pass (p=0.571) | pass (p=0.348) | 499999 | conditional, loss fraction 2e-06 |
| C1 | CANDIDATE | libstdcxx | p=0.802 (not in F5) | p=0.662 (not in F5) | p=0.296 (not in F5) | p=0.924 (not in F5) | n.a. | p=0.744 (not in F5) | p=0.882 (not in F5) | p=0.858 (not in F5) | p=0.234 (not in F5) | p=0.571 (not in F5) | p=0.348 (not in F5) | 499999 | conditional, loss fraction 2e-06 |
| C1 | LEGACY | libcxx | p=0.947 (not in F5) | p=0.271 (not in F5) | p=0.746 (not in F5) | p=0.402 (not in F5) | n.a. | p=0.654 (not in F5) | p=1.94e-14 (not in F5) | p=0.46 (not in F5) | p=2.6e-74 (not in F5) | p=1.56e-91 (not in F5) | n.a. | 499748 | conditional, loss fraction 0.000504 |
| C1 | LEGACY | libstdcxx | p=0.769 (not in F5) | p=0.901 (not in F5) | p=0.255 (not in F5) | p=0.576 (not in F5) | n.a. | p=0.644 (not in F5) | p=7.29e-15 (not in F5) | p=0.303 (not in F5) | p=4.68e-65 (not in F5) | p=1.48e-149 (not in F5) | n.a. | 499672 | conditional, loss fraction 0.000656 |
| C2 | CANDIDATE | libcxx | pass (p=0.802) | pass (p=0.662) | pass (p=0.296) | pass (p=0.924) | n.a. | pass (p=0.744) | pass (p=0.882) | pass (p=0.858) | pass (p=0.234) | pass (p=0.571) | pass (p=0.348) | 499999 | conditional, loss fraction 2e-06 |
| C2 | CANDIDATE | libstdcxx | p=0.802 (not in F5) | p=0.662 (not in F5) | p=0.296 (not in F5) | p=0.924 (not in F5) | n.a. | p=0.744 (not in F5) | p=0.882 (not in F5) | p=0.858 (not in F5) | p=0.234 (not in F5) | p=0.571 (not in F5) | p=0.348 (not in F5) | 499999 | conditional, loss fraction 2e-06 |
| C2 | LEGACY | libcxx | p=0.947 (not in F5) | p=0.271 (not in F5) | p=0.746 (not in F5) | p=0.402 (not in F5) | n.a. | p=0.654 (not in F5) | p=1.94e-14 (not in F5) | p=0.46 (not in F5) | p=2.6e-74 (not in F5) | p=1.56e-91 (not in F5) | n.a. | 499748 | conditional, loss fraction 0.000504 |
| C2 | LEGACY | libstdcxx | p=0.769 (not in F5) | p=0.901 (not in F5) | p=0.255 (not in F5) | p=0.576 (not in F5) | n.a. | p=0.644 (not in F5) | p=7.29e-15 (not in F5) | p=0.303 (not in F5) | p=4.68e-65 (not in F5) | p=1.48e-149 (not in F5) | n.a. | 499672 | conditional, loss fraction 0.000656 |
| C3 | CANDIDATE | libcxx | pass (p=0.945) | pass (p=0.888) | pass (p=0.399) | pass (p=0.415) | n.a. | pass (p=0.909) | pass (p=0.167) | pass (p=0.14) | **FAIL** (p=2.92e-09) | n.a. | n.a. | 499912 | conditional, loss fraction 0.000176 |
| C3 | CANDIDATE | libstdcxx | p=0.945 (not in F5) | p=0.888 (not in F5) | p=0.399 (not in F5) | p=0.415 (not in F5) | n.a. | p=0.909 (not in F5) | p=0.167 (not in F5) | p=0.14 (not in F5) | **FAIL** (p=2.92e-09) | n.a. | n.a. | 499912 | conditional, loss fraction 0.000176 |
| C3 | LEGACY | libcxx | p=0.203 (not in F5) | p=0.596 (not in F5) | p=0.251 (not in F5) | p=0.529 (not in F5) | n.a. | p=0.316 (not in F5) | p=0 (not in F5) | p=0 (not in F5) | n.a. | n.a. | n.a. | 497213 | conditional, loss fraction 0.00557 |
| C3 | LEGACY | libstdcxx | p=0.474 (not in F5) | p=0.516 (not in F5) | p=0.422 (not in F5) | p=0.116 (not in F5) | n.a. | p=0.157 (not in F5) | p=0 (not in F5) | p=0 (not in F5) | n.a. | n.a. | n.a. | 497206 | conditional, loss fraction 0.00559 |
| C4 | CANDIDATE | libcxx | pass (p=0.933) | pass (p=0.585) | pass (p=0.7) | pass (p=0.869) | n.a. | pass (p=0.955) | **FAIL** (p=3.35e-28) | pass (p=1) | **FAIL** (p=3.09e-63) | n.a. | n.a. | 499583 | conditional, loss fraction 0.000834 |
| C4 | CANDIDATE | libstdcxx | p=0.933 (not in F5) | p=0.585 (not in F5) | p=0.7 (not in F5) | p=0.869 (not in F5) | n.a. | p=0.955 (not in F5) | **FAIL** (p=3.35e-28) | p=1 (not in F5) | **FAIL** (p=3.09e-63) | n.a. | n.a. | 499583 | conditional, loss fraction 0.000834 |
| C4 | LEGACY | libcxx | p=0.986 (not in F5) | p=0.202 (not in F5) | p=0.83 (not in F5) | p=0.373 (not in F5) | n.a. | p=0 (not in F5) | n.a. | p=0 (not in F5) | n.a. | n.a. | n.a. | 487835 | conditional, loss fraction 0.0243 |
| C4 | LEGACY | libstdcxx | p=0.521 (not in F5) | p=0.949 (not in F5) | p=0.573 (not in F5) | p=0.812 (not in F5) | n.a. | p=0 (not in F5) | n.a. | p=0 (not in F5) | n.a. | n.a. | n.a. | 487929 | conditional, loss fraction 0.0241 |
| C5 | CANDIDATE | libcxx | n.a. | n.a. | n.a. | pass (p=0.872) | pass (p=0.908) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C5 | CANDIDATE | libstdcxx | n.a. | n.a. | n.a. | p=0.872 (not in F5) | p=0.908 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C5 | LEGACY | libcxx | n.a. | n.a. | n.a. | p=0.765 (not in F5) | p=0.164 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C5 | LEGACY | libstdcxx | n.a. | n.a. | n.a. | p=0.418 (not in F5) | p=0.66 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C6 | CANDIDATE | libcxx | n.a. | n.a. | n.a. | pass (p=0.819) | pass (p=0.968) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C6 | CANDIDATE | libstdcxx | n.a. | n.a. | n.a. | p=0.819 (not in F5) | p=0.968 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C6 | LEGACY | libcxx | n.a. | n.a. | n.a. | p=0.66 (not in F5) | p=0.592 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |
| C6 | LEGACY | libstdcxx | n.a. | n.a. | n.a. | p=0.0237 (not in F5) | p=0.599 (not in F5) | n.a. | n.a. | n.a. | n.a. | n.a. | n.a. | 500000 | unconditional |

Only the CANDIDATE cells of the primary environment (`libcxx`) enter the F5 decision: within version 2.0.0 the stream is a function of the engine alone, so the other environment's rows are the same draws and would enter Holm twice. The LEGACY rows are the comparator and are printed with their p-values but are not gated -- where the released 1.0.0 form discards a non-negligible share of the intended draws, the draws it keeps are no longer distributed as the target, and a test detecting that is the conditioning result rather than a defect in this analysis.

F5 global: FAILS. Holm jointly over 50 tests across all cells; 36 cells are conditional on success and are labelled as such

