# JTJ1001 major-revision working documents

Hand-written planning material for the major revision of *Sampling the Bi-Kappa Distribution*
(APS Open Science, JTJ1001). **Nothing here is generated** — `doxygen Doxyfile` writes only into
`docs/api/` and never touches this directory.

These are internal working documents, not manuscript text and not a rebuttal draft.

## Layout

| Path | Contents |
|---|---|
| `planning/reviewer_response_matrix.md` | The authoritative record. Reviewer-by-reviewer analysis (R1.1–R1.7, R2.A1–R2.C3), the experiment programme, the contribution ledger, and the **alignment protocol + log** (§9). |
| `planning/introduction_rewrite_proposal.md` | Draft Introduction and contribution statement, plus the low-κ risk audit that preceded Experiment 4. |
| `literature/step1_claim_audit.md` | Every planned literature claim against its primary source, with a PASS / MODIFY / DROP / BLOCKED verdict. §7 is the second round. |
| `experiments/` | Cross-experiment notes. The experiments themselves live in `experiments/` at the repository root, each with its own `README.md` and `results/`. |

## Where the evidence actually lives

| | Location | Status |
|---|---|---|
| Experiment 1 — radial / directional / anisotropic / frame validation | `experiments/exp1_radial_directional/` | **complete** |
| Experiment 2 — capped vs uncapped characterization | `experiments/exp2_cap_characterization/` | **complete** |
| Experiment 3 — performance benchmark | — | not started; prerequisite in matrix §9.5 |
| Experiment 4 — finite-precision / low-κ audit | `experiments/exp4_precision/` | **complete** |
| Primary sources read | `paper/reference/README.md` | 20 sources |

## Reading order

1. `planning/reviewer_response_matrix.md` §1 — revision context and what the referees asked.
2. `literature/step1_claim_audit.md` — what the literature actually supports. **Read before
   writing any prose**; it is the document that most often reverses an intended claim.
3. `planning/reviewer_response_matrix.md` §9 — current alignment state, contribution ledger,
   and open blockers.
4. `planning/introduction_rewrite_proposal.md` — only after the above.

## Keeping these in sync

Follow the alignment protocol in `planning/reviewer_response_matrix.md` §9.1, and add a row to
the §9.2 log. The two standing rules:

- **Superseded recommendations are marked withdrawn in place, never deleted silently** — a
  returning referee who read the earlier reasoning would otherwise see an unexplained reversal.
  Live example: `Y = T/(1+T)` → `W = 1/(1+T)`.
- **Every quantitative statement travels with its parameter range** or it does not travel.
