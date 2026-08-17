# bi-kappa

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20617011.svg)](https://doi.org/10.5281/zenodo.20617011)

Monte Carlo particle samplers for plasma **velocity** and **position** distributions, in C++ and Python.

- **C++** (`cpp/`) — header-only samplers: bi-kappa & bi-Maxwellian velocities, plus rejection samplers for any speed-squared (|v|²) or spatial density you define.
- **Python** (`python/`) — equivalent general samplers + Jupyter notebooks to visualize the output.

📖 **API reference:** [kehengzhu.github.io/bi-kappa-distribution-sampling](https://kehengzhu.github.io/bi-kappa-distribution-sampling/). Locally, open [`docs/index.html`](docs/index.html) in a browser.

---

## Quick start — C++

```bash
cd cpp && make && ./main.exe     # runs tests, then writes samples_*.txt
```

To use a sampler in your own code, every generator follows the same three steps — **construct → `define(...)` → call**:

```cpp
#include "bi_kappa_distribution.H"

bi_kappa_distribution<double> dist;
dist.define(/*kappa*/ 2.0, /*theta_perp*/ 1.0, /*theta_par*/ 2.0,
            /*ub*/ {0, 0, 1}, /*cap*/ 20.0, /*seed*/ 12345);

auto v = dist();     // v = {vx, vy, vz}, sampled in the global frame
```

*Prerequisite: any C++11 compiler. Set `CXX` in `cpp/GNUmakefile` if not `g++`.*

## The five samplers

| Header | Samples | define(...) takes |
|---|---|---|
| `bi_kappa_distribution` | 3D bi-kappa velocity | `kappa, theta_perp, theta_par, ub, cap[, seed]` |
| `bi_maxwellian_distribution` | 3D bi-Maxwellian velocity | `theta_perp, theta_par, ub, cap` |
| `general_velocity_generator` | speed from your `g(w)`, `w = \|v\|²`, isotropic direction | `g, v²_min, v²_max` |
| `field_aligned_velocity_generator` | parallel speed from your `g(w)`, `w = v_par²`, along `ub`; Maxwellian perpendicular | `g, v²_min, v²_max, theta_perp, ub, sign` |
| `general_position_generator` | position from your density `rho(x)` | `dimension, lower, upper, rho` |

They all share the same behavior:

- **Construct → `define(...)` → call.** Calling before `define(...)` throws.
- **Seeding:** call `dist.seed(s)` any time. `bi_kappa_distribution` also accepts the seed as the last `define(...)` argument. A **negative** `s` draws an unpredictable seed from `std::random_device`.
- **No-arg call** `dist()` uses the sampler's own RNG. **Bring-your-own** `dist(gen)` accepts an external `std::mt19937`.
- **`ub`** is the magnetic-field direction (need not be unit length). Output is rotated into the global frame; the default `{0,0,1}` returns field-aligned components directly.
- **`cap`** selects which of **two distinct target distributions** you sample.
  - A **finite** `cap` (default `20`) rejects any component with `|v_i| / theta_i > cap` and
    resamples, throwing after 10⁶ failed tries. The resulting samples follow the bi-Kappa
    distribution **conditioned on** all normalized components lying inside that box — a different,
    bounded law with its own normalization, not the bi-Kappa distribution itself. Its moments are
    finite for every `kappa`, because it is a different distribution, not because truncation
    regularizes the original one.
  - **`bi_kappa_distribution<double>::no_cap()`** disables the cap. Samples then follow the full
    bi-Kappa law, and the draw uses a fixed sequence of high-level variates with no outer
    acceptance–rejection loop. Second moments diverge for `kappa <= 3/2`; that is a property of
    the distribution, not a defect.
  - `dist.param().capped()` reports which mode is active.
- **Change one parameter** without re-specifying the rest: `dist.kappa(3.0)`, `dist.ub({0,1,0})`, etc.

For the full API — every overload, parameter constraint, and the Python classes — see the [API reference](https://kehengzhu.github.io/bi-kappa-distribution-sampling/). A worked example using all five samplers lives in [`cpp/main.cpp`](cpp/main.cpp).

## Quick start — Python

```bash
cd python
uv sync                # or: pip install numpy scipy pandas matplotlib ipykernel
python general_generators.py     # writes samples_general_{velocity,position}.txt and samples_field_aligned.txt
```

`general_generators.py` provides `GeneralVelocityGenerator`, `FieldAlignedVelocityGenerator`, and `GeneralPositionGenerator` — the Python equivalents of the C++ general samplers.

## Visualize

Open [`python/visualize_samples.ipynb`](python/visualize_samples.ipynb) and run all cells. Point `workspace_root` at your sample files:

- `"./"` — files in `python/`
- `"../cpp/"` — files written by the C++ driver

## Tests

`./main.exe` runs `run_all_tests()` before sampling and **exits non-zero on any failure**. Coverage: frame-transform correctness, `define(...)`/getter behavior, the velocity cap, seed reproducibility, the radial formation (`test_radius_formation`), and uncapped/capped target-law semantics (`test_no_cap_semantics`).

## Validation experiments

`experiments/` holds the standalone validation studies behind the manuscript. Each is
self-contained and reproducible from one documented command, with committed scripts, fixed
seeds, a manifest, and an environment record.

| Directory | What it establishes |
|---|---|
| `exp1_radial_directional/` | The central radial law `T = R² ~ β'(3/2, κ−1/2)`, directional uniformity, radial–direction independence, anisotropy, and arbitrary **B**-frame invariance. 1.35×10⁷ draws, κ = 0.51…10. |
| `exp2_cap_characterization/` | How the optional component-wise cap changes the sampled law: it is the bi-Kappa **conditioned on the box**, not the same distribution. The rejected fraction equals the total-variation distortion, decays only as `λ^−(2κ−1)`, and the box breaks axisymmetry about **B**. |
| `exp4_precision/` | Supported numerical range as a function of κ, precision and standard library. `double` is clean for κ ≥ 0.55, `float` for κ ≥ 0.75. |

Bulk per-draw output (`raw/*.bin`) is regenerable and deliberately not tracked; manifests,
compact summaries and manuscript-quality figures are.

## Documentation

The full API reference is generated from the in-source comments with [Doxygen](https://www.doxygen.nl/) and lives in [`docs/api/`](docs/api/).

**Read it locally** — open [`docs/api/index.html`](docs/api/index.html) in a browser. (Note: clicking that file on github.com shows the HTML *source*, not the rendered page — use GitHub Pages below for a live site.)

**Regenerate after changing code or comments:**
```bash
doxygen Doxyfile      # run from the repo root; rewrites docs/api/ only
```
Install Doxygen first if needed: `brew install doxygen`.

`docs/` is laid out as:

| Path | Contents | Generated? |
|---|---|---|
| `docs/index.html` | redirect that keeps the Pages root URL working | hand-written |
| `docs/api/` | Doxygen API reference | **yes — overwritten by `doxygen Doxyfile`** |
| `docs/revision/` | manuscript-revision working documents (see its `README.md`) | hand-written |

## License

Released under the MIT License. See [`LICENSE`](LICENSE) for the full text.
