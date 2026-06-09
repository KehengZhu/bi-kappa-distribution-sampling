# bi-kappa

Monte Carlo particle samplers for plasma **velocity** and **position** distributions, in C++ and Python.

- **C++** (`cpp/`) — header-only samplers: bi-kappa & bi-Maxwellian velocities, plus rejection samplers for any energy or spatial density you define.
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

## The four samplers

| Header | Samples | define(...) takes |
|---|---|---|
| `bi_kappa_distribution` | 3D bi-kappa velocity | `kappa, theta_perp, theta_par, ub, cap[, seed]` |
| `bi_maxwellian_distribution` | 3D bi-Maxwellian velocity | `theta_perp, theta_par, ub, cap` |
| `general_velocity_generator` | speed from your `f(E)`, isotropic direction | `f, E_min, E_max, mass` |
| `general_position_generator` | position from your density `rho(x)` | `dimension, lower, upper, rho` |

They all share the same behavior:

- **Construct → `define(...)` → call.** Calling before `define(...)` throws.
- **Seeding:** call `dist.seed(s)` any time. `bi_kappa_distribution` also accepts the seed as the last `define(...)` argument. A **negative** `s` draws an unpredictable seed from `std::random_device`.
- **No-arg call** `dist()` uses the sampler's own RNG. **Bring-your-own** `dist(gen)` accepts an external `std::mt19937`.
- **`ub`** is the magnetic-field direction (need not be unit length). Output is rotated into the global frame; the default `{0,0,1}` returns field-aligned components directly.
- **`cap`** (default `20`) rejects any component with `|v_i| / theta_i > cap` and resamples; throws after 10⁶ failed tries.
- **Change one parameter** without re-specifying the rest: `dist.kappa(3.0)`, `dist.ub({0,1,0})`, etc.

For the full API — every overload, parameter constraint, and the Python classes — see the [API reference](https://kehengzhu.github.io/bi-kappa-distribution-sampling/). A worked example using all four samplers lives in [`cpp/main.cpp`](cpp/main.cpp).

## Quick start — Python

```bash
cd python
uv sync                # or: pip install numpy scipy pandas matplotlib ipykernel
python general_generators.py     # writes samples_general_{velocity,position}.txt
```

`general_generators.py` provides `GeneralVelocityGenerator` and `GeneralPositionGenerator` — the Python equivalents of the C++ general samplers.

## Visualize

Open [`python/visualize_samples.ipynb`](python/visualize_samples.ipynb) and run all cells. Point `workspace_root` at your sample files:

- `"./"` — files in `python/`
- `"../cpp/"` — files written by the C++ driver

## Tests

`./main.exe` runs `run_all_tests()` before sampling and **exits non-zero on any failure**. Coverage: frame-transform correctness, `define(...)`/getter behavior, the velocity cap, and seed reproducibility.

## Documentation

The full API reference is generated from the in-source comments with [Doxygen](https://www.doxygen.nl/) and lives in [`docs/`](docs/).

**Read it locally** — open [`docs/index.html`](docs/index.html) in a browser. (Note: clicking that file on github.com shows the HTML *source*, not the rendered page — use GitHub Pages below for a live site.)

**Regenerate after changing code or comments:**
```bash
doxygen Doxyfile      # run from the repo root; rewrites docs/
```
Install Doxygen first if needed: `brew install doxygen`.

## License

Released under the MIT License. See [`LICENSE`](LICENSE) for the full text.
