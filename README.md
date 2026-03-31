# bi-kappa

Monte Carlo particle sampling project with:

- C++ generators for velocity and position distributions
- Built-in transform tests for field-aligned to global-frame rotation
- Python generators for custom velocity/position sampling
- Jupyter notebook visualizations for generated samples

## Project Layout

- `cpp/`
	- `bi_kappa_distribution.H`: bi-kappa velocity sampler (with frame transform)
	- `bi_maxwellian_distribution.H`: bi-Maxwellian velocity sampler (with frame transform)
	- `general_velocity_generator.H`: rejection sampler for user-defined energy PDF
	- `general_position_generator.H`: rejection sampler for user-defined spatial density
	- `main.cpp`: transform self-tests + sample generation driver
	- `GNUmakefile`: Makefile
- `python/`
	- `general_generators.py`: Python implementations of general velocity/position generators
	- `visualize_samples.ipynb`: analysis and visualization notebook
	- `pyproject.toml`: Python dependencies

## C++ Generators

### bi_kappa_distribution

- Samples 3D velocities from a bi-kappa distribution
- Parameters include `kappa`, `theta_perp`, `theta_par`, and field direction `ub`
- Applies a local-field-frame to global-frame transform via
	`rotate_from_fieldAligned_frame(...)`

### bi_maxwellian_distribution

- Samples 3D velocities from an anisotropic Gaussian (bi-Maxwellian)
- Parameters include `theta_perp`, `theta_par`, and field direction `ub`
- Also applies the same frame transform to map local components into global coordinates

### general_velocity_generator

- Samples energy using rejection sampling from user-defined `f(E)`
- Converts sampled energy to speed and assigns isotropic direction in 3D

### general_position_generator

- Samples positions from user-defined density using rejection sampling in a box domain
- Supports 1D/2D/3D, with output padded to 3 components when needed

## C++ Usage Examples

```cpp
#include <array>
#include <cmath>
#include <random>
#include <vector>
#include "bi_kappa_distribution.H"
#include "bi_maxwellian_distribution.H"
#include "general_velocity_generator.H"
#include "general_position_generator.H"

int main() {
    typedef double Real;
    std::mt19937 gen(12345u);

    // 1) bi_kappa_distribution
    bi_kappa_distribution<Real>::point_type ub1 = {0.0, 0.0, 1.0};
    bi_kappa_distribution<Real> bikappa(2.0, 1.0, 2.0, ub1);
    bi_kappa_distribution<Real>::point_type vk = bikappa(gen);

    // 2) bi_maxwellian_distribution
    bi_maxwellian_distribution<Real>::point_type ub2 = {1.0, 0.0, 0.0};
    bi_maxwellian_distribution<Real> bimaxwell(1.0, 2.0, ub2);
    bi_maxwellian_distribution<Real>::point_type vm = bimaxwell(gen);

    // 3) general_velocity_generator: f(E) = exp(-E), E in [0, 20], mass = 1
    auto energyPdf = [](Real E) -> Real { return std::exp(-E); };
    general_velocity_generator<Real> velGen(energyPdf, 0.0, 20.0, 1.0);
    general_velocity_generator<Real>::point_type vg = velGen(gen);

    // 4) general_position_generator: rho(x, y) = 1 + sin(x) sin(y), domain [-pi, pi]^2
    typedef general_position_generator<Real>::point_type position_point_type;
    auto rho = [](const position_point_type &x) -> Real {
        return 1.0 + std::sin(x[0]) * std::sin(x[1]);
    };
    const Real pi = 3.14159265358979323846;
    general_position_generator<Real> posGen(2, {-pi, -pi}, {pi, pi}, rho);
    position_point_type x = posGen(gen);

    return 0;
}
```

### Notes for the Examples

- `ub` can be non-normalized; it is normalized internally by the transform.
- If `ub` is left at the default zero vector, transform is skipped and local-frame components are returned directly.
- `general_position_generator` returns a `std::vector<Real>` of length equal to `dimension`.
- Use your project `Real` type if needed (for example from Chombo), or replace `double` above.

## C++ Transform Tests

`cpp/main.cpp` runs transform self-tests before generating samples.

Current checks include:

- component decomposition consistency (parallel and perpendicular parts)
- axis mapping for known `ub` direction
- dot-product preservation (orthogonality property)
- linearity check `T(a+b)=T(a)+T(b)` (T denotes transform)

If a test fails, the executable exits with non-zero status.

## Build and Run C++

### 1) Prerequisites

- Compiler that supports C++11.
- In cpp/GNUmakefile, set `CXX=/path/to/your/compiler`(g++ by default).

### 2) Build

```bash
cd cpp
make
```

This produces an executable:

- `main.exe`

### 3) Run

```bash
cd cpp
./main.exe
```

Outputs include:

- `samples_bikappa.txt`
- `samples_bimaxwellian.txt`
- `samples_general_velocity.txt`
- `samples_general_position.txt`

## Python Generators

`python/general_generators.py` provides:

- `GeneralVelocityGenerator`
- `GeneralPositionGenerator`

and a `main()` that writes:

- `samples_general_velocity.txt`
- `samples_general_position.txt`

from the current working directory.

## Python Setup

From `python/`, dependencies are defined in `pyproject.toml`.

Option A (uv):

```bash
cd python
uv sync
```

Option B (pip):

```bash
cd python
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install numpy pandas matplotlib scipy ipykernel
```

## Run Python Generator Script

```bash
cd python
python general_generators.py
```

## Run Visualization Notebook

Notebook:

- `python/visualize_samples.ipynb`

Start Jupyter:

```bash
cd python
jupyter lab
```

Then open `visualize_samples.ipynb` and run all cells.

### Sample File Location for Notebook

The notebook currently uses:

```python
workspace_root = "./"
# workspace_root = "../cpp/"
```

Choose one workflow:

- Keep `workspace_root = "./"` and place sample files in `python/`
- Or set `workspace_root = "../cpp/"` to read samples generated by C++ in `cpp/`

## Typical End-to-End Workflow

1. Build and run C++ in `cpp/` to run transform tests and generate sample files.
2. Open `python/visualize_samples.ipynb`.
3. Set `workspace_root` to match where sample files exist.
4. Run all notebook cells to generate plots and diagnostics.

