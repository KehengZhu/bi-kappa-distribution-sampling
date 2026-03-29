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
	- `GNUmakefile`: Chombo-based build entry
- `python/`
	- `general_generators.py`: Python implementations of general velocity/position generators
	- `visualize_samples.ipynb`: analysis and visualization notebook
	- `pyproject.toml`: Python dependencies

## C++ Generators

### BiKappaDistribution

- Samples 3D velocities from a bi-kappa distribution
- Parameters include `kappa`, `theta_perp`, `theta_para`, and field direction `ub`
- Applies a local-field-frame to global-frame transform via
	`rotate_from_fieldAligned_frame(...)`

### BiMaxwellianDistribution

- Samples 3D velocities from an anisotropic Gaussian (bi-Maxwellian)
- Parameters include `theta_perp`, `theta_para`, and field direction `ub`
- Also applies the same frame transform to map local components into global coordinates

### GeneralVelocityGenerator

- Samples energy using rejection sampling from user-defined `f(E)`
- Converts sampled energy to speed and assigns isotropic direction in 3D

### GeneralPositionGenerator

- Samples positions from user-defined density using rejection sampling in a box domain
- Supports 1D/2D/3D, with output padded to 3 components when needed

## C++ Usage Examples

```cpp
#include <random>
#include "bi_kappa_distribution.H"
#include "bi_maxwellian_distribution.H"
#include "general_velocity_generator.H"
#include "general_position_generator.H"

int main() {
	typedef double Real;
    typedef std::array<Real, 3> Point3D;
	std::mt19937 gen(12345u);

	// 1) BiKappaDistribution
	BiKappaDistribution<Real>::Point3D ub1 = {0.0, 0.0, 1.0};
	BiKappaDistribution<Real> bikappa(2.0, 1.0, 2.0, ub1);
	Point3D vk = bikappa(gen); // generate one 3D vector sample

	// 2) BiMaxwellianDistribution
	BiMaxwellianDistribution<Real>::Point3D ub2 = {1.0, 0.0, 0.0};
	BiMaxwellianDistribution<Real> bimaxwell(1.0, 2.0, ub2);
	Point3D vm = bimaxwell(gen);

	// 3) GeneralVelocityGenerator: f(E) = exp(-E), E in [0, 20], mass = 1
	auto energyPdf = [](Real E) -> Real { return std::exp(-E); };
	GeneralVelocityGenerator<Real> velGen(energyPdf, 0.0, 20.0, 1.0);
	Point3D vg = velGen(gen);

	// 4) GeneralPositionGenerator: rho(x, y) = 1 + sin(x) sin(y), domain [-pi, pi]^2
	auto rho = [](const GeneralPositionGenerator<Real>::PositionVector &x) -> Real {
		return 1.0 + std::sin(x[0]) * std::sin(x[1]);
	};
	const Real pi = 3.14159265358979323846;
	GeneralPositionGenerator<Real> posGen(2, {-pi, -pi}, {pi, pi}, rho);
	GeneralPositionGenerator<Real>::PositionVector x = posGen(gen); // generate one 2D vector sample

	return 0;
}
```

### Notes for the Examples

- `ub` can be non-normalized; it is normalized internally by the transform.
- `GeneralPositionGenerator` returns a 3-component vector; in 2D, `z` is padded with `0`.
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

- Chombo installed locally
- C++ toolchain matching your Chombo build

In `cpp/GNUmakefile`, set:

```make
CHOMBO_HOME := /path/to/Chombo-3.2/lib
```

### 2) Build

```bash
cd cpp
make
```

This produces an executable similar to:

- `main3d.Linux.mpicxx.gfortran.OPT.MPI.ex`

### 3) Run

```bash
cd cpp
./main3d.Linux.mpicxx.gfortran.OPT.MPI.ex
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

