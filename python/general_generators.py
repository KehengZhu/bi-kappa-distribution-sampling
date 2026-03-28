"""General particle generators for velocity and position sampling.

This module provides two rejection-sampling based generators:
- GeneralVelocityGenerator: sample 3D velocity from an arbitrary energy PDF f(E)
- GeneralPositionGenerator: sample 1D/2D/3D position from an arbitrary density rho(x)

Example:
    import random

    def energy_pdf(e):
        return e * (2.718281828459045 ** (-e))

    vg = GeneralVelocityGenerator(
        energy_pdf=energy_pdf,
        energy_min=0.0,
        energy_max=20.0,
        particle_mass=1.0,
    )

    rng = random.Random(42)
    v = vg.sample(rng)
    print("Velocity:", v)

    def rho_2d(pos):
        x, y = pos
        return 1.0 + 0.25 * x * x + 0.1 * y

    pg = GeneralPositionGenerator(
        dimension=2,
        lower_bounds=[-1.0, -1.0],
        upper_bounds=[1.0, 1.0],
        density_function=rho_2d,
    )

    p = pg.sample(rng)
    print("Position:", p)
"""

import math
import random
from typing import Callable, List, Sequence


class GeneralVelocityGenerator:
    """Generate isotropic 3D velocities from an arbitrary energy distribution.

    The sampled velocity satisfies E = m * |V|^2 / 2.

    Args:
        energy_pdf: Callable f(E) >= 0 on [energy_min, energy_max].
        energy_min: Lower energy bound (inclusive).
        energy_max: Upper energy bound (inclusive).
        particle_mass: Particle mass (must be positive).
        probe_points: Number of probe points to estimate PDF upper bound.
        max_reject_tries: Max iterations for rejection sampling.
    """

    def __init__(
        self,
        energy_pdf: Callable[[float], float],
        energy_min: float,
        energy_max: float,
        particle_mass: float,
        probe_points: int = 1024,
        max_reject_tries: int = 100000,
    ) -> None:
        if not callable(energy_pdf):
            raise ValueError("energy_pdf must be callable")
        if not (energy_min < energy_max):
            raise ValueError("energy_min must be strictly less than energy_max")
        if particle_mass <= 0.0:
            raise ValueError("particle_mass must be positive")
        if probe_points <= 1:
            raise ValueError("probe_points must be greater than 1")
        if max_reject_tries <= 0:
            raise ValueError("max_reject_tries must be positive")

        self.energy_pdf = energy_pdf
        self.energy_min = float(energy_min)
        self.energy_max = float(energy_max)
        self.particle_mass = float(particle_mass)
        self.max_reject_tries = int(max_reject_tries)
        self._pdf_upper_bound = self._estimate_pdf_upper_bound(probe_points)

        if self._pdf_upper_bound <= 0.0:
            raise ValueError("energy_pdf must be positive somewhere on [energy_min, energy_max]")

    def _eval_pdf(self, energy: float) -> float:
        value = float(self.energy_pdf(energy))
        if not math.isfinite(value):
            raise ValueError("energy_pdf returned non-finite value")
        if value < 0.0:
            raise ValueError("energy_pdf returned a negative value")
        return value

    def _estimate_pdf_upper_bound(self, probe_points: int) -> float:
        max_value = 0.0
        width = self.energy_max - self.energy_min
        for i in range(probe_points):
            ratio = i / float(probe_points - 1)
            energy = self.energy_min + ratio * width
            value = self._eval_pdf(energy)
            if value > max_value:
                max_value = value
        return 1.05 * max_value

    def _sample_energy(self, rng: random.Random) -> float:
        for _ in range(self.max_reject_tries):
            energy = rng.uniform(self.energy_min, self.energy_max)
            y = rng.uniform(0.0, self._pdf_upper_bound)
            if y <= self._eval_pdf(energy):
                return energy
        raise RuntimeError("failed to sample energy: rejection sampling exceeded max_reject_tries")

    def sample(self, rng: random.Random = None) -> List[float]:
        """Sample a single 3D velocity [vx, vy, vz]."""
        if rng is None:
            rng = random.Random()

        energy = self._sample_energy(rng)
        speed = math.sqrt(2.0 * energy / self.particle_mass)

        cos_theta = rng.uniform(-1.0, 1.0)
        sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
        phi = rng.uniform(0.0, 2.0 * math.pi)

        vx = speed * sin_theta * math.cos(phi)
        vy = speed * sin_theta * math.sin(phi)
        vz = speed * cos_theta
        return [vx, vy, vz]


class GeneralPositionGenerator:
    """Generate positions in 1D/2D/3D from an arbitrary spatial density.

    Args:
        dimension: Number of dimensions (1, 2, or 3).
        lower_bounds: Lower bounds for each axis.
        upper_bounds: Upper bounds for each axis.
        density_function: Callable rho(position) >= 0.
        probe_points: Number of random probe points for estimating max density.
        max_reject_tries: Max iterations for rejection sampling.
    """

    def __init__(
        self,
        dimension: int,
        lower_bounds: Sequence[float],
        upper_bounds: Sequence[float],
        density_function: Callable[[Sequence[float]], float],
        probe_points: int = 4096,
        max_reject_tries: int = 200000,
    ) -> None:
        if dimension not in (1, 2, 3):
            raise ValueError("dimension must be 1, 2, or 3")
        if not callable(density_function):
            raise ValueError("density_function must be callable")
        if len(lower_bounds) != dimension or len(upper_bounds) != dimension:
            raise ValueError("bound vector lengths must match dimension")
        if probe_points <= 0:
            raise ValueError("probe_points must be positive")
        if max_reject_tries <= 0:
            raise ValueError("max_reject_tries must be positive")

        self.dimension = int(dimension)
        self.lower_bounds = [float(v) for v in lower_bounds]
        self.upper_bounds = [float(v) for v in upper_bounds]
        self.density_function = density_function
        self.max_reject_tries = int(max_reject_tries)

        for i in range(self.dimension):
            if not (self.lower_bounds[i] < self.upper_bounds[i]):
                raise ValueError("each lower bound must be strictly less than upper bound")

        self._density_upper_bound = self._estimate_density_upper_bound(probe_points)
        if self._density_upper_bound <= 0.0:
            raise ValueError("density_function must be positive somewhere in the domain")

    def _eval_density(self, position: Sequence[float]) -> float:
        value = float(self.density_function(position))
        if not math.isfinite(value):
            raise ValueError("density_function returned non-finite value")
        if value < 0.0:
            raise ValueError("density_function returned negative value")
        return value

    def _sample_uniform_point(self, rng: random.Random) -> List[float]:
        return [
            rng.uniform(self.lower_bounds[i], self.upper_bounds[i])
            for i in range(self.dimension)
        ]

    def _estimate_density_upper_bound(self, probe_points: int) -> float:
        probe_rng = random.Random(1337)
        max_value = 0.0
        for _ in range(probe_points):
            point = self._sample_uniform_point(probe_rng)
            value = self._eval_density(point)
            if value > max_value:
                max_value = value
        return 1.05 * max_value

    def sample(self, rng: random.Random = None) -> List[float]:
        """Sample a single position vector with length = dimension."""
        if rng is None:
            rng = random.Random()

        for _ in range(self.max_reject_tries):
            point = self._sample_uniform_point(rng)
            y = rng.uniform(0.0, self._density_upper_bound)
            if y <= self._eval_density(point):
                return point

        raise RuntimeError("failed to sample position: rejection sampling exceeded max_reject_tries")
