"""NumPy versions of the general velocity and position generators.

GeneralVelocityGenerator, FieldAlignedVelocityGenerator and
GeneralPositionGenerator are small NumPy ports of the C++ classes
general_velocity_generator, field_aligned_velocity_generator and
general_position_generator in cpp/, meant for prototyping a density in
Python.  They take the same parameters and use the same rejection envelope (the
largest density found at a set of probe points, times 1.05).  They are not a
validated equivalent of the C++ classes: the C++ regression suite does not
cover them, they check fewer invalid inputs, and they do not reproduce the C++
random streams.

Each generator is called with a numpy.random.Generator and returns one draw:
``v = gen(rng)``.  Running this file writes three example sample files to the
current directory.
"""

import math
from typing import Callable, List

import numpy as np

# Fewer probes can miss the peak of g and bias the samples silently (as in the C++ headers).
MIN_PROBE_POINTS = 64


class GeneralVelocityGenerator:
    """Samples velocities from a user-defined speed-squared distribution g(w), w = |v|^2.

    The speed is |v| = sqrt(w); the direction is isotropic in 3D. (For a
    distribution originally posed in energy, w = |v|^2 = 2E/m, so the mass drops
    out of this generator entirely.)  If speed_sq_min == speed_sq_max, every sample has
    the single speed sqrt(speed_sq_min)."""

    def __init__(self, speed_sq_pdf: Callable, speed_sq_min: float, speed_sq_max: float,
                 probe_points: int = 1024, max_reject_tries: int = 100000):
        if speed_sq_min < 0.0:
            raise ValueError("speed_sq_min must be non-negative (w = |v|^2 >= 0)")
        if not (speed_sq_min <= speed_sq_max):
            raise ValueError("speed_sq_max must not be less than speed_sq_min")
        if probe_points < MIN_PROBE_POINTS:
            raise ValueError("probe_points must be at least 64")
        self.speed_sq_pdf = speed_sq_pdf
        self.speed_sq_min = speed_sq_min
        self.speed_sq_max = speed_sq_max
        self.max_reject_tries = max_reject_tries

        # Estimate upper bound of the PDF (not needed for a single-valued range).
        self.pdf_upper_bound = (self._estimate_pdf_upper_bound(probe_points)
                                if speed_sq_min < speed_sq_max else 0.0)

    def _estimate_pdf_upper_bound(self, probe_points: int) -> float:
        """Estimate the maximum value of the speed-squared PDF, with a 5% margin."""
        grid = np.linspace(self.speed_sq_min, self.speed_sq_max, probe_points)
        pdf_values = np.array([self.speed_sq_pdf(w) for w in grid])
        return np.max(pdf_values) * 1.05

    def sample_speed_sq(self, rng: np.random.Generator) -> float:
        """Sample w = |v|^2 using rejection sampling."""
        if not self.speed_sq_min < self.speed_sq_max:
            return self.speed_sq_min
        for _ in range(self.max_reject_tries):
            w = rng.uniform(self.speed_sq_min, self.speed_sq_max)
            u = rng.uniform(0, self.pdf_upper_bound)
            if u <= self.speed_sq_pdf(w):
                return w
        raise RuntimeError("Failed to sample speed_sq after max_reject_tries attempts")

    def __call__(self, rng: np.random.Generator) -> np.ndarray:
        """Generate a random 3D velocity vector."""
        speed = np.sqrt(self.sample_speed_sq(rng))

        # Isotropic direction on sphere
        cos_theta = rng.uniform(-1.0, 1.0)
        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta**2))
        phi = rng.uniform(0.0, 2.0 * np.pi)
        
        vx = speed * sin_theta * np.cos(phi)
        vy = speed * sin_theta * np.sin(phi)
        vz = speed * cos_theta
        
        return np.array([vx, vy, vz])


class GeneralPositionGenerator:
    """Samples positions in 1D/2D/3D from a user density rho(x) by rejection sampling.

    The domain is the axis-aligned box lower_bounds <= x <= upper_bounds, a hard
    support.  density_function takes a list of length ``dimension`` and returns
    a finite, non-negative value, positive somewhere in the box; it need not be
    normalized.  The rejection envelope is the largest density found at
    probe_points uniform points, drawn from a fixed seed so that two identical
    constructions give the same envelope, times a 5% margin.  A peak narrow
    enough to fall between the probes is undersampled without warning.  Each
    call returns an array of length ``dimension``.
    """

    PROBE_SEED = 1337        # fixed, as in the C++ class, so the envelope is reproducible

    def __init__(self, dimension: int, lower_bounds: List[float], upper_bounds: List[float],
                 density_function: Callable, probe_points: int = 4096, max_reject_tries: int = 200000):
        if isinstance(dimension, bool) or dimension not in (1, 2, 3):
            raise ValueError("dimension must be 1, 2, or 3")
        if not callable(density_function):
            raise ValueError("density_function must be callable")
        lower = np.asarray(lower_bounds, dtype=float)
        upper = np.asarray(upper_bounds, dtype=float)
        if lower.shape != (dimension,) or upper.shape != (dimension,):
            raise ValueError("lower_bounds and upper_bounds must each have `dimension` entries")
        if not (np.isfinite(lower).all() and np.isfinite(upper).all()):
            raise ValueError("the bounds must be finite")
        if not (lower < upper).all():
            raise ValueError("each lower bound must be strictly less than its upper bound")
        if probe_points <= 0:
            raise ValueError("probe_points must be positive")
        if max_reject_tries <= 0:
            raise ValueError("max_reject_tries must be positive")

        self.dimension = dimension
        self.lower_bounds = lower
        self.upper_bounds = upper
        self.density_function = density_function
        self.max_reject_tries = max_reject_tries

        self.density_upper_bound = self._estimate_density_upper_bound(probe_points)
        if not self.density_upper_bound > 0.0:
            raise ValueError("density_function must be positive somewhere inside the domain")

    def _density(self, x: np.ndarray) -> float:
        rho = float(self.density_function(x.tolist()))
        if not math.isfinite(rho):
            raise ValueError("density_function returned a non-finite value")
        if rho < 0.0:
            raise ValueError("density_function returned a negative value; expected a "
                             "non-negative density")
        return rho

    def _estimate_density_upper_bound(self, probe_points: int) -> float:
        """Largest density at probe_points uniform points in the box, with a 5% margin."""
        probe_rng = np.random.default_rng(self.PROBE_SEED)
        width = self.upper_bounds - self.lower_bounds
        max_density = 0.0
        for _ in range(probe_points):
            x = self.lower_bounds + probe_rng.uniform(0.0, 1.0, self.dimension) * width
            max_density = max(max_density, self._density(x))
        return max_density * 1.05

    def __call__(self, rng: np.random.Generator) -> np.ndarray:
        """Generate one position, an array of length ``dimension``."""
        width = self.upper_bounds - self.lower_bounds
        for _ in range(self.max_reject_tries):
            x = self.lower_bounds + rng.uniform(0.0, 1.0, self.dimension) * width
            u = rng.uniform(0.0, self.density_upper_bound)
            if u <= self._density(x):
                return x

        raise RuntimeError("Failed to sample position after max_reject_tries attempts")


class FieldAlignedVelocityGenerator:
    """Samples field-aligned velocities: an arbitrary parallel speed-squared
    distribution g(w_par), w_par = v_par^2, along B, with a Maxwellian
    perpendicular spread (thermal speed theta_perp). (For a distribution
    originally posed in energy, w_par = v_par^2 = 2E_par/m, so the mass drops
    out.)  If speed_sq_min == speed_sq_max, every sample has the single parallel speed
    sqrt(speed_sq_min).  Each perpendicular component is capped at
    max_normalized_velocity * theta_perp (default 20, as in the C++ header); pass
    math.inf for no cap."""

    def __init__(self, speed_sq_pdf: Callable, speed_sq_min: float, speed_sq_max: float,
                 theta_perp: float, ub: List[float], parallel_sign: int,
                 probe_points: int = 1024, max_reject_tries: int = 100000,
                 max_normalized_velocity: float = 20.0):
        if speed_sq_min < 0.0:
            raise ValueError("speed_sq_min must be non-negative (w = v_par^2 >= 0)")
        if not (speed_sq_min <= speed_sq_max):
            raise ValueError("speed_sq_max must not be less than speed_sq_min")
        if probe_points < MIN_PROBE_POINTS:
            raise ValueError("probe_points must be at least 64")
        if not max_normalized_velocity > 0.0:
            raise ValueError("max_normalized_velocity must be positive, or math.inf for no cap")
        if theta_perp <= 0.0:
            raise ValueError("theta_perp must be positive")
        if parallel_sign not in (1, -1):
            raise ValueError("parallel_sign must be exactly +1 or -1")
        ub = np.asarray(ub, dtype=float)
        if np.dot(ub, ub) == 0.0:
            raise ValueError("ub (magnetic field) must be a non-zero vector")

        self.speed_sq_pdf = speed_sq_pdf
        self.speed_sq_min = speed_sq_min
        self.speed_sq_max = speed_sq_max
        self.theta_perp = theta_perp
        self.parallel_sign = float(parallel_sign)
        self.max_normalized_velocity = max_normalized_velocity
        self.max_reject_tries = max_reject_tries

        # Field-aligned orthonormal frame (e1, e2 perpendicular; e3 = B/|B|).
        self.e1, self.e2, self.e3 = self._build_frame(ub)

        # Estimate upper bound of the parallel speed-squared PDF.
        self.pdf_upper_bound = (self._estimate_pdf_upper_bound(probe_points)
                                if speed_sq_min < speed_sq_max else 0.0)

    @staticmethod
    def _build_frame(ub: np.ndarray):
        """Right-handed orthonormal frame with e3 aligned to ub."""
        e3 = ub / np.linalg.norm(ub)
        # Pick a seed vector least aligned with e3 to stay well-conditioned.
        seed = np.zeros(3)
        seed[np.argmin(np.abs(e3))] = 1.0
        e1 = np.cross(e3, seed)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        return e1, e2, e3

    def _estimate_pdf_upper_bound(self, probe_points: int) -> float:
        grid = np.linspace(self.speed_sq_min, self.speed_sq_max, probe_points)
        pdf_values = np.array([self.speed_sq_pdf(w) for w in grid])
        return np.max(pdf_values) * 1.05

    def sample_speed_sq(self, rng: np.random.Generator) -> float:
        """Sample w_par = v_par^2 using rejection sampling."""
        if not self.speed_sq_min < self.speed_sq_max:
            return self.speed_sq_min
        for _ in range(self.max_reject_tries):
            w = rng.uniform(self.speed_sq_min, self.speed_sq_max)
            u = rng.uniform(0, self.pdf_upper_bound)
            if u <= self.speed_sq_pdf(w):
                return w
        raise RuntimeError("Failed to sample speed_sq after max_reject_tries attempts")

    def __call__(self, rng: np.random.Generator) -> np.ndarray:
        """Generate a random 3D velocity vector in the global frame."""
        v_par = self.parallel_sign * np.sqrt(self.sample_speed_sq(rng))

        # Perpendicular 2D Maxwellian: N(0, theta_perp^2/2) per component, with the pair
        # redrawn while either component exceeds the cap.
        sigma_perp = self.theta_perp / np.sqrt(2.0)
        perp_max = self.max_normalized_velocity * self.theta_perp
        for _ in range(self.max_reject_tries):
            v_perp1 = sigma_perp * rng.standard_normal()
            v_perp2 = sigma_perp * rng.standard_normal()
            if abs(v_perp1) <= perp_max and abs(v_perp2) <= perp_max:
                break
        else:
            raise RuntimeError("Failed to sample the perpendicular components under "
                               "max_normalized_velocity after max_reject_tries attempts")

        return v_perp1 * self.e1 + v_perp2 * self.e2 + v_par * self.e3


def main():
    """Write example samples from the three generators to text files."""
    n_particle = 200000
    rng = np.random.default_rng(20030410)

    theta_perp = 1.0

    # g(w) = sqrt(w) exp(-w/theta^2) with w = |v|^2 is the speed-squared density of an
    # isotropic Maxwellian of thermal speed theta: dw = 2|v| d|v| turns the sqrt(w) into
    # the |v|^2 of the Maxwell speed pdf, so each Cartesian component is N(0, theta^2/2).
    print("=== Example 1: GeneralVelocityGenerator g(w)=sqrt(w)exp(-w/theta^2), w=|v|^2 ===")
    theta_iso = theta_perp
    speed_sq_pdf = lambda w: np.sqrt(w) * np.exp(-w / theta_iso**2)
    vel_gen = GeneralVelocityGenerator(speed_sq_pdf, 0.0, 25.0)

    samples_velocity = []
    for i in range(n_particle):
        v = vel_gen(rng)
        samples_velocity.append(v)
        if (i + 1) % 100000 == 0:
            print(f"  generated {i + 1} samples")
    
    # Write to file
    with open("samples_general_velocity.txt", "w") as f:
        for v in samples_velocity:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
    print(f"Wrote {n_particle} samples to samples_general_velocity.txt\n")
    
    print("=== Example 2: GeneralPositionGenerator rho=1+sin(x)sin(y) ===")
    rho = lambda x: 1.0 + np.sin(x[0]) * np.sin(x[1])

    # 2D domain [-pi, pi]^2
    pos_gen = GeneralPositionGenerator(2, [-np.pi, -np.pi], [np.pi, np.pi], rho)
    
    samples_position = []
    for i in range(n_particle):
        x = pos_gen(rng)
        samples_position.append(x)
        if (i + 1) % 100000 == 0:
            print(f"  generated {i + 1} samples")
    
    # Write to file
    with open("samples_general_position.txt", "w") as f:
        for x in samples_position:
            f.write(f"{x[0]} {x[1]}\n")
    print(f"Wrote {n_particle} samples to samples_general_position.txt\n")

    # The simplest possible input: g(w) = 1 makes w = v_par^2 uniform on [0, w_max],
    # and dw = 2 v_par d(v_par) turns that into the linear ramp
    # p(v_par) = 2 v_par / w_max on [0, sqrt(w_max)]; B is tilted off the axes to
    # exercise the rotation.
    print("=== Example 3: FieldAlignedVelocityGenerator g(w)=1, B=(1,1,1), + ===")
    w_max_par = 4.0    # v_par spans [0, 2]
    parallel_speed_sq_pdf = lambda _w: 1.0
    fa_gen = FieldAlignedVelocityGenerator(parallel_speed_sq_pdf, 0.0, w_max_par,
                                           theta_perp=theta_perp, ub=[1.0, 1.0, 1.0],
                                           parallel_sign=+1)

    samples_field_aligned = []
    for i in range(n_particle):
        v = fa_gen(rng)
        samples_field_aligned.append(v)
        if (i + 1) % 100000 == 0:
            print(f"  generated {i + 1} samples")

    # Write to file
    with open("samples_field_aligned.txt", "w") as f:
        for v in samples_field_aligned:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
    print(f"Wrote {n_particle} samples to samples_field_aligned.txt")


if __name__ == "__main__":
    main()
