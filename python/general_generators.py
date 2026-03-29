import numpy as np
from typing import Callable, List


class GeneralVelocityGenerator:
    """Samples velocities from a user-defined energy distribution f(E)."""
    
    def __init__(self, energy_pdf: Callable, energy_min: float, energy_max: float, 
                 particle_mass: float, probe_points: int = 1024, max_reject_tries: int = 100000):
        self.energy_pdf = energy_pdf
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.particle_mass = particle_mass
        self.max_reject_tries = max_reject_tries
        
        # Estimate upper bound of the PDF
        self.pdf_upper_bound = self._estimate_pdf_upper_bound(probe_points)
    
    def _estimate_pdf_upper_bound(self, probe_points: int) -> float:
        """Estimate the maximum value of the energy PDF in the given range."""
        energies = np.linspace(self.energy_min, self.energy_max, probe_points)
        pdf_values = np.array([self.energy_pdf(E) for E in energies])
        return np.max(pdf_values)
    
    def sample_energy(self, rng: np.random.Generator) -> float:
        """Sample energy using rejection sampling."""
        for _ in range(self.max_reject_tries):
            E = rng.uniform(self.energy_min, self.energy_max)
            u = rng.uniform(0, self.pdf_upper_bound)
            if u <= self.energy_pdf(E):
                return E
        raise RuntimeError("Failed to sample energy after max_reject_tries attempts")
    
    def __call__(self, rng: np.random.Generator) -> np.ndarray:
        """Generate a random 3D velocity vector."""
        energy = self.sample_energy(rng)
        speed = np.sqrt(2.0 * energy / self.particle_mass)
        
        # Isotropic direction on sphere
        cos_theta = rng.uniform(-1.0, 1.0)
        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta**2))
        phi = rng.uniform(0.0, 2.0 * np.pi)
        
        vx = speed * sin_theta * np.cos(phi)
        vy = speed * sin_theta * np.sin(phi)
        vz = speed * cos_theta
        
        return np.array([vx, vy, vz])


class GeneralPositionGenerator:
    """Samples positions in 1D/2D/3D from a user density function using rejection sampling."""
    
    def __init__(self, dimension: int, lower_bounds: List[float], upper_bounds: List[float],
                 density_function: Callable, probe_points: int = 4096, max_reject_tries: int = 200000):
        self.dimension = dimension
        self.lower_bounds = np.array(lower_bounds)
        self.upper_bounds = np.array(upper_bounds)
        self.density_function = density_function
        self.max_reject_tries = max_reject_tries
        
        # Estimate upper bound of the density
        self.density_upper_bound = self._estimate_density_upper_bound(probe_points)
    
    def _estimate_density_upper_bound(self, probe_points: int) -> float:
        """Estimate the maximum value of the density function in the domain."""
        max_density = 0.0
        for _ in range(probe_points):
            x = self.lower_bounds + np.random.rand(self.dimension) * (self.upper_bounds - self.lower_bounds)
            rho = self.density_function(x.tolist())
            max_density = max(max_density, rho)
        return max_density
    
    def __call__(self, rng: np.random.Generator) -> np.ndarray:
        """Generate a random position vector using rejection sampling."""
        for _ in range(self.max_reject_tries):
            # Sample uniformly in the domain
            x = self.lower_bounds + rng.uniform(0, 1, self.dimension) * (self.upper_bounds - self.lower_bounds)
            u = rng.uniform(0, self.density_upper_bound)
            
            if u <= self.density_function(x.tolist()):
                # Pad to 3D if necessary
                if self.dimension < 3:
                    x = np.pad(x, (0, 3 - self.dimension), mode='constant', constant_values=0)
                return x
        
        raise RuntimeError("Failed to sample position after max_reject_tries attempts")


def main():
    """Generate samples from General Velocity and General Position distributions."""
    n_particle = 1000000
    rng = np.random.default_rng()
    
    print("=== Example 1: GeneralVelocityGenerator f(E)=exp(-E) ===")
    energy_pdf = lambda E: np.exp(-E)
    vel_gen = GeneralVelocityGenerator(energy_pdf, 0.0, 20.0, 1.0)
    
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
    rho = lambda x: (1+np.sin(x[0])*np.sin(x[1])*np.sin(x[2])) if len(x) >= 3 else (1+np.sin(x[0])*np.sin(x[1]))
    
    # 3D domain
    from numpy import pi
    pos_gen = GeneralPositionGenerator(2, [-pi, -pi], [pi, pi], rho)
    
    samples_position = []
    for i in range(n_particle):
        x = pos_gen(rng)
        samples_position.append(x)
        if (i + 1) % 100000 == 0:
            print(f"  generated {i + 1} samples")
    
    # Write to file
    with open("samples_general_position.txt", "w") as f:
        for x in samples_position:
            f.write(f"{x[0]} {x[1]} {x[2]}\n")
    print(f"Wrote {n_particle} samples to samples_general_position.txt")


if __name__ == "__main__":
    main()
