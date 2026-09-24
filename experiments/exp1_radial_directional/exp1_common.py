"""Shared definitions for the Experiment 1 cell-count test.

The direction cells, the equal-probability radius shells, and the reconstruction of the
C++ field-aligned basis, used by exp1_cells.py, exp1_tail.py and exp1_frame.py.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

# Independence test binning: radial quartiles against 40 direction cells.  The solid-angle
# element is d(cos theta) d(phi), so equal intervals of cos(theta) crossed with equal
# intervals of phi give cells of equal solid angle -- no sphere-pixelization scheme needed.
# Binning both angles is what makes this a test of the whole direction: with cos(theta)
# alone, a dependence between radius and azimuth would pass unseen.
N_COS_BINS = 5
N_PHI_BINS = 8

# The manuscript's test: Pearson's chi^2 over cells of equal probability under the target.
# N_SHELLS radius shells of equal probability, bounded by the exact deciles of R, crossed
# with the N_COS_BINS x N_PHI_BINS equal-solid-angle direction cells above.  Because the
# radius and the direction are independent under the target and the direction is uniform,
# every one of the N_SHELLS * 40 cells has probability 1/400.  The shells are compared in
# log R against log edges, so no bounded transform of R is formed and nothing underflows.
N_SHELLS = 10


def log_shell_edges(kappa: float) -> np.ndarray:
    """log of the radii that split R into N_SHELLS shells of equal probability.

    P(T <= t) = P(W >= 1/(1+t)) with W = 1/(1+T) ~ Beta(kappa-1/2, 3/2); the edge at
    cumulative probability q therefore has W_q = Beta.ppf(1-q) and R_q^2 = (1-W_q)/W_q.
    W_q stays well inside the double range for every kappa tested (R_q <= 1.4e50 at 0.51).
    """
    q = np.arange(1, N_SHELLS) / N_SHELLS
    w = stats.beta.ppf(1.0 - q, kappa - 0.5, 1.5)
    return 0.5 * (np.log1p(-w) - np.log(w))


def cell_pvalues(counts: np.ndarray) -> tuple[float, float]:
    """Pearson chi^2 p-values of (radius shells alone, all cells) against equal occupancy."""
    return (float(stats.chisquare(counts.sum(axis=1)).pvalue),
            float(stats.chisquare(counts.ravel()).pvalue))


def field_basis(ub: np.ndarray) -> np.ndarray:
    """Rebuild the orthonormal field-aligned basis exactly as bi_kappa_distribution does.

    Returns M with columns (e1, e2, e3) so that v_global = M @ v_local.  Reproducing the
    C++ construction here rather than assuming a canonical basis is the point: it is what
    makes the recovered local coordinates a genuine test of the shipped rotation.
    """
    ub = np.asarray(ub, dtype=float)
    e3 = ub / np.linalg.norm(ub)

    e2 = np.ones(3)
    maxcomp = 0
    if abs(e3[1]) > abs(e3[maxcomp]):
        maxcomp = 1
    if abs(e3[2]) > abs(e3[maxcomp]):
        maxcomp = 2
    e2[maxcomp] = 1.0 - e3.sum() / e3[maxcomp]
    e2 /= np.linalg.norm(e2)

    e1 = np.cross(e2, e3)
    return np.column_stack([e1, e2, e3])


def is_axis_aligned_z(ub: np.ndarray) -> bool:
    """Mirror the C++ short-circuit: for ub == +z the sampler skips the rotation entirely."""
    eps = np.finfo(float).eps

    def approx(a: float, b: float) -> bool:
        return abs(a - b) <= eps * max(1.0, abs(a), abs(b))

    return approx(ub[0], 0.0) and approx(ub[1], 0.0) and approx(ub[2], 1.0)
