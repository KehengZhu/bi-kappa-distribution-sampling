"""Demonstrate the directional bias from sampling the polar angle Theta
uniformly instead of cos(Theta).

A common mistake when drawing a direction on the unit sphere is to generate
the polar angle Theta uniformly. Because the spherical area element is
dOmega = sin(Theta) dTheta dPhi = -d(cos Theta) dPhi, a uniform-in-Theta
draw oversamples the poles and undersamples the equator. The correct,
area-uniform draw takes cos(Theta) ~ Uniform[-1, 1].

Figure design follows the nature-figure skill: one hero quantitative panel
(the cos-Theta marginal that proves uniformity) supported by two clean 3-D
direction clouds, a restrained red/blue signal palette, direct labels, and
editable-text PDF export.

Produces paper/overleaf/figures/sphere-sampling-bias.pdf.
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)

# Match the manuscript's other figures (matplotlib defaults: DejaVu Sans,
# tab10 colors, red dashed reference lines, alpha-0.5 histograms, light grid).
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'mathtext.fontset': 'dejavusans',
    'pdf.fonttype': 42,         # editable TrueType text in PDF
    'svg.fonttype': 'none',
    'font.size': 12,
})

# Same roles/colors as the Q--Q plots: blue = correct/target,
# orange = the deviating comparison, red dashed = analytic reference.
C_BAD = 'tab:orange'   # uniform-Theta (biased)
C_GOOD = 'tab:blue'    # uniform-cos(Theta) (correct)
C_REF = 'red'          # uniform-on-sphere reference (dashed)
C_AXIS = '#444444'     # neutral pole axis


def sample_directions_uniform_theta(n):
    """WRONG: polar angle Theta drawn uniformly in [0, pi]."""
    theta = rng.uniform(0.0, np.pi, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_theta = np.sin(theta)
    return np.column_stack((sin_theta * np.cos(phi),
                            sin_theta * np.sin(phi),
                            np.cos(theta)))


def sample_directions_uniform_costheta(n):
    """CORRECT: cos(Theta) drawn uniformly in [-1, 1] (uniform on sphere)."""
    cos_theta = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    return np.column_stack((sin_theta * np.cos(phi),
                            sin_theta * np.sin(phi),
                            cos_theta))


def _sphere_panel(ax, pts, color, title):
    """Clean 3-D direction cloud: no panes, ticks, or grid (skill stance).

    No translucent ``plot_surface`` body here on purpose: mplot3d sorts whole
    artists by a single average depth (no per-pixel z-buffer), so a
    semi-transparent surface composites *over* the scatter and looks like a
    frosted overlay. A thin wireframe plus depth-shaded points reads as a
    sphere without that artifact.
    """
    # Faint lat/long wireframe for spherical context (drawn first = behind).
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, rcount=10, ccount=16, color='0.82',
                      linewidth=0.3, alpha=0.8)

    # depthshade dims the far hemisphere, which reads as genuine 3-D.
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=2.0, alpha=0.6,
               color=color, edgecolors='none', depthshade=True,
               rasterized=True)

    ax.set_box_aspect((1, 1, 1), zoom=1.3)
    ax.view_init(elev=20, azim=-60)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.grid(False)
    ax.set_axis_off()
    # Vertical pole axis as a visual anchor.
    ax.plot([0, 0], [0, 0], [-1.02, 1.02], color=C_AXIS, lw=0.9, alpha=0.7,
            zorder=5)
    ax.text(0, 0, 1.06, 'pole', ha='center', va='bottom', fontsize=8.5,
            color=C_AXIS)
    ax.text(0, 0, -1.18, 'pole', ha='center', va='top', fontsize=8.5,
            color=C_AXIS)
    ax.set_title(title, y=1.03, fontsize=12)


def main():
    n_scatter = 4000
    n_hist = 400000

    dir_bad = sample_directions_uniform_theta(n_scatter)
    dir_good = sample_directions_uniform_costheta(n_scatter)
    cos_bad = sample_directions_uniform_theta(n_hist)[:, 2]
    cos_good = sample_directions_uniform_costheta(n_hist)[:, 2]

    fig = plt.figure(figsize=(12.0, 3.9), dpi=200)
    # Hero layout: two subordinate spheres + one wide quantitative hero.
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.5],
                          wspace=0.2, left=0.0, right=0.985,
                          bottom=0.16, top=0.90)

    ax0 = fig.add_subplot(gs[0, 0], projection='3d')
    _sphere_panel(ax0, dir_bad, C_BAD, r'Uniform $\Theta$ (biased)')

    ax1 = fig.add_subplot(gs[0, 1], projection='3d')
    _sphere_panel(ax1, dir_good, C_GOOD, r'Uniform $\cos\Theta$ (correct)')

    # ---- Hero panel: marginal density of cos(Theta) ----
    ax2 = fig.add_subplot(gs[0, 2])
    bins = np.linspace(-1.0, 1.0, 49)
    ax2.hist(cos_bad, bins=bins, density=True, alpha=0.5,
             color=C_BAD, edgecolor='black', linewidth=0.4,
             label=r'Uniform $\Theta$ (biased)')
    ax2.hist(cos_good, bins=bins, density=True, alpha=0.5,
             color=C_GOOD, edgecolor='black', linewidth=0.4,
             label=r'Uniform $\cos\Theta$ (correct)')
    ax2.axhline(0.5, color=C_REF, linestyle='--', linewidth=1.5,
                label='Uniform on sphere')

    ax2.set_xlim(-1, 1)
    ax2.set_ylim(0, 2.4)
    ax2.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax2.set_xlabel(r'$\cos\Theta = U_z\,/\,|\mathbf{U}|$')
    ax2.set_ylabel('Probability density')
    ax2.set_title(r'Marginal distribution of $\cos\Theta$', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper center', frameon=False, fontsize=10)

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.normpath(os.path.join(
        here, '..', 'paper', 'overleaf', 'figures', 'sphere-sampling-bias.pdf'))
    fig.savefig(out, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()
