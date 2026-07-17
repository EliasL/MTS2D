#!/usr/bin/env python3
"""Plot the macroscopic compatibility test for the two history measures."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


HERE = Path(__file__).resolve().parent
SIMULATION_SCRIPTS = HERE.parents[2] / "SimulationScripts"
sys.path.insert(0, str(SIMULATION_SCRIPTS))

from MTMath.poincareEnergy import C2PoincareDisk  # noqa: E402


DATA = np.genfromtxt(
    HERE / "Ftotal_vs_Tfull_70x70_simple_shear.csv",
    delimiter=",",
    names=True,
)
BLUE = "#1769aa"
ORANGE = "#d95f02"
BLACK = "#252525"


def matrix_series(prefix: str) -> np.ndarray:
    return np.stack(
        [
            DATA[f"{prefix}11_mean"],
            DATA[f"{prefix}12_mean"],
            DATA[f"{prefix}21_mean"],
            DATA[f"{prefix}22_mean"],
        ],
        axis=1,
    ).reshape(-1, 2, 2)


def poincare(matrices: np.ndarray) -> np.ndarray:
    result = []
    for matrix in matrices:
        x, y = C2PoincareDisk(matrix.T @ matrix)
        result.append((float(x), float(y)))
    return np.asarray(result)


def draw_sheared_mesh(ax: plt.Axes, gamma: float = 1.0, n: int = 5) -> None:
    transform = np.array([[1.0, gamma], [0.0, 1.0]])
    values = np.linspace(0.0, 1.0, n + 1)
    for value in values:
        horizontal = np.array([[0.0, value], [1.0, value]]) @ transform.T
        vertical = np.array([[value, 0.0], [value, 1.0]]) @ transform.T
        ax.plot(horizontal[:, 0], horizontal[:, 1], color=BLUE, linewidth=0.9)
        ax.plot(vertical[:, 0], vertical[:, 1], color=BLUE, linewidth=0.9)
    for i in range(n):
        for j in range(n):
            diagonal = np.array(
                [[i / n, j / n], [(i + 1) / n, (j + 1) / n]]
            ) @ transform.T
            ax.plot(diagonal[:, 0], diagonal[:, 1], color=BLUE, linewidth=0.55)
    boundary = np.array([[0, 0], [1, 0], [2, 1], [1, 1], [0, 0]], dtype=float)
    ax.plot(boundary[:, 0], boundary[:, 1], color=BLACK, linewidth=1.8)
    ax.annotate(
        "top displacement = 1",
        xy=(1.55, 1.04), xytext=(0.55, 1.24),
        arrowprops={"arrowstyle": "->", "color": BLACK, "linewidth": 1.0},
        ha="center", fontsize=8,
    )
    ax.text(
        1.02, -0.17,
        r"$\overline{F}_{12}=1$ at $\gamma=1$",
        ha="center", fontsize=9,
    )
    ax.set_xlim(-0.08, 2.12)
    ax.set_ylim(-0.25, 1.34)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Imposed periodic simple shear", fontsize=10)


gamma = DATA["gamma"]
mathcal_f = matrix_series("Ftotal")
t_total = matrix_series("Tfull")
expected = np.stack(
    [np.ones_like(gamma), gamma, np.zeros_like(gamma), np.ones_like(gamma)],
    axis=1,
).reshape(-1, 2, 2)

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.75), constrained_layout=True)
draw_sheared_mesh(axes[0])

axes[1].plot(gamma, gamma, color=BLACK, linestyle=":", linewidth=2.1,
             label=r"required: $\langle F_{12}\rangle=\gamma$")
axes[1].plot(gamma, mathcal_f[:, 0, 1], color=BLUE, linewidth=2.1,
             label=r"$\langle\mathcal{F}_{12}\rangle$")
axes[1].plot(gamma, t_total[:, 0, 1], color=ORANGE, linewidth=2.1,
             label=r"$\langle(T_{\rm total})_{12}\rangle$")
axes[1].scatter([1.0, 1.0], [mathcal_f[-1, 0, 1], t_total[-1, 0, 1]],
                color=[BLUE, ORANGE], s=35, zorder=5)
axes[1].annotate("1.005", (1.0, mathcal_f[-1, 0, 1]), xytext=(-37, 5),
                 textcoords="offset points", color=BLUE, fontsize=8)
axes[1].annotate("0.335", (1.0, t_total[-1, 0, 1]), xytext=(-37, -12),
                 textcoords="offset points", color=ORANGE, fontsize=8)
axes[1].set_xlabel(r"imposed shear $\gamma$")
axes[1].set_ylabel("mean total shear")
axes[1].set_title("Compatibility with the boundary map", fontsize=10)
axes[1].grid(alpha=0.2)
axes[1].legend(frameon=False, fontsize=7.5, loc="upper left")

expected_pc = poincare(expected)
mathcal_f_pc = poincare(mathcal_f)
t_total_pc = poincare(t_total)
axes[2].add_patch(Circle((0, 0), 1, fill=False, edgecolor="0.45", linewidth=1.0))
axes[2].plot(expected_pc[:, 0], expected_pc[:, 1], color=BLACK,
             linestyle=":", linewidth=2.1, label="required shear path")
axes[2].plot(mathcal_f_pc[:, 0], mathcal_f_pc[:, 1], color=BLUE,
             linewidth=2.1, label=r"mean $\mathcal{F}$")
axes[2].plot(t_total_pc[:, 0], t_total_pc[:, 1], color=ORANGE,
             linewidth=2.1, label=r"mean $T_{\rm total}$")
for coordinates, color in ((expected_pc, BLACK), (mathcal_f_pc, BLUE), (t_total_pc, ORANGE)):
    axes[2].scatter(*coordinates[-1], color=color, s=34, zorder=5)
axes[2].annotate(r"$\gamma=1$", expected_pc[-1], xytext=(-29, 8),
                 textcoords="offset points", fontsize=8, color=BLACK)
axes[2].annotate(r"$T_{\rm total}$", t_total_pc[-1], xytext=(5, 1),
                 textcoords="offset points", fontsize=8, color=ORANGE)
axes[2].axhline(0, color="0.88", linewidth=0.7)
axes[2].axvline(0, color="0.88", linewidth=0.7)
axes[2].set_xlim(-0.34, 0.18)
axes[2].set_ylim(-0.05, 0.55)
axes[2].set_aspect("equal")
axes[2].set_xlabel(r"Poincaré $x_p$")
axes[2].set_ylabel(r"Poincaré $y_p$")
axes[2].set_title("Trajectory of the mean matrices", fontsize=10)
axes[2].legend(frameon=False, fontsize=7.5, loc="upper left")

fig.suptitle(
    r"A total deformation history must reproduce the imposed macroscopic map",
    fontsize=12,
)
fig.savefig(HERE / "macroscopic_compatibility_example.png", dpi=300)
fig.savefig(HERE / "macroscopic_compatibility_example.pdf")
plt.close(fig)
