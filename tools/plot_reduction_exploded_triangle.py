from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

MPL_CONFIG_DIR = Path("tmp/matplotlib").resolve()
MPL_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CONFIG_DIR))

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class GhostNode:
    label: str
    reference_id: int
    ref_pos: tuple[float, float]
    pos: tuple[float, float]


NODES = (
    GhostNode(
        "ghost[0]",
        reference_id=1,
        ref_pos=(-0.5, -0.5),
        pos=(5.99003, 1.1387799999999999),
    ),
    GhostNode(
        "ghost[1]",
        reference_id=0,
        ref_pos=(0.5, -0.5),
        pos=(5.9866200000000003, 1.1463699999999999),
    ),
    GhostNode(
        "ghost[2]",
        reference_id=3,
        ref_pos=(-0.5, 0.5),
        pos=(5.6848400000000003, 2.3721800000000002),
    ),
)


def coords(field: str) -> np.ndarray:
    return np.array([getattr(node, field) for node in NODES], dtype=float)


def signed_area(points: np.ndarray) -> float:
    a, b, c = points
    ab = b - a
    ac = c - a
    return 0.5 * (ab[0] * ac[1] - ab[1] * ac[0])


def draw_triangle(ax: plt.Axes, points: np.ndarray, title: str) -> None:
    closed = np.vstack([points, points[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="tab:blue", linewidth=2.0)
    ax.fill(points[:, 0], points[:, 1], color="tab:blue", alpha=0.12)
    ax.scatter(points[:, 0], points[:, 1], color="tab:blue", s=48, zorder=3)

    for i, (node, point) in enumerate(zip(NODES, points)):
        ax.text(
            point[0],
            point[1],
            f"  g{i}\n  node {node.reference_id}",
            fontsize=9,
            ha="left",
            va="bottom",
        )

    for i, j in ((0, 1), (1, 2), (2, 0)):
        midpoint = 0.5 * (points[i] + points[j])
        length = np.linalg.norm(points[i] - points[j])
        ax.text(
            midpoint[0],
            midpoint[1],
            f"{length:.6g}",
            fontsize=8,
            color="tab:red",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
        )

    area = signed_area(points)
    ax.set_title(f"{title}\narea={area:.16g}")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    margin = 0.12 * max(np.ptp(points[:, 0]), np.ptp(points[:, 1]), 1e-6)
    ax.set_xlim(points[:, 0].min() - margin, points[:, 0].max() + margin)
    ax.set_ylim(points[:, 1].min() - margin, points[:, 1].max() + margin)


def main() -> None:
    output_dir = Path("tmp")
    output_dir.mkdir(exist_ok=True)

    ref = coords("ref_pos")
    cur = coords("pos")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)
    draw_triangle(axes[0], ref, "Reference triangle")
    draw_triangle(axes[1], cur, "Current triangle")

    fig.suptitle(
        "Reduction explosion element eIndex=1: reference vs current triangle",
        fontsize=13,
    )

    png = output_dir / "reduction_exploded_triangle.png"
    svg = output_dir / "reduction_exploded_triangle.svg"
    fig.savefig(png, dpi=220)
    fig.savefig(svg)
    print(png.resolve())
    print(svg.resolve())


if __name__ == "__main__":
    main()
