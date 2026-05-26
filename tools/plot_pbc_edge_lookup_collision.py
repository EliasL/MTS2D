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
class Ghost:
    ref_id: int
    ghost_id: tuple[int, int]
    shift: tuple[int, int]
    pos: tuple[float, float]


@dataclass(frozen=True)
class Element:
    label: str
    angle_node: int
    ghosts: tuple[Ghost, Ghost, Ghost]
    color: str

    @property
    def angle_edge_indices(self) -> tuple[int, int]:
        return tuple(sorted(((self.angle_node + 1) % 3, (self.angle_node + 2) % 3)))

    @property
    def angle_edge_node_ids(self) -> tuple[int, int]:
        i, j = self.angle_edge_indices
        return tuple(sorted((self.ghosts[i].ref_id, self.ghosts[j].ref_id)))

    @property
    def points(self) -> np.ndarray:
        return np.array([g.pos for g in self.ghosts], dtype=float)


ELEMENTS = (
    Element(
        "existing eIndex=2",
        angle_node=2,
        color="tab:blue",
        ghosts=(
            Ghost(ref_id=1, ghost_id=(1, 0), shift=(0, 0), pos=(1, 1)),
            Ghost(ref_id=0, ghost_id=(2, 0), shift=(2, 0), pos=(2, 0)),
            Ghost(ref_id=3, ghost_id=(1, 1), shift=(0, 0), pos=(2, 1)),
        ),
    ),
    Element(
        "existing eIndex=7",
        angle_node=2,
        color="tab:orange",
        ghosts=(
            Ghost(ref_id=0, ghost_id=(2, 2), shift=(2, 2), pos=(4, 2)),
            Ghost(ref_id=1, ghost_id=(1, 2), shift=(0, 2), pos=(3, 3)),
            Ghost(ref_id=2, ghost_id=(2, 1), shift=(2, 0), pos=(3, 2)),
        ),
    ),
    Element(
        "new eIndex=0",
        angle_node=0,
        color="tab:green",
        ghosts=(
            Ghost(ref_id=3, ghost_id=(-1, 1), shift=(-2, 0), pos=(0, 1)),
            Ghost(ref_id=0, ghost_id=(0, 0), shift=(0, 0), pos=(0, 0)),
            Ghost(ref_id=1, ghost_id=(1, 0), shift=(0, 0), pos=(1, 1)),
        ),
    ),
)


def draw_element(ax: plt.Axes, element: Element) -> None:
    pts = element.points
    closed = np.vstack([pts, pts[0]])
    ax.plot(closed[:, 0], closed[:, 1], color=element.color, linewidth=1.8)
    ax.fill(pts[:, 0], pts[:, 1], color=element.color, alpha=0.12)

    edge_i, edge_j = element.angle_edge_indices
    edge = pts[[edge_i, edge_j]]
    ax.plot(edge[:, 0], edge[:, 1], color=element.color, linewidth=5.0,
            solid_capstyle="round")

    angle = pts[element.angle_node]
    ax.scatter([angle[0]], [angle[1]], color=element.color, marker="x", s=90,
               linewidths=2.5)

    center = pts.mean(axis=0)
    ax.text(center[0], center[1], element.label, color=element.color,
            ha="center", va="center", fontsize=9,
            bbox={"facecolor": "white", "edgecolor": element.color, "alpha": 0.82})

    for k, ghost in enumerate(element.ghosts):
        x, y = ghost.pos
        ax.scatter([x], [y], color=element.color, s=42)
        ax.text(
            x + 0.035,
            y + 0.045,
            f"g{k}: node {ghost.ref_id}\n"
            f"gid {ghost.ghost_id}\n"
            f"shift {ghost.shift}",
            fontsize=8,
            color=element.color,
            ha="left",
            va="bottom",
        )


def main() -> None:
    output_dir = Path("tmp")
    output_dir.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5, 6.5), constrained_layout=True)
    for element in ELEMENTS:
        draw_element(ax, element)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("current x")
    ax.set_ylabel("current y")
    ax.set_title(
        "PBC edge lookup collision: three ghost placements share EdgeKey(node 0, node 1)"
    )
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.set_xlim(-0.35, 4.55)
    ax.set_ylim(-0.35, 3.45)

    explanation = (
        "Thick segments are the current angle edges used by addToEdgeTwinLookup.\n"
        "All three edge keys are identical because the key uses real node ids only: "
        "EdgeKey(0, 1).\n"
        "The plotted ghost ids/periodic shifts show that these are different "
        "periodic placements."
    )
    ax.text(
        0.01,
        0.99,
        explanation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )

    png = output_dir / "pbc_edge_lookup_collision.png"
    svg = output_dir / "pbc_edge_lookup_collision.svg"
    fig.savefig(png, dpi=220)
    fig.savefig(svg)
    print(png.resolve())
    print(svg.resolve())


if __name__ == "__main__":
    main()
