#!/usr/bin/env python3
"""Compare spatially averaged Ftotal and full T components across VTUs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np


COMPONENTS = ("11", "12", "21", "22")
LOAD_PATTERN = re.compile(r"_load=([-+0-9.eE]+)_")
STEP_PATTERN = re.compile(r"\.([0-9]+)\.vtu$")


def frame_key(path: Path) -> int:
    match = STEP_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Cannot extract frame number from {path}")
    return int(match.group(1))


def load_from_path(path: Path) -> float:
    match = LOAD_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Cannot extract load from {path}")
    return float(match.group(1))


def scalar_cell_field(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.cell_data:
        raise KeyError(f"Missing VTU cell field {name!r}")
    blocks = mesh.cell_data[name]
    if len(blocks) != 1:
        raise ValueError(f"Expected one cell block for {name!r}, got {len(blocks)}")
    return np.asarray(blocks[0], dtype=float).reshape(-1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtu_directory", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    files = sorted(args.vtu_directory.glob("*.vtu"), key=frame_key)
    if not files:
        raise FileNotFoundError(f"No VTU files found in {args.vtu_directory}")

    rows: list[dict[str, float]] = []
    for path in files:
        mesh = meshio.read(path)
        row: dict[str, float] = {"gamma": load_from_path(path)}
        ftotal = {
            component: scalar_cell_field(mesh, f"Ftotal{component}")
            for component in COMPONENTS
        }
        elastic = {
            component: scalar_cell_field(mesh, f"F_E{component}")
            for component in COMPONENTS
        }
        plastic_history = {
            component: scalar_cell_field(mesh, f"T{component}")
            for component in COMPONENTS
        }
        # Full current approach: T_full = F_E T = F_E F_P H.
        tfull = {
            "11": elastic["11"] * plastic_history["11"]
            + elastic["12"] * plastic_history["21"],
            "12": elastic["11"] * plastic_history["12"]
            + elastic["12"] * plastic_history["22"],
            "21": elastic["21"] * plastic_history["11"]
            + elastic["22"] * plastic_history["21"],
            "22": elastic["21"] * plastic_history["12"]
            + elastic["22"] * plastic_history["22"],
        }
        for component in COMPONENTS:
            row[f"Ftotal{component}"] = float(np.mean(ftotal[component]))
            row[f"Tfull{component}"] = float(np.mean(tfull[component]))
        rows.append(row)

    gamma = np.asarray([row["gamma"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for ax, component in zip(axes.flat, COMPONENTS):
        ftotal = np.asarray([row[f"Ftotal{component}"] for row in rows])
        tfull_value = np.asarray([row[f"Tfull{component}"] for row in rows])
        ax.plot(
            gamma,
            ftotal,
            color="#1769aa",
            linewidth=2.2,
            marker="o",
            markersize=3,
            label=rf"$\langle \mathcal{{F}}_{{{component}}}\rangle$",
        )
        ax.plot(
            gamma,
            tfull_value,
            color="#d95f02",
            linewidth=2.0,
            linestyle="--",
            marker="s",
            markersize=2.7,
            label=rf"$\langle T^{{\mathrm{{full}}}}_{{{component}}}\rangle$",
        )
        ax.set_title(rf"Component ${component[0]}{component[1]}$")
        ax.set_ylabel("element-wise mean")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)

    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\gamma$")
    fig.suptitle(
        r"Reconstructed full deformation $\mathcal{F}$ vs. $T^{\mathrm{full}}=F_E T$",
        fontsize=15,
    )
    fig.tight_layout()
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    output_csv = args.output_csv or args.output_png.with_suffix(".csv")
    fieldnames = ["gamma"] + [
        f"{matrix_name}{component}_mean"
        for component in COMPONENTS
        for matrix_name in ("Ftotal", "Tfull")
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "gamma": row["gamma"],
                    **{
                        f"{matrix_name}{component}_mean": row[
                            f"{matrix_name}{component}"
                        ]
                        for component in COMPONENTS
                        for matrix_name in ("Ftotal", "Tfull")
                    },
                }
            )

    print(f"Wrote {args.output_png}")
    print(f"Wrote {output_csv}")


if __name__ == "__main__":
    main()
