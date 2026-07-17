#!/usr/bin/env python3
"""Compare the system-spanning avalanches in the two 70x70 shear runs."""

from __future__ import annotations

import csv
import glob
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SIMULATION_SCRIPTS = HERE.parents[2] / "SimulationScripts"
sys.path.insert(0, str(SIMULATION_SCRIPTS))

# Reuse the simulation plotting/data layer so VTU connectivity, periodic-cell
# outlines, and element-energy conventions agree with the existing videos.
from Plotting.dataFunctions import VTUData  # noqa: E402
from Plotting.pyplotFunctions import plot_mesh  # noqa: E402


CURRENT_DIR = next((HERE / "current").glob("simpleShear*"))
MATHCALF_DIR = next((HERE / "mathcalF").glob("simpleShear*"))

EVENTS = {
    "current": {
        "directory": CURRENT_DIR,
        "loads": (0.66, 0.67, 0.68),
        "event_load": 0.67,
        "label": r"Current $T$ trajectory",
    },
    "mathcalF": {
        "directory": MATHCALF_DIR,
        "loads": (0.72, 0.73, 0.74),
        "event_load": 0.73,
        "label": r"Direct $\mathcal{F}$ trajectory",
    },
}


def vtu_at(directory: Path, load: float) -> Path:
    pattern = str(directory / "data" / f"*load={load:g}_*.vtu")
    matches = glob.glob(pattern)
    if len(matches) != 1:
        raise RuntimeError(f"Expected one VTU for gamma={load:g}, found {matches}")
    return Path(matches[0])


def macro_row(directory: Path, load: float) -> pd.Series:
    frame = pd.read_csv(directory / "macroData.csv")
    rows = frame[np.isclose(frame["load"], load)]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one macro row for gamma={load:g}")
    return rows.iloc[0]


def set_cell_view(ax: plt.Axes, gamma: float, size: float = 70.0) -> None:
    pad = 1.5
    ax.set_xlim(-pad, size * (1.0 + gamma) + pad)
    ax.set_ylim(-pad, size + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def overlay_mesh(ax: plt.Axes, data: VTUData) -> None:
    nodes = data.get_nodes()
    triangulation = mtri.Triangulation(
        nodes[:, 0], nodes[:, 1], triangles=data.get_connectivity()
    )
    ax.triplot(
        triangulation,
        color="black",
        linewidth=0.075,
        alpha=0.20,
        rasterized=True,
    )


def make_energy_figure() -> None:
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "figure.dpi": 170,
            "savefig.dpi": 300,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 6.0), constrained_layout=True)
    energy_limits = (0.0, 0.12)
    stage_names = ("Immediately before", "Avalanche step", "Immediately after")

    for row_index, event in enumerate(EVENTS.values()):
        for column_index, (load, stage) in enumerate(zip(event["loads"], stage_names)):
            ax = axes[row_index, column_index]
            path = vtu_at(event["directory"], load)
            data = VTUData(path)
            macro = macro_row(event["directory"], load)

            plot_mesh(
                str(path),
                e_lims=energy_limits,
                mesh_property="energy",
                ax=ax,
                add_colorbar=False,
                add_rombus=True,
            )
            overlay_mesh(ax, data)
            set_cell_view(ax, load)

            ax.set_title(
                stage
                + rf"  $\gamma={load:.2f}$"
                + "\n"
                + rf"$\langle E\rangle={macro['avg_energy']:.4f}$, "
                + rf"$\langle\sigma_{{12}}\rangle={macro['avg_sigma12']:.3f}$, "
                + rf"$N_{{\rm changed}}={int(macro['nr_elements_with_m3_change'])}$"
            )

        axes[row_index, 0].set_ylabel(event["label"], labelpad=12, fontsize=11)

    scalar_map = mpl.cm.ScalarMappable(
        norm=mpl.colors.Normalize(*energy_limits), cmap="coolwarm"
    )
    scalar_map.set_array([])
    colorbar = fig.colorbar(
        scalar_map,
        ax=axes,
        orientation="horizontal",
        fraction=0.055,
        pad=0.035,
        aspect=55,
    )
    colorbar.set_label(r"Element energy $E_i$ (values above 0.12 saturated)")
    fig.suptitle(
        "The same class of system-spanning avalanche occurs in both trajectories",
        fontsize=13,
    )
    for extension in ("png", "pdf"):
        fig.savefig(HERE / f"avalanche_energy_mesh_comparison.{extension}", bbox_inches="tight")
    plt.close(fig)


def cell_field(data: VTUData, name: str) -> np.ndarray:
    return np.asarray(data.get_cell_data(name), dtype=float).ravel()


def t_full_12(data: VTUData) -> np.ndarray:
    return (
        cell_field(data, "F_E11") * cell_field(data, "T12")
        + cell_field(data, "F_E12") * cell_field(data, "T22")
    )


def t_full_21(data: VTUData) -> np.ndarray:
    return (
        cell_field(data, "F_E21") * cell_field(data, "T11")
        + cell_field(data, "F_E22") * cell_field(data, "T21")
    )


def plot_scalar_mesh(
    ax: plt.Axes,
    data: VTUData,
    values: np.ndarray,
    norm: mpl.colors.Normalize,
) -> None:
    nodes = data.get_nodes()
    connectivity = data.get_connectivity()
    triangulation = mtri.Triangulation(
        nodes[:, 0], nodes[:, 1], triangles=connectivity
    )
    ax.tripcolor(
        triangulation,
        facecolors=values,
        cmap="coolwarm",
        norm=norm,
        edgecolors="none",
        rasterized=True,
    )
    ax.triplot(
        triangulation,
        color="black",
        linewidth=0.07,
        alpha=0.16,
        rasterized=True,
    )


def make_history_figure(component: str = "12") -> None:
    if component not in {"12", "21"}:
        raise ValueError(f"Unsupported component: {component}")
    loads = EVENTS["mathcalF"]["loads"]
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 6.0), constrained_layout=True)
    limit = 1.8
    norm = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    imposed_value = (lambda load: load) if component == "12" else (lambda load: 0.0)
    if component == "12":
        row_specs = (
            (r"Direct $\mathcal{F}_{12}-\gamma$", "Ftotal12"),
            (r"Reconstructed $(T_{\mathrm{full}})_{12}-\gamma$", "Tfull12"),
        )
    else:
        row_specs = (
            (r"Direct $\mathcal{F}_{21}$", "Ftotal21"),
            (r"Reconstructed $(T_{\mathrm{full}})_{21}$", "Tfull21"),
        )
    stage_names = ("Immediately before", "Avalanche step", "Immediately after")

    for column_index, (load, stage) in enumerate(zip(loads, stage_names)):
        path = vtu_at(MATHCALF_DIR, load)
        data = VTUData(path)
        fields = {
            "Ftotal12": cell_field(data, "Ftotal12"),
            "Tfull12": t_full_12(data),
            "Ftotal21": cell_field(data, "Ftotal21"),
            "Tfull21": t_full_21(data),
        }
        for row_index, (row_label, field_name) in enumerate(row_specs):
            ax = axes[row_index, column_index]
            field = fields[field_name]
            plot_scalar_mesh(ax, data, field - imposed_value(load), norm)
            set_cell_view(ax, load)
            ax.set_title(
                stage
                + rf"  $\gamma={load:.2f}$"
                + "\n"
                + rf"element mean $={np.mean(field):.3f}$"
            )
            if column_index == 0:
                ax.set_ylabel(row_label, labelpad=12, fontsize=11)

    scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")
    scalar_map.set_array([])
    colorbar = fig.colorbar(
        scalar_map,
        ax=axes,
        orientation="horizontal",
        fraction=0.055,
        pad=0.035,
        aspect=55,
        extend="both",
    )
    if component == "12":
        colorbar.set_label(r"Local shear history relative to imposed shear")
        title = (
            r"On the identical mesh event, $\mathcal{F}$ preserves total shear while "
            r"$T_{\mathrm{full}}$ reclassifies part of it"
        )
        suffix = ""
    else:
        colorbar.set_label(r"Local transverse shear history (imposed value is zero)")
        title = (
            r"The same event creates a large $21$ contribution in $T_{\mathrm{full}}$, "
            r"but not in $\mathcal{F}$"
        )
        suffix = "_21"
    fig.suptitle(title, fontsize=13)
    for extension in ("png", "pdf"):
        fig.savefig(
            HERE / f"avalanche_history_interpretation{suffix}.{extension}",
            bbox_inches="tight",
        )
    plt.close(fig)


def write_event_summary() -> None:
    rows = []
    for method, event in EVENTS.items():
        before_load, event_load, after_load = event["loads"]
        before = macro_row(event["directory"], before_load)
        event_row = macro_row(event["directory"], event_load)
        rows.append(
            {
                "method": method,
                "before_gamma": before_load,
                "event_gamma": event_load,
                "energy_before": before["avg_energy"],
                "energy_after_event": event_row["avg_energy"],
                "total_energy_change_event_step": event_row["total_energy_change"],
                "sigma12_before": before["avg_sigma12"],
                "sigma12_after_event": event_row["avg_sigma12"],
                "sigma12_change_event_step": event_row[
                    "avg_sigma12_change_from_init"
                ],
                "elements_with_m3_change": int(
                    event_row["nr_elements_with_m3_change"]
                ),
                "fraction_of_9800_elements_changed": event_row[
                    "nr_elements_with_m3_change"
                ]
                / 9800.0,
                "participation_fraction": event_row["participationFraction"],
                "m3_participation_fraction": event_row[
                    "m3_participationFraction"
                ],
                "accepted_edge_flips": int(event_row["nr_edge_flips"]),
                "after_gamma": after_load,
            }
        )

    with (HERE / "avalanche_event_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def write_history_difference_summary() -> None:
    rows = []
    for load in EVENTS["mathcalF"]["loads"]:
        data = VTUData(vtu_at(MATHCALF_DIR, load))
        direct = np.stack(
            [
                cell_field(data, "Ftotal11"),
                cell_field(data, "Ftotal12"),
                cell_field(data, "Ftotal21"),
                cell_field(data, "Ftotal22"),
            ],
            axis=1,
        ).reshape(-1, 2, 2)
        elastic = np.stack(
            [
                cell_field(data, "F_E11"),
                cell_field(data, "F_E12"),
                cell_field(data, "F_E21"),
                cell_field(data, "F_E22"),
            ],
            axis=1,
        ).reshape(-1, 2, 2)
        plastic_history = np.stack(
            [
                cell_field(data, "T11"),
                cell_field(data, "T12"),
                cell_field(data, "T21"),
                cell_field(data, "T22"),
            ],
            axis=1,
        ).reshape(-1, 2, 2)
        reconstructed = elastic @ plastic_history
        difference = np.linalg.norm(direct - reconstructed, axis=(1, 2))
        rows.append(
            {
                "gamma": load,
                "mean_direct_F12": np.mean(direct[:, 0, 1]),
                "mean_reconstructed_Tfull12": np.mean(reconstructed[:, 0, 1]),
                "mean_direct_F21": np.mean(direct[:, 1, 0]),
                "mean_reconstructed_Tfull21": np.mean(reconstructed[:, 1, 0]),
                "mean_frobenius_difference": np.mean(difference),
                "q95_frobenius_difference": np.quantile(difference, 0.95),
                "fraction_difference_above_0.5": np.mean(difference > 0.5),
                "fraction_difference_above_1": np.mean(difference > 1.0),
                "fraction_difference_above_2": np.mean(difference > 2.0),
            }
        )

    with (HERE / "avalanche_history_difference_summary.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    make_energy_figure()
    make_history_figure("12")
    make_history_figure("21")
    write_event_summary()
    write_history_difference_summary()
