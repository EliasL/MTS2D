#!/usr/bin/env python3
"""Inspect the late path separation of element 4250 in the mathcal-F run."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle
import meshio
import numpy as np


HERE = Path(__file__).resolve().parent
SIMULATION_SCRIPTS = HERE.parents[2] / "SimulationScripts"
sys.path.insert(0, str(SIMULATION_SCRIPTS))

from MTMath.poincareEnergy import C2PoincareDisk  # noqa: E402
from Plotting.dataFunctions import infer_strain_from_vtu, resolve_vtu_files  # noqa: E402
from Plotting.element_tracking_animation import _neighborhood_rings  # noqa: E402


ELEMENT = 4250
RUN_DIR = next((HERE / "mathcalF").glob("simpleShear*"))
OUTPUT_STEM = HERE / "element_4250_final_event"
SELECTED_LOADS = (0.90, 0.91, 0.95, 0.96, 0.97, 0.98)


def cell_fields(mesh: meshio.Mesh) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(values[0], dtype=float).reshape(-1)
        for name, values in mesh.cell_data.items()
    }


def matrix(fields: dict[str, np.ndarray], prefix: str) -> np.ndarray:
    return np.array(
        [
            [fields[f"{prefix}11"][ELEMENT], fields[f"{prefix}12"][ELEMENT]],
            [fields[f"{prefix}21"][ELEMENT], fields[f"{prefix}22"][ELEMENT]],
        ]
    )


def poincare_coordinate(value: np.ndarray) -> np.ndarray:
    x, y = C2PoincareDisk(value.T @ value)
    return np.array([float(x), float(y)])


def frame_data() -> list[dict[str, object]]:
    paths = resolve_vtu_files(RUN_DIR / "collection.pvd")
    frames: list[dict[str, object]] = []
    for path in paths:
        load = float(infer_strain_from_vtu(path))
        if load < 0.85 - 1e-12:
            continue
        mesh = meshio.read(path)
        fields = cell_fields(mesh)
        triangles = np.asarray(mesh.cells_dict["triangle"], dtype=int)
        point_refs = np.asarray(mesh.point_data["refIndex"], dtype=int).reshape(-1)
        target = triangles[ELEMENT]
        plastic_branch = matrix(fields, "T")
        elastic = matrix(fields, "F_E")
        frames.append(
            {
                "load": load,
                "path": Path(path),
                "mesh": mesh,
                "triangles": triangles,
                "point_refs": point_refs,
                "target": target,
                "target_refs": tuple(int(value) for value in point_refs[target]),
                "T": plastic_branch,
                "Ttotal": elastic @ plastic_branch,
                "Ftotal": matrix(fields, "Ftotal"),
                "F": matrix(fields, "F"),
            }
        )
    return frames


def plot_path_and_topology(frames: list[dict[str, object]]) -> None:
    mpl.rcParams.update({"font.size": 9, "savefig.dpi": 300})
    fig = plt.figure(figsize=(15, 7.3), constrained_layout=True)
    grid = fig.add_gridspec(2, 6, height_ratios=(1.3, 1.0))
    ax = fig.add_subplot(grid[0, :])
    ax.add_patch(Circle((0, 0), 1, fill=False, edgecolor="black", linewidth=1.2))
    ax.axhline(0, color="0.88", linewidth=0.8)
    ax.axvline(0, color="0.88", linewidth=0.8)

    loads = np.array([frame["load"] for frame in frames])
    path_specs = (
        ("T", "black", "--", r"plastic branch $T$"),
        ("Ftotal", "#2166ac", "-", r"transported history $\mathcal{F}$"),
        ("F", "#d95f02", ":", r"current local triangle $F$"),
    )
    coordinates: dict[str, np.ndarray] = {}
    for name, color, style, label in path_specs:
        xy = np.stack([poincare_coordinate(frame[name]) for frame in frames])
        coordinates[name] = xy
        ax.plot(xy[:, 0], xy[:, 1], style, color=color, linewidth=2, label=label)
        ax.scatter(xy[:, 0], xy[:, 1], color=color, s=15, zorder=4)

    for load in SELECTED_LOADS:
        index = int(np.flatnonzero(np.isclose(loads, load))[0])
        for name, color, _, _ in path_specs:
            xy = coordinates[name][index]
            ax.scatter(*xy, color=color, edgecolor="white", linewidth=0.7, s=52, zorder=6)
        xy = coordinates["Ftotal"][index]
        offset = (5, 8) if load not in (0.95, 0.96) else (5, -13)
        ax.annotate(
            rf"$\gamma={load:.2f}$",
            xy=xy,
            xytext=offset,
            textcoords="offset points",
            color="#2166ac",
        )

    ax.set_xlim(-0.3, 0.58)
    ax.set_ylim(-0.35, 0.92)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x_p$")
    ax.set_ylabel(r"$y_p$")
    ax.legend(loc="lower left", frameon=False)
    ax.set_title("Element 4250: late Poincaré paths and the corresponding edge flips")

    selected_frames = [
        frames[int(np.flatnonzero(np.isclose(loads, load))[0])]
        for load in SELECTED_LOADS
    ]
    extents = []
    neighborhood_data = []
    for frame in selected_frames:
        points = np.asarray(frame["mesh"].points[:, :2], dtype=float)
        triangles = frame["triangles"]
        target = frame["target"]
        second_ring, first_ring = _neighborhood_rings(triangles, ELEMENT)
        center = points[target].mean(axis=0)
        outer = [points[triangles[index]] - center for index in sorted(second_ring)]
        neighbors = [points[triangles[index]] - center for index in sorted(first_ring)]
        target_polygon = points[target] - center
        all_points = np.concatenate(outer + neighbors + [target_polygon])
        extents.append(np.max(np.abs(all_points)))
        neighborhood_data.append((outer, neighbors, target_polygon, center))
    half_width = 1.07 * max(extents)

    for column, (frame, neighborhood) in enumerate(zip(selected_frames, neighborhood_data)):
        local_ax = fig.add_subplot(grid[1, column])
        outer, neighbors, target_polygon, center = neighborhood
        local_ax.add_collection(
            PolyCollection(outer, facecolors="0.93", edgecolors="0.65", linewidths=0.55)
        )
        local_ax.add_collection(
            PolyCollection(
                neighbors,
                facecolors="#d9e8f5",
                edgecolors="0.35",
                linewidths=0.7,
            )
        )
        local_ax.add_collection(
            PolyCollection(
                [target_polygon],
                facecolors="#d73027",
                edgecolors="black",
                linewidths=1.2,
            )
        )
        mesh = frame["mesh"]
        points = np.asarray(mesh.points[:, :2], dtype=float)
        point_refs = frame["point_refs"]
        triangles = frame["triangles"]
        second_ring, first_ring = _neighborhood_rings(triangles, ELEMENT)
        visible_cells = sorted(second_ring | first_ring | {ELEMENT})
        visible_nodes = np.unique(triangles[visible_cells])
        local_points = points[visible_nodes] - center
        local_ax.scatter(local_points[:, 0], local_points[:, 1], s=7, color="black", zorder=5)
        target_nodes = set(int(index) for index in frame["target"])
        for node_index, local_point in zip(visible_nodes, local_points):
            if int(node_index) not in target_nodes:
                continue
            local_ax.annotate(
                str(int(point_refs[node_index])),
                xy=local_point,
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=8,
                fontweight="bold",
                zorder=7,
            )
        local_ax.set_xlim(-half_width, half_width)
        local_ax.set_ylim(-half_width, half_width)
        local_ax.set_aspect("equal")
        local_ax.set_axis_off()
        refs = ",".join(str(value) for value in frame["target_refs"])
        local_ax.set_title(rf"$\gamma={frame['load']:.2f}$" + "\n" + f"nodes [{refs}]")

    fig.savefig(OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_poincare_topology.png"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_poincare_topology.pdf"), bbox_inches="tight")
    plt.close(fig)


def node_fractional_positions(frame: dict[str, object], node_ids: tuple[int, ...]) -> np.ndarray:
    mesh = frame["mesh"]
    points = np.asarray(mesh.points[:, :2], dtype=float)
    refs = frame["point_refs"]
    origin = points[np.flatnonzero(refs == 0)[0]]
    load = float(frame["load"])
    box = np.array([[70.0, load * 70.0], [0.0, 70.0]])
    fractional = (points - origin) @ np.linalg.inv(box).T
    result = np.empty((len(node_ids), 2), dtype=float)
    for index, node_id in enumerate(node_ids):
        candidates = fractional[refs == node_id]
        if len(candidates) == 0:
            raise ValueError(f"Node {node_id} is absent at gamma={load}")
        result[index] = candidates[0] - np.floor(candidates[0])
    return result


def unwrap_periodic(values: np.ndarray) -> np.ndarray:
    result = values.copy()
    for frame_index in range(1, len(result)):
        delta = result[frame_index] - result[frame_index - 1]
        result[frame_index] -= np.round(delta)
    return result


def plot_node_trajectories(frames: list[dict[str, object]]) -> None:
    node_ids = tuple(
        sorted({node for frame in frames for node in frame["target_refs"]})
    )
    loads = np.array([frame["load"] for frame in frames])
    positions = np.stack([node_fractional_positions(frame, node_ids) for frame in frames])
    positions = unwrap_periodic(positions) * 70.0
    displacements = positions - positions[0]

    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.2), sharex=True, constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0, 1, len(node_ids)))
    for node_index, (node_id, color) in enumerate(zip(node_ids, colors)):
        axes[0].plot(loads, displacements[:, node_index, 0], "-o", ms=3, color=color, label=str(node_id))
        axes[1].plot(loads, displacements[:, node_index, 1], "-o", ms=3, color=color)
    for ax in axes:
        for event_load in (0.91, 0.95, 0.96, 0.98):
            ax.axvline(event_load, color="0.7", linewidth=0.8, linestyle="--")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"de-sheared node displacement $\Delta q_x$")
    axes[1].set_ylabel(r"de-sheared node displacement $\Delta q_y$")
    axes[1].set_xlabel(r"imposed shear $\gamma$")
    axes[0].legend(title="reference node", ncol=len(node_ids), frameon=False, loc="upper left")
    axes[0].set_title("Physical nodes involved in the late element-4250 reconnections")
    fig.savefig(OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_node_trajectories.png"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_node_trajectories.pdf"), bbox_inches="tight")
    plt.close(fig)

    csv_path = OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_data.csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = [
            "gamma",
            "target_reference_nodes",
            "T_poincare_x",
            "T_poincare_y",
            "Ttotal_poincare_x",
            "Ttotal_poincare_y",
            "Ftotal_poincare_x",
            "Ftotal_poincare_y",
            "F_poincare_x",
            "F_poincare_y",
        ]
        for node_id in node_ids:
            header.extend([f"node_{node_id}_dqx", f"node_{node_id}_dqy"])
        writer.writerow(header)
        for frame_index, frame in enumerate(frames):
            row = [
                frame["load"],
                " ".join(str(value) for value in frame["target_refs"]),
                *poincare_coordinate(frame["T"]),
                *poincare_coordinate(frame["Ttotal"]),
                *poincare_coordinate(frame["Ftotal"]),
                *poincare_coordinate(frame["F"]),
            ]
            row.extend(displacements[frame_index].reshape(-1))
            writer.writerow(row)


def plot_fixed_node_triangle_test(frames: list[dict[str, object]]) -> None:
    """Track the same three physical nodes through the final edge flip."""
    node_ids = (2125, 2058, 2194)
    loads = np.array([frame["load"] for frame in frames])
    fractional = np.stack(
        [node_fractional_positions(frame, node_ids) for frame in frames]
    )
    fractional = unwrap_periodic(fractional)

    edge_matrices = []
    for load, node_positions in zip(loads, fractional):
        box = np.array([[70.0, load * 70.0], [0.0, 70.0]])
        physical = node_positions @ box.T
        edge_matrices.append(
            np.column_stack(
                (physical[1] - physical[0], physical[2] - physical[0])
            )
        )
    edge_matrices = np.asarray(edge_matrices)
    reference_index = int(np.flatnonzero(np.isclose(loads, 0.96))[0])
    reference_edges = edge_matrices[reference_index]
    increments = edge_matrices @ np.linalg.inv(reference_edges)
    coordinates = np.stack([poincare_coordinate(value) for value in increments])

    shown_loads = (0.96, 0.97, 0.98, 1.00)
    colors = ("#4daf4a", "#377eb8", "#e41a1c", "#984ea3")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), constrained_layout=True)

    reference_triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    for load, color in zip(shown_loads, colors):
        index = int(np.flatnonzero(np.isclose(loads, load))[0])
        transformed = reference_triangle @ increments[index].T
        closed = np.vstack((transformed, transformed[0]))
        axes[0].plot(
            closed[:, 0],
            closed[:, 1],
            "-o",
            color=color,
            linewidth=2,
            markersize=4,
            label=rf"$\gamma={load:.2f}$",
        )
    axes[0].set_aspect("equal")
    axes[0].axhline(0, color="0.85", linewidth=0.8)
    axes[0].axvline(0, color="0.85", linewidth=0.8)
    axes[0].set_xlabel("relative edge coordinate 1")
    axes[0].set_ylabel("relative edge coordinate 2")
    axes[0].set_title("Same physical nodes: 2125, 2058, 2194")
    axes[0].legend(frameon=False)

    late = loads >= 0.96 - 1e-12
    axes[1].plot(
        coordinates[late, 0],
        coordinates[late, 1],
        color="0.25",
        linewidth=1.5,
    )
    for load, color in zip(shown_loads, colors):
        index = int(np.flatnonzero(np.isclose(loads, load))[0])
        axes[1].scatter(*coordinates[index], color=color, s=58, zorder=4)
        axes[1].annotate(
            rf"$\gamma={load:.2f}$",
            xy=coordinates[index],
            xytext=(5, 6),
            textcoords="offset points",
            color=color,
        )
    axes[1].axhline(0, color="0.85", linewidth=0.8)
    axes[1].axvline(0, color="0.85", linewidth=0.8)
    axes[1].set_aspect("equal")
    axes[1].set_xlabel(r"incremental $x_p$ from $\gamma=0.96$")
    axes[1].set_ylabel(r"incremental $y_p$ from $\gamma=0.96$")
    axes[1].set_title("Fixed-node deformation continues after the final flip")

    fig.suptitle(
        "Material-node test across the element-4250 topology change at "
        r"$\gamma=0.98$"
    )
    fig.savefig(
        OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_fixed_node_triangle.png"),
        bbox_inches="tight",
    )
    fig.savefig(
        OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_fixed_node_triangle.pdf"),
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_total_history_counterexample(frames: list[dict[str, object]]) -> None:
    """Contrast the two total histories with the motion of fixed material nodes."""
    loads = np.array([frame["load"] for frame in frames])
    late = loads >= 0.90 - 1e-12
    history_coordinates = {
        name: np.stack([poincare_coordinate(frame[name]) for frame in frames])
        for name in ("Ttotal", "Ftotal")
    }

    node_ids = (2125, 2058, 2194)
    fractional = np.stack(
        [node_fractional_positions(frame, node_ids) for frame in frames]
    )
    fractional = unwrap_periodic(fractional)
    edge_matrices = []
    for load, node_positions in zip(loads, fractional):
        box = np.array([[70.0, load * 70.0], [0.0, 70.0]])
        physical = node_positions @ box.T
        edge_matrices.append(
            np.column_stack(
                (physical[1] - physical[0], physical[2] - physical[0])
            )
        )
    edge_matrices = np.asarray(edge_matrices)
    reference_index = int(np.flatnonzero(np.isclose(loads, 0.96))[0])
    increments = edge_matrices @ np.linalg.inv(edge_matrices[reference_index])
    node_coordinates = np.stack([poincare_coordinate(value) for value in increments])

    fig, axes = plt.subplots(1, 2, figsize=(9.1, 3.45), constrained_layout=True)
    history_styles = (
        ("Ftotal", "#1769aa", "-", r"transported $\mathcal{F}$"),
        ("Ttotal", "#d95f02", "--", r"reconstructed $T_{\mathrm{total}}$"),
    )
    for name, color, linestyle, label in history_styles:
        xy = history_coordinates[name]
        axes[0].plot(
            xy[late, 0], xy[late, 1], linestyle, color=color,
            linewidth=2.2, marker="o", markersize=3.6, label=label,
        )
    for load in (0.96, 0.97, 0.98):
        index = int(np.flatnonzero(np.isclose(loads, load))[0])
        for name, color, _, _ in history_styles:
            xy = history_coordinates[name][index]
            axes[0].scatter(*xy, color=color, edgecolor="white", linewidth=0.7, s=48, zorder=4)
        axes[0].annotate(
            rf"$\gamma={load:.2f}$",
            xy=history_coordinates["Ftotal"][index],
            xytext=(5, 4), textcoords="offset points", color="#1769aa", fontsize=8,
        )
    axes[0].scatter(0, 0, color="black", marker="+", s=45, zorder=5)
    axes[0].set_xlim(0.12, 0.53)
    axes[0].set_ylim(0.32, 0.86)
    axes[0].set_aspect("equal")
    axes[0].grid(alpha=0.18)
    axes[0].set_xlabel(r"Poincaré $x_p$")
    axes[0].set_ylabel(r"Poincaré $y_p$")
    axes[0].set_title("Stored total histories")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")

    fixed_late = loads >= 0.96 - 1e-12
    axes[1].plot(
        node_coordinates[fixed_late, 0], node_coordinates[fixed_late, 1],
        color="0.2", linewidth=2, marker="o", markersize=4,
    )
    node_colors = {0.96: "#4daf4a", 0.97: "#377eb8", 0.98: "#e41a1c", 1.00: "#984ea3"}
    for load, color in node_colors.items():
        index = int(np.flatnonzero(np.isclose(loads, load))[0])
        axes[1].scatter(*node_coordinates[index], color=color, s=52, zorder=4)
        axes[1].annotate(
            rf"$\gamma={load:.2f}$", xy=node_coordinates[index],
            xytext=(5, 5), textcoords="offset points", color=color, fontsize=8,
        )
    axes[1].axhline(0, color="0.82", linewidth=0.8)
    axes[1].axvline(0, color="0.82", linewidth=0.8)
    axes[1].set_aspect("equal")
    axes[1].set_xlabel(r"incremental $x_p$ from $\gamma=0.96$")
    axes[1].set_ylabel(r"incremental $y_p$ from $\gamma=0.96$")
    axes[1].set_title("Same physical nodes [2125, 2058, 2194]")

    fig.suptitle(
        r"Element 4250: $T_{\mathrm{total}}$ retreats while material deformation continues",
        fontsize=12,
    )
    stem = OUTPUT_STEM.with_name(OUTPUT_STEM.name + "_Ttotal_counterexample")
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    data = frame_data()
    plot_path_and_topology(data)
    plot_node_trajectories(data)
    plot_fixed_node_triangle_test(data)
    plot_total_history_counterexample(data)
