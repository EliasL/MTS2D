#!/usr/bin/env python3
"""Plot edge-flip J(theta) scan CSVs.

Expected input layout:

  test_data/edge_flip_j_scan/
    scenario_name/
      metadata.csv
      option_a_element_1.csv
      option_a_element_2.csv
      option_b_element_1.csv
      option_b_element_2.csv
    logged_or_single_element_scenario/
      self_element.csv

Outputs:

  test_data/edge_flip_j_scan/all_scenarios_option_a.png
  test_data/edge_flip_j_scan/all_scenarios_option_b.png
  test_data/edge_flip_j_scan/<scenario>/options_ab_elements.png
  test_data/edge_flip_j_scan/<scenario>/j_breakdown_e_components.png
  test_data/edge_flip_j_scan/<scenario>/j_breakdown_e_squared_components.png
  test_data/edge_flip_j_scan/<scenario>/j_breakdown_rotation_energy.png
  test_data/edge_flip_j_scan/<scenario>/j_breakdown_stress_components.png
  test_data/edge_flip_j_scan/<single_element_scenario>/self_element.png
  test_data/edge_flip_j_scan/<scenario>/options_ab_elements_animation.mp4
  test_data/edge_flip_j_scan/<scenario>/option_a_element_1_elastic_reduction_animation.mp4
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache" / "matplotlib")
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Polygon
except ModuleNotFoundError as exc:
    script_path = Path(__file__).resolve()
    venv_python = script_path.parent / ".venv" / "bin" / "python"
    if venv_python.exists() and Path(sys.executable).resolve() != venv_python.resolve():
        os.execv(str(venv_python), [str(venv_python), str(script_path), *sys.argv[1:]])

    raise SystemExit(
        "This script requires matplotlib. If you installed it in tools/.venv, "
        "run `tools/.venv/bin/python tools/plot_edge_flip_j_scans.py` or activate "
        "the environment with `source tools/.venv/bin/activate` before running the "
        "script. To install it, run "
        "`uv pip install -r tools/requirements.txt --python tools/.venv/bin/python`."
    ) from exc

try:
    from tqdm.auto import tqdm
except ModuleNotFoundError:
    tqdm = None


SCENARIO_ORDER = [
    "integer_vertical_shear",
    "integer_horizontal_shear",
    "double_integer_horizontal_shear",
    "half_horizontal_shear",
    "half_pure_shear",
    "one_pure_shear",
    "logged_element_3923",
]
DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "test_data" / "edge_flip_j_scan"
DEFAULT_ANIMATION_SCENARIO = "half_horizontal_shear"
ANIMATION_PANEL_ORDER = (
    "option_a_element_1",
    "option_a_element_2",
    "option_b_element_1",
    "option_b_element_2",
)
CURRENT_NODE_LABELS = ("a", "b", "c")
REFERENCE_NODE_LABELS = ("A", "B", "C")
MATRIX_COMPONENTS = ((0, 0), (0, 1), (1, 0), (1, 1))
CANONICAL_REFERENCE_NODES = np.array(
    [
        [-0.5, -0.5],
        [0.5, -0.5],
        [-0.5, 0.5],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class ScanRow:
    theta: float
    theta_degrees: float
    theta_old: float
    theta_old_degrees: float
    J: float
    energy_old: float
    energy_new: float
    valid: bool
    best: bool
    Q_old: tuple[tuple[float, float], tuple[float, float]]
    Q_new: tuple[tuple[float, float], tuple[float, float]]
    F_old: tuple[tuple[float, float], tuple[float, float]]
    F_new: tuple[tuple[float, float], tuple[float, float]]
    P_new: tuple[tuple[float, float], tuple[float, float]]
    sigma_old: tuple[tuple[float, float], tuple[float, float]]
    sigma_new: tuple[tuple[float, float], tuple[float, float]]
    E_old: tuple[tuple[float, float], tuple[float, float]]
    E_new: tuple[tuple[float, float], tuple[float, float]]
    delta_E: tuple[tuple[float, float], tuple[float, float]]


@dataclass(frozen=True)
class ScanSeries:
    key: str
    label: str
    x: list[float]
    y: list[float]
    best: list[bool]
    rows: list[ScanRow]
    color: str
    linestyle: str
    linewidth: float


def pretty_label(name: str) -> str:
    return name.replace("_", " ")


def log_progress(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def parse_bool_flag(raw: str, field_name: str, path: Path) -> bool:
    if raw == "0":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{path}: expected {field_name} to be 0 or 1, got {raw!r}")


def parse_matrix2d(
    row: dict[str, str],
    prefix: str,
    path: Path,
) -> tuple[tuple[float, float], tuple[float, float]]:
    try:
        return (
            (float(row[f"{prefix}_00"]), float(row[f"{prefix}_01"])),
            (float(row[f"{prefix}_10"]), float(row[f"{prefix}_11"])),
        )
    except KeyError as exc:
        raise ValueError(f"{path}: missing expected column {exc.args[0]!r}") from exc


def load_scan_rows(path: Path) -> list[ScanRow]:
    rows: list[ScanRow] = []

    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            valid = parse_bool_flag(row["valid"], "valid", path)
            rows.append(
                ScanRow(
                    theta=float(row["theta"]),
                    theta_degrees=float(row["theta_degrees"]),
                    theta_old=float(row["theta_old"]),
                    theta_old_degrees=float(row["theta_old_degrees"]),
                    J=float(row["J"]) if valid else math.nan,
                    energy_old=float(row["energy_old"]),
                    energy_new=float(row["energy_new"]),
                    valid=valid,
                    best=parse_bool_flag(row["is_best"], "is_best", path),
                    Q_old=parse_matrix2d(row, "Q_old", path),
                    Q_new=parse_matrix2d(row, "Q_new", path),
                    F_old=parse_matrix2d(row, "F_old", path),
                    F_new=parse_matrix2d(row, "F_new", path),
                    P_new=parse_matrix2d(row, "P_new", path),
                    sigma_old=parse_matrix2d(row, "sigma_old", path),
                    sigma_new=parse_matrix2d(row, "sigma_new", path),
                    E_old=parse_matrix2d(row, "E_old", path),
                    E_new=parse_matrix2d(row, "E_new", path),
                    delta_E=parse_matrix2d(row, "E_new_minus_E_old", path),
                )
            )

    if not rows:
        raise ValueError(f"{path}: no scan rows found")
    return rows


def load_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            key = row["key"]
            if key in metadata:
                raise ValueError(f"{path}: duplicate metadata key {key!r}")
            metadata[key] = row["value"]
    if not metadata:
        raise ValueError(f"{path}: no metadata rows found")
    return metadata


def load_series(
    path: Path,
    key: str,
    label: str,
    color: str,
    linestyle: str,
    linewidth: float,
) -> ScanSeries:
    rows = load_scan_rows(path)
    return ScanSeries(
        key=key,
        label=label,
        x=[row.theta_degrees for row in rows],
        y=[row.J for row in rows],
        best=[row.best for row in rows],
        rows=rows,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
    )


def single_scenario_series_map(scenario_dir: Path) -> dict[str, ScanSeries]:
    return {item.key: item for item in collect_single_scenario(scenario_dir)}


def scenario_dirs(root: Path) -> list[Path]:
    candidates = [p for p in root.iterdir() if p.is_dir()]
    order = {name: i for i, name in enumerate(SCENARIO_ORDER)}
    return sorted(candidates, key=lambda p: (order.get(p.name, 10_000), p.name))


def finite_limits(series: list[list[float]]) -> tuple[float, float] | None:
    values = [v for ys in series for v in ys if math.isfinite(v)]
    if not values:
        return None
    ymin = min(values)
    ymax = max(values)
    if ymin == ymax:
        pad = max(1.0, abs(ymin) * 0.05)
    else:
        pad = 0.05 * (ymax - ymin)
    return ymin - pad, ymax + pad


def x_limits(series: list[ScanSeries]) -> tuple[float, float]:
    values = [x for item in series for x in item.x]
    if not values:
        raise ValueError("Cannot determine x limits for an empty series list")
    xmin = min(values)
    xmax = max(values)
    if xmin == xmax:
        pad = 1.0
    else:
        pad = 0.02 * (xmax - xmin)
    return xmin - pad, xmax + pad


def collect_all_scenarios(root: Path, option: str) -> list[ScanSeries]:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    series: list[ScanSeries] = []

    for idx, scenario_dir in enumerate(scenario_dirs(root)):
        color = colors[idx % len(colors)]
        scenario_label = pretty_label(scenario_dir.name)
        for element, linestyle, linewidth in ((1, "-", 2.8), (2, "--", 1.5)):
            key = f"{scenario_dir.name}_option_{option}_element_{element}"
            path = scenario_dir / f"option_{option}_element_{element}.csv"
            if not path.exists():
                continue
            series.append(
                load_series(
                    path,
                    key=key,
                    label=f"{scenario_label} e{element}",
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                )
            )
    return series


def collect_single_scenario(scenario_dir: Path) -> list[ScanSeries]:
    styles = {
        ("a", 1): ("tab:blue", "-", 2.8, "option A e1"),
        ("a", 2): ("tab:blue", "--", 1.5, "option A e2"),
        ("b", 1): ("tab:orange", "-", 2.8, "option B e1"),
        ("b", 2): ("tab:orange", "--", 1.5, "option B e2"),
    }
    series: list[ScanSeries] = []

    for option in ("a", "b"):
        for element in (1, 2):
            key = f"option_{option}_element_{element}"
            path = scenario_dir / f"{key}.csv"
            if not path.exists():
                continue
            color, linestyle, linewidth, label = styles[(option, element)]
            series.append(
                load_series(
                    path,
                    key=key,
                    label=label,
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                )
            )

    self_path = scenario_dir / "self_element.csv"
    if self_path.exists():
        series.append(
            load_series(
                self_path,
                key="self_element",
                label="self element",
                color="tab:green",
                linestyle="-",
                linewidth=2.4,
            )
        )
    return series


def mark_best(ax: plt.Axes, item: ScanSeries) -> None:
    best_x = [
        xi for xi, yi, bi in zip(item.x, item.y, item.best) if bi and math.isfinite(yi)
    ]
    best_y = [yi for yi, bi in zip(item.y, item.best) if bi and math.isfinite(yi)]
    if best_x:
        ax.scatter(best_x, best_y, color=item.color, s=18, marker="o", zorder=4)


def plot_series(
    output: Path,
    title: str,
    series: list[ScanSeries],
    dpi: int,
    figsize: tuple[float, float],
) -> None:
    if not series:
        return

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    draw_order = sorted(series, key=lambda item: item.linestyle == "--")
    for item in draw_order:
        ax.plot(
            item.x,
            item.y,
            color=item.color,
            linestyle=item.linestyle,
            linewidth=item.linewidth,
            zorder=2 if item.linestyle == "--" else 1,
            label=item.label,
        )

    for item in series:
        mark_best(ax, item)

    limits = finite_limits([item.y for item in series])
    if limits is not None:
        ax.set_ylim(*limits)
    ax.set_title(title)
    ax.set_xlabel("theta (degrees)")
    ax.set_ylabel("J")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize="small", ncols=2)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def plot_all_scenarios(root: Path, option: str, dpi: int) -> Path:
    output = root / f"all_scenarios_option_{option}.png"
    plot_series(
        output,
        f"Option {option.upper()}: J(theta) across scenarios",
        collect_all_scenarios(root, option),
        dpi,
        (11, 6.5),
    )
    return output


def plot_single_scenario(scenario_dir: Path, dpi: int) -> Path:
    has_self_only = (scenario_dir / "self_element.csv").exists() and not any(
        (scenario_dir / f"option_{option}_element_{element}.csv").exists()
        for option in ("a", "b")
        for element in (1, 2)
    )
    output = scenario_dir / ("self_element.png" if has_self_only else "options_ab_elements.png")
    title = (
        f"{pretty_label(scenario_dir.name)}: self element"
        if has_self_only
        else f"{pretty_label(scenario_dir.name)}: options A/B"
    )
    plot_series(
        output,
        title,
        collect_single_scenario(scenario_dir),
        dpi,
        (10, 6),
    )
    return output


def add_best_theta_lines(ax: plt.Axes, item: ScanSeries) -> None:
    for theta_degrees, is_best in zip(item.x, item.best):
        if is_best:
            ax.axvline(theta_degrees, color="0.7", linewidth=1.0, alpha=0.45, zorder=0)


def row_matrix_component_values(
    rows: list[ScanRow],
    field_name: str,
    i: int,
    j: int,
) -> list[float]:
    return [getattr(row, field_name)[i][j] for row in rows]


def row_scalar_values(rows: list[ScanRow], field_name: str) -> list[float]:
    return [getattr(row, field_name) for row in rows]


def breakdown_target_series(scenario_dir: Path) -> ScanSeries | None:
    series_by_key = single_scenario_series_map(scenario_dir)
    for key in ("option_a_element_1", "self_element"):
        if key in series_by_key:
            return series_by_key[key]
    return None


def plot_component_pairs(
    output: Path,
    title: str,
    item: ScanSeries,
    left_field: str,
    right_field: str,
    left_label: str,
    right_label: str,
    dpi: int,
) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, constrained_layout=True)

    for ax, (i, j) in zip(axes.flat, MATRIX_COMPONENTS):
        ax.plot(
            item.x,
            row_matrix_component_values(item.rows, left_field, i, j),
            color="tab:blue",
            linewidth=2.0,
            label=left_label,
        )
        ax.plot(
            item.x,
            row_matrix_component_values(item.rows, right_field, i, j),
            color="tab:orange",
            linewidth=1.8,
            label=right_label,
        )
        add_best_theta_lines(ax, item)
        ax.set_title(f"component ({i},{j})")
        ax.grid(True, alpha=0.25)

    axes[0, 0].legend(fontsize="small")
    axes[1, 0].set_xlabel("theta (degrees)")
    axes[1, 1].set_xlabel("theta (degrees)")
    axes[0, 0].set_ylabel("value")
    axes[1, 0].set_ylabel("value")
    fig.suptitle(title)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_squared_delta_e_components(scenario_dir: Path, item: ScanSeries, dpi: int) -> Path:
    output = scenario_dir / "j_breakdown_e_squared_components.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, constrained_layout=True)

    for ax, (i, j) in zip(axes.flat, MATRIX_COMPONENTS):
        delta_values = row_matrix_component_values(item.rows, "delta_E", i, j)
        ax.plot(
            item.x,
            [value * value for value in delta_values],
            color="tab:red",
            linewidth=2.0,
        )
        add_best_theta_lines(ax, item)
        ax.set_title(f"squared difference component ({i},{j})")
        ax.grid(True, alpha=0.25)

    axes[1, 0].set_xlabel("theta (degrees)")
    axes[1, 1].set_xlabel("theta (degrees)")
    axes[0, 0].set_ylabel("squared difference")
    axes[1, 0].set_ylabel("squared difference")
    fig.suptitle(f"{pretty_label(scenario_dir.name)}: {item.label} E difference squared")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_rotation_energy_breakdown(scenario_dir: Path, item: ScanSeries, dpi: int) -> Path:
    output = scenario_dir / "j_breakdown_rotation_energy.png"
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True, constrained_layout=True)

    axes[0].plot(item.x, item.y, color="black", linewidth=2.2, label="J")
    add_best_theta_lines(axes[0], item)
    axes[0].set_ylabel("J")
    axes[0].legend(fontsize="small")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(item.x, item.x, color="tab:blue", linewidth=1.8, label="theta new")
    axes[1].plot(
        item.x,
        row_scalar_values(item.rows, "theta_old_degrees"),
        color="tab:orange",
        linewidth=2.0,
        label="theta old",
    )
    add_best_theta_lines(axes[1], item)
    axes[1].set_ylabel("degrees")
    axes[1].legend(fontsize="small")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(
        item.x,
        row_scalar_values(item.rows, "energy_old"),
        color="tab:blue",
        linewidth=2.0,
        label="energy old",
    )
    axes[2].plot(
        item.x,
        row_scalar_values(item.rows, "energy_new"),
        color="tab:orange",
        linewidth=1.8,
        label="energy new",
    )
    add_best_theta_lines(axes[2], item)
    axes[2].set_ylabel("energy")
    axes[2].set_xlabel("theta (degrees)")
    axes[2].legend(fontsize="small")
    axes[2].grid(True, alpha=0.25)

    fig.suptitle(f"{pretty_label(scenario_dir.name)}: {item.label} J, rotation, energy")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_stress_breakdown(scenario_dir: Path, item: ScanSeries, dpi: int) -> Path:
    output = scenario_dir / "j_breakdown_stress_components.png"
    fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True, constrained_layout=True)

    axes[0].plot(item.x, item.y, color="black", linewidth=2.2, label="J")
    add_best_theta_lines(axes[0], item)
    axes[0].set_ylabel("J")
    axes[0].legend(fontsize="small")
    axes[0].grid(True, alpha=0.25)

    for ax, (i, j) in zip(axes[1:], MATRIX_COMPONENTS):
        ax.plot(
            item.x,
            row_matrix_component_values(item.rows, "sigma_old", i, j),
            color="tab:blue",
            linewidth=2.0,
            label="sigma old",
        )
        ax.plot(
            item.x,
            row_matrix_component_values(item.rows, "sigma_new", i, j),
            color="tab:orange",
            linewidth=1.8,
            label="sigma new",
        )
        add_best_theta_lines(ax, item)
        ax.set_ylabel(f"sigma ({i},{j})")
        ax.grid(True, alpha=0.25)

    axes[1].legend(fontsize="small")
    axes[-1].set_xlabel("theta (degrees)")
    fig.suptitle(f"{pretty_label(scenario_dir.name)}: {item.label} J and stress")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_j_breakdown_scenario(scenario_dir: Path, dpi: int) -> list[Path]:
    item = breakdown_target_series(scenario_dir)
    if item is None:
        return []

    outputs = [
        plot_component_pairs(
            scenario_dir / "j_breakdown_e_components.png",
            f"{pretty_label(scenario_dir.name)}: {item.label} E components",
            item,
            "E_old",
            "E_new",
            "E old",
            "E new",
            dpi,
        ),
        plot_squared_delta_e_components(scenario_dir, item, dpi),
        plot_rotation_energy_breakdown(scenario_dir, item, dpi),
        plot_stress_breakdown(scenario_dir, item, dpi),
    ]
    return outputs


def rotation_matrix(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=float)


def center_nodes(nodes: np.ndarray) -> np.ndarray:
    return nodes - nodes.mean(axis=0, keepdims=True)


def closed_polygon(nodes: np.ndarray) -> np.ndarray:
    return np.vstack((nodes, nodes[0]))


def rotated_reference_nodes(theta: float) -> np.ndarray:
    return center_nodes(CANONICAL_REFERENCE_NODES @ rotation_matrix(theta).T)


def current_element_nodes(row: ScanRow) -> np.ndarray:
    F_new = np.array(row.F_new, dtype=float)
    current_edge_matrix = F_new @ rotation_matrix(row.theta)
    nodes = np.array(
        [
            [0.0, 0.0],
            [current_edge_matrix[0, 0], current_edge_matrix[1, 0]],
            [current_edge_matrix[0, 1], current_edge_matrix[1, 1]],
        ],
        dtype=float,
    )
    return center_nodes(nodes)


def p_new_transformed_reference_nodes(row: ScanRow) -> np.ndarray:
    p_new = np.array(row.P_new, dtype=float)
    return center_nodes(rotated_reference_nodes(row.theta) @ p_new.T)


def geometry_limits(series: list[ScanSeries]) -> tuple[tuple[float, float], tuple[float, float]]:
    min_x = math.inf
    max_x = -math.inf
    min_y = math.inf
    max_y = -math.inf

    for item in series:
        for row in item.rows:
            if not row.valid:
                continue
            for nodes in (current_element_nodes(row), rotated_reference_nodes(row.theta)):
                min_x = min(min_x, float(np.min(nodes[:, 0])))
                max_x = max(max_x, float(np.max(nodes[:, 0])))
                min_y = min(min_y, float(np.min(nodes[:, 1])))
                max_y = max(max_y, float(np.max(nodes[:, 1])))

    if not math.isfinite(min_x) or not math.isfinite(min_y):
        raise ValueError("Could not determine geometry limits for animation")

    span = max(max_x - min_x, max_y - min_y)
    pad = 0.12 * span if span > 0.0 else 0.2
    return (min_x - pad, max_x + pad), (min_y - pad, max_y + pad)


def reduction_geometry_limits(item: ScanSeries) -> tuple[tuple[float, float], tuple[float, float]]:
    min_x = math.inf
    max_x = -math.inf
    min_y = math.inf
    max_y = -math.inf

    for row in item.rows:
        if not row.valid:
            continue
        for nodes in (
            current_element_nodes(row),
            rotated_reference_nodes(row.theta),
            p_new_transformed_reference_nodes(row),
        ):
            min_x = min(min_x, float(np.min(nodes[:, 0])))
            max_x = max(max_x, float(np.max(nodes[:, 0])))
            min_y = min(min_y, float(np.min(nodes[:, 1])))
            max_y = max(max_y, float(np.max(nodes[:, 1])))

    if not math.isfinite(min_x) or not math.isfinite(min_y):
        raise ValueError("Could not determine geometry limits for reduction animation")

    span = max(max_x - min_x, max_y - min_y)
    pad = 0.14 * span if span > 0.0 else 0.2
    return (min_x - pad, max_x + pad), (min_y - pad, max_y + pad)


def validate_animation_series(
    scenario_dir: Path,
    series: list[ScanSeries],
    expected_theta_samples: int,
) -> None:
    if len(series) != len(ANIMATION_PANEL_ORDER):
        raise ValueError(
            f"{scenario_dir}: expected {len(ANIMATION_PANEL_ORDER)} scan series for "
            "animation"
        )

    reference_rows = series[0].rows
    if len(reference_rows) != expected_theta_samples:
        raise ValueError(
            f"{scenario_dir}: metadata expects {expected_theta_samples} theta samples "
            f"but {series[0].key} has {len(reference_rows)} rows"
        )

    for item in series[1:]:
        if len(item.rows) != len(reference_rows):
            raise ValueError(
                f"{scenario_dir}: {item.key} has {len(item.rows)} rows, expected "
                f"{len(reference_rows)}"
            )
        for idx, (lhs, rhs) in enumerate(zip(reference_rows, item.rows)):
            if not math.isclose(lhs.theta, rhs.theta, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(
                    f"{scenario_dir}: theta mismatch at row {idx} between "
                    f"{series[0].key} and {item.key}"
                )


def format_value(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.5g}"


def build_animation_series(scenario_dir: Path) -> list[ScanSeries]:
    series_by_key = single_scenario_series_map(scenario_dir)
    missing = [key for key in ANIMATION_PANEL_ORDER if key not in series_by_key]
    if missing:
        raise ValueError(
            f"{scenario_dir}: missing scan CSVs required for animation: "
            + ", ".join(missing)
        )
    return [series_by_key[key] for key in ANIMATION_PANEL_ORDER]


def build_animation_item(scenario_dir: Path, key: str) -> ScanSeries:
    series_by_key = single_scenario_series_map(scenario_dir)
    if key not in series_by_key:
        raise ValueError(f"{scenario_dir}: missing scan CSV required for animation: {key}")
    return series_by_key[key]


def make_animation_writer(fps: int) -> animation.FFMpegWriter:
    if not animation.FFMpegWriter.isAvailable():
        raise RuntimeError("FFMpegWriter is required to save the animation as MP4")
    return animation.FFMpegWriter(
        fps=fps,
        codec="libx264",
        extra_args=[
            "-vf",
            "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "-pix_fmt",
            "yuv420p",
        ],
    )


def save_animation_frames(
    fig: plt.Figure,
    output: Path,
    writer: animation.FFMpegWriter,
    dpi: int,
    frame_count: int,
    update_frame: Callable[[int], list[object]],
    description: str,
) -> None:
    if frame_count <= 0:
        raise ValueError(f"{description}: expected at least one animation frame")

    progress = None
    if tqdm is not None:
        progress = tqdm(total=frame_count, desc=description, unit="frame")
    else:
        log_progress(f"{description}: rendering {frame_count} frames")

    try:
        with writer.saving(fig, str(output), dpi=dpi):
            for frame_number in range(frame_count):
                update_frame(frame_number)
                fig.canvas.draw()
                writer.grab_frame()
                if progress is not None:
                    progress.update(1)
                elif (
                    frame_number == 0
                    or frame_number + 1 == frame_count
                    or (frame_number + 1) % 25 == 0
                ):
                    log_progress(
                        f"{description}: frame {frame_number + 1}/{frame_count}"
                    )
    finally:
        if progress is not None:
            progress.close()


def plot_scenario_animation(
    scenario_dir: Path,
    dpi: int,
    fps: int,
    frame_step: int,
) -> Path:
    if fps <= 0:
        raise ValueError(f"Animation fps must be positive, got {fps}")
    if frame_step <= 0:
        raise ValueError(f"Animation frame step must be positive, got {frame_step}")

    metadata = load_metadata(scenario_dir / "metadata.csv")
    try:
        expected_theta_samples = int(metadata["theta_samples"])
    except KeyError as exc:
        raise ValueError(f"{scenario_dir}: metadata.csv is missing theta_samples") from exc

    series = build_animation_series(scenario_dir)
    validate_animation_series(scenario_dir, series, expected_theta_samples)

    frame_indices = list(range(0, len(series[0].rows), frame_step))
    if frame_indices[-1] != len(series[0].rows) - 1:
        frame_indices.append(len(series[0].rows) - 1)

    writer = make_animation_writer(fps)
    output = scenario_dir / "options_ab_elements_animation.mp4"
    fig, axes = plt.subplot_mosaic(
        [
            ["plot", "plot", "option_a_element_1", "option_a_element_2"],
            ["plot", "plot", "option_b_element_1", "option_b_element_2"],
        ],
        figsize=(15.5, 8.0),
        constrained_layout=True,
    )

    plot_ax = axes["plot"]
    draw_order = sorted(series, key=lambda item: item.linestyle == "--")
    for item in draw_order:
        plot_ax.plot(
            item.x,
            item.y,
            color=item.color,
            linestyle=item.linestyle,
            linewidth=item.linewidth,
            zorder=2 if item.linestyle == "--" else 1,
            label=item.label,
        )
    for item in series:
        mark_best(plot_ax, item)

    y_limits = finite_limits([item.y for item in series])
    if y_limits is not None:
        plot_ax.set_ylim(*y_limits)
    plot_ax.set_xlim(*x_limits(series))
    plot_ax.set_title(f"{pretty_label(scenario_dir.name)}: options A/B")
    plot_ax.set_xlabel("theta (degrees)")
    plot_ax.set_ylabel("J")
    plot_ax.grid(True, alpha=0.25)
    plot_ax.legend(fontsize="small", loc="upper right")

    theta_indicator = plot_ax.axvline(
        series[0].rows[0].theta_degrees,
        color="0.25",
        linewidth=1.3,
        alpha=0.8,
        zorder=3,
    )
    summary_text = plot_ax.text(
        0.02,
        0.98,
        "",
        transform=plot_ax.transAxes,
        va="top",
        ha="left",
        fontsize="small",
        bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "0.8"},
    )

    (xlim, ylim) = geometry_limits(series)
    x_offset = 0.03 * (xlim[1] - xlim[0])
    y_offset = 0.03 * (ylim[1] - ylim[0])
    marker_artists: dict[str, plt.Line2D] = {}
    current_patches: dict[str, Polygon] = {}
    reference_patches: dict[str, Polygon] = {}
    panel_texts: dict[str, plt.Text] = {}
    current_node_texts: dict[str, list[plt.Text]] = {}
    reference_node_texts: dict[str, list[plt.Text]] = {}

    for item in series:
        geom_ax = axes[item.key]
        geom_ax.set_title(item.label)
        geom_ax.set_xlim(*xlim)
        geom_ax.set_ylim(*ylim)
        geom_ax.set_aspect("equal", adjustable="box")
        geom_ax.grid(True, alpha=0.25)
        geom_ax.axhline(0.0, color="0.85", linewidth=0.8, zorder=0)
        geom_ax.axvline(0.0, color="0.85", linewidth=0.8, zorder=0)

        current_patch = Polygon(
            closed_polygon(current_element_nodes(item.rows[0])),
            closed=True,
            facecolor=item.color,
            edgecolor=item.color,
            alpha=0.18,
            linewidth=2.2,
        )
        reference_patch = Polygon(
            closed_polygon(rotated_reference_nodes(item.rows[0].theta)),
            closed=True,
            fill=False,
            edgecolor=item.color,
            linestyle="--",
            linewidth=2.0,
        )
        geom_ax.add_patch(current_patch)
        geom_ax.add_patch(reference_patch)

        current_labels = [
            geom_ax.text(
                0.0,
                0.0,
                label,
                va="bottom",
                ha="left",
                fontsize="small",
                color=item.color,
                fontweight="bold",
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
            )
            for label in CURRENT_NODE_LABELS
        ]
        reference_labels = [
            geom_ax.text(
                0.0,
                0.0,
                label,
                va="bottom",
                ha="right",
                fontsize="small",
                color=item.color,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
            )
            for label in REFERENCE_NODE_LABELS
        ]

        panel_text = geom_ax.text(
            0.03,
            0.97,
            "",
            transform=geom_ax.transAxes,
            va="top",
            ha="left",
            fontsize="x-small",
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "0.8"},
        )

        if item.key == ANIMATION_PANEL_ORDER[0]:
            geom_ax.text(
                0.03,
                0.03,
                "fill: current\ndashed: rotated reference",
                transform=geom_ax.transAxes,
                va="bottom",
                ha="left",
                fontsize="x-small",
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
            )

        (marker,) = plot_ax.plot(
            [],
            [],
            marker="o",
            markersize=7,
            linestyle="None",
            color=item.color,
            zorder=5,
        )
        marker_artists[item.key] = marker
        current_patches[item.key] = current_patch
        reference_patches[item.key] = reference_patch
        panel_texts[item.key] = panel_text
        current_node_texts[item.key] = current_labels
        reference_node_texts[item.key] = reference_labels

    def update(frame_number: int) -> list[object]:
        row_index = frame_indices[frame_number]
        theta_degrees = series[0].rows[row_index].theta_degrees
        theta_indicator.set_xdata([theta_degrees, theta_degrees])

        summary_lines = [f"theta = {theta_degrees:.1f} degrees"]
        artists: list[object] = [theta_indicator, summary_text]

        for item in series:
            row = item.rows[row_index]
            marker = marker_artists[item.key]
            current_patch = current_patches[item.key]
            reference_patch = reference_patches[item.key]
            panel_text = panel_texts[item.key]
            current_labels = current_node_texts[item.key]
            reference_labels = reference_node_texts[item.key]

            if row.valid and math.isfinite(row.J):
                current_nodes = current_element_nodes(row)
                reference_nodes = rotated_reference_nodes(row.theta)
                marker.set_data([row.theta_degrees], [row.J])
                current_patch.set_xy(closed_polygon(current_nodes))
                reference_patch.set_xy(closed_polygon(reference_nodes))
                current_patch.set_visible(True)
                reference_patch.set_visible(True)
                for label_artist, node in zip(current_labels, current_nodes):
                    label_artist.set_position((node[0] + x_offset, node[1] + y_offset))
                    label_artist.set_visible(True)
                for label_artist, node in zip(reference_labels, reference_nodes):
                    label_artist.set_position((node[0] - x_offset, node[1] + y_offset))
                    label_artist.set_visible(True)
                panel_text.set_text(f"J = {format_value(row.J)}")
                summary_lines.append(f"{item.label}: J = {format_value(row.J)}")
            else:
                marker.set_data([], [])
                current_patch.set_visible(False)
                reference_patch.set_visible(False)
                for label_artist in current_labels:
                    label_artist.set_visible(False)
                for label_artist in reference_labels:
                    label_artist.set_visible(False)
                panel_text.set_text("invalid")
                summary_lines.append(f"{item.label}: invalid")

            artists.extend(
                [
                    marker,
                    current_patch,
                    reference_patch,
                    panel_text,
                    *current_labels,
                    *reference_labels,
                ]
            )

        summary_text.set_text("\n".join(summary_lines))
        return artists

    save_animation_frames(
        fig=fig,
        output=output,
        writer=writer,
        dpi=dpi,
        frame_count=len(frame_indices),
        update_frame=update,
        description=f"{scenario_dir.name}: options A/B animation",
    )
    plt.close(fig)
    return output


def plot_option_a_elastic_reduction_animation(
    scenario_dir: Path,
    dpi: int,
    fps: int,
    frame_step: int,
) -> Path:
    if fps <= 0:
        raise ValueError(f"Animation fps must be positive, got {fps}")
    if frame_step <= 0:
        raise ValueError(f"Animation frame step must be positive, got {frame_step}")

    metadata = load_metadata(scenario_dir / "metadata.csv")
    try:
        expected_theta_samples = int(metadata["theta_samples"])
    except KeyError as exc:
        raise ValueError(f"{scenario_dir}: metadata.csv is missing theta_samples") from exc

    item = build_animation_item(scenario_dir, "option_a_element_1")
    if len(item.rows) != expected_theta_samples:
        raise ValueError(
            f"{scenario_dir}: metadata expects {expected_theta_samples} theta samples "
            f"but {item.key} has {len(item.rows)} rows"
        )

    frame_indices = list(range(0, len(item.rows), frame_step))
    if frame_indices[-1] != len(item.rows) - 1:
        frame_indices.append(len(item.rows) - 1)

    writer = make_animation_writer(fps)
    output = scenario_dir / "option_a_element_1_elastic_reduction_animation.mp4"
    fig, axes = plt.subplot_mosaic(
        [["geometry", "components"]],
        figsize=(13.0, 6.5),
        constrained_layout=True,
    )

    geom_ax = axes["geometry"]
    component_ax = axes["components"]
    (xlim, ylim) = reduction_geometry_limits(item)
    x_offset = 0.03 * (xlim[1] - xlim[0])
    y_offset = 0.03 * (ylim[1] - ylim[0])

    geom_ax.set_title(f"{pretty_label(scenario_dir.name)}: option A e1")
    geom_ax.set_xlim(*xlim)
    geom_ax.set_ylim(*ylim)
    geom_ax.set_aspect("equal", adjustable="box")
    geom_ax.grid(True, alpha=0.25)
    geom_ax.axhline(0.0, color="0.85", linewidth=0.8, zorder=0)
    geom_ax.axvline(0.0, color="0.85", linewidth=0.8, zorder=0)

    current_patch = Polygon(
        closed_polygon(current_element_nodes(item.rows[0])),
        closed=True,
        facecolor="tab:blue",
        edgecolor="tab:blue",
        alpha=0.16,
        linewidth=2.2,
        label="current",
    )
    reference_patch = Polygon(
        closed_polygon(rotated_reference_nodes(item.rows[0].theta)),
        closed=True,
        fill=False,
        edgecolor="0.25",
        linestyle="--",
        linewidth=2.0,
        label="rotating reference",
    )
    p_new_patch = Polygon(
        closed_polygon(p_new_transformed_reference_nodes(item.rows[0])),
        closed=True,
        fill=False,
        edgecolor="tab:orange",
        linewidth=2.2,
        label="P new applied to reference",
    )
    geom_ax.add_patch(current_patch)
    geom_ax.add_patch(reference_patch)
    geom_ax.add_patch(p_new_patch)
    geom_ax.legend(fontsize="small", loc="upper right")

    current_labels = [
        geom_ax.text(
            0.0,
            0.0,
            label,
            va="bottom",
            ha="left",
            fontsize="small",
            color="tab:blue",
            fontweight="bold",
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
        )
        for label in CURRENT_NODE_LABELS
    ]
    reference_labels = [
        geom_ax.text(
            0.0,
            0.0,
            label,
            va="bottom",
            ha="right",
            fontsize="small",
            color="0.25",
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
        )
        for label in REFERENCE_NODE_LABELS
    ]
    p_new_labels = [
        geom_ax.text(
            0.0,
            0.0,
            f"{label}'",
            va="top",
            ha="left",
            fontsize="small",
            color="tab:orange",
            fontweight="bold",
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 0.15},
        )
        for label in REFERENCE_NODE_LABELS
    ]
    info_text = geom_ax.text(
        0.03,
        0.97,
        "",
        transform=geom_ax.transAxes,
        va="top",
        ha="left",
        fontsize="x-small",
        bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "0.8"},
    )

    component_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    component_markers: list[plt.Line2D] = []
    for color, (i, j) in zip(component_colors, MATRIX_COMPONENTS):
        component_ax.plot(
            item.x,
            row_matrix_component_values(item.rows, "E_new", i, j),
            color=color,
            linewidth=1.9,
            label=f"E new ({i},{j})",
        )
        (marker,) = component_ax.plot(
            [],
            [],
            marker="o",
            markersize=6,
            linestyle="None",
            color=color,
        )
        component_markers.append(marker)

    add_best_theta_lines(component_ax, item)
    theta_indicator = component_ax.axvline(
        item.rows[0].theta_degrees,
        color="0.25",
        linewidth=1.2,
        alpha=0.8,
        zorder=0,
    )
    component_ax.set_title("E new components")
    component_ax.set_xlabel("theta (degrees)")
    component_ax.set_ylabel("value")
    component_ax.grid(True, alpha=0.25)
    component_ax.legend(fontsize="small", loc="upper right")

    def update(frame_number: int) -> list[object]:
        row_index = frame_indices[frame_number]
        row = item.rows[row_index]
        theta_indicator.set_xdata([row.theta_degrees, row.theta_degrees])
        artists: list[object] = [
            theta_indicator,
            current_patch,
            reference_patch,
            p_new_patch,
            info_text,
            *current_labels,
            *reference_labels,
            *p_new_labels,
            *component_markers,
        ]

        if row.valid:
            current_nodes = current_element_nodes(row)
            reference_nodes = rotated_reference_nodes(row.theta)
            p_new_nodes = p_new_transformed_reference_nodes(row)
            current_patch.set_xy(closed_polygon(current_nodes))
            reference_patch.set_xy(closed_polygon(reference_nodes))
            p_new_patch.set_xy(closed_polygon(p_new_nodes))
            current_patch.set_visible(True)
            reference_patch.set_visible(True)
            p_new_patch.set_visible(True)

            for label_artist, node in zip(current_labels, current_nodes):
                label_artist.set_position((node[0] + x_offset, node[1] + y_offset))
                label_artist.set_visible(True)
            for label_artist, node in zip(reference_labels, reference_nodes):
                label_artist.set_position((node[0] - x_offset, node[1] + y_offset))
                label_artist.set_visible(True)
            for label_artist, node in zip(p_new_labels, p_new_nodes):
                label_artist.set_position((node[0] + x_offset, node[1] - y_offset))
                label_artist.set_visible(True)

            p_new = np.array(row.P_new, dtype=float)
            info_text.set_text(
                "\n".join(
                    [
                        f"theta = {row.theta_degrees:.1f} degrees",
                        f"J = {format_value(row.J)}",
                        "P new =",
                        f"[[{format_value(p_new[0, 0])}, {format_value(p_new[0, 1])}],",
                        f" [{format_value(p_new[1, 0])}, {format_value(p_new[1, 1])}]]",
                    ]
                )
            )

            for marker, (i, j) in zip(component_markers, MATRIX_COMPONENTS):
                marker.set_data([row.theta_degrees], [row.E_new[i][j]])
        else:
            current_patch.set_visible(False)
            reference_patch.set_visible(False)
            p_new_patch.set_visible(False)
            for label_artist in current_labels:
                label_artist.set_visible(False)
            for label_artist in reference_labels:
                label_artist.set_visible(False)
            for label_artist in p_new_labels:
                label_artist.set_visible(False)
            info_text.set_text(f"theta = {row.theta_degrees:.1f} degrees\ninvalid row")
            for marker in component_markers:
                marker.set_data([], [])

        return artists

    save_animation_frames(
        fig=fig,
        output=output,
        writer=writer,
        dpi=dpi,
        frame_count=len(frame_indices),
        update_frame=update,
        description=f"{scenario_dir.name}: option A e1 elastic reduction",
    )
    plt.close(fig)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Directory containing scenario CSV folders.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG output DPI.")
    parser.add_argument(
        "--animation-scenario",
        default=DEFAULT_ANIMATION_SCENARIO,
        help=(
            "Scenario name to animate. Defaults to half_horizontal_shear. "
            "Use --skip-animation to disable animation output."
        ),
    )
    parser.add_argument(
        "--animation-fps",
        type=int,
        default=24,
        help="Frames per second for the animation.",
    )
    parser.add_argument(
        "--animation-frame-step",
        type=int,
        default=3,
        help="Use every Nth theta sample in the animation.",
    )
    parser.add_argument(
        "--animation-dpi",
        type=int,
        default=120,
        help="Animation output DPI.",
    )
    parser.add_argument(
        "--skip-animation",
        action="store_true",
        help="Skip the scenario animation file.",
    )
    args = parser.parse_args()

    root = args.root
    if not root.exists():
        raise SystemExit(f"Input folder does not exist: {root}")

    scenario_paths = scenario_dirs(root)
    log_progress(f"Loading edge-flip scan plots from {root}")
    log_progress("Plotting combined option A/B scenario figures")

    outputs = [
        plot_all_scenarios(root, "a", args.dpi),
        plot_all_scenarios(root, "b", args.dpi),
    ]
    log_progress("Wrote combined option A/B plots")

    for index, scenario_dir in enumerate(scenario_paths, start=1):
        log_progress(
            f"[{index}/{len(scenario_paths)}] Plotting scenario curves for {scenario_dir.name}"
        )
        outputs.append(plot_single_scenario(scenario_dir, args.dpi))

    for index, scenario_dir in enumerate(scenario_paths, start=1):
        log_progress(
            f"[{index}/{len(scenario_paths)}] Plotting J breakdowns for {scenario_dir.name}"
        )
        outputs.extend(plot_j_breakdown_scenario(scenario_dir, args.dpi))

    if not args.skip_animation:
        scenario_dir = root / args.animation_scenario
        if not scenario_dir.exists():
            raise SystemExit(f"Animation scenario folder does not exist: {scenario_dir}")
        log_progress(f"Rendering animations for {scenario_dir.name}")
        outputs.append(
            plot_scenario_animation(
                scenario_dir,
                dpi=args.animation_dpi,
                fps=args.animation_fps,
                frame_step=args.animation_frame_step,
            )
        )
        log_progress(f"Finished options A/B animation for {scenario_dir.name}")
        outputs.append(
            plot_option_a_elastic_reduction_animation(
                scenario_dir,
                dpi=args.animation_dpi,
                fps=args.animation_fps,
                frame_step=args.animation_frame_step,
            )
        )
        log_progress(f"Finished elastic-reduction animation for {scenario_dir.name}")

    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
