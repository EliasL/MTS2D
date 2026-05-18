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
  test_data/edge_flip_j_scan/<single_element_scenario>/self_element.png
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache" / "matplotlib")
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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


@dataclass(frozen=True)
class ScanSeries:
    label: str
    x: list[float]
    y: list[float]
    best: list[bool]
    color: str
    linestyle: str
    linewidth: float


def pretty_label(name: str) -> str:
    return name.replace("_", " ")


def load_scan(path: Path) -> tuple[list[float], list[float], list[bool]]:
    theta_degrees: list[float] = []
    j_values: list[float] = []
    is_best: list[bool] = []

    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            theta_degrees.append(float(row["theta_degrees"]))
            valid = row["valid"] == "1"
            j_values.append(float(row["J"]) if valid else math.nan)
            is_best.append(row["is_best"] == "1")

    return theta_degrees, j_values, is_best


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


def collect_all_scenarios(root: Path, option: str) -> list[ScanSeries]:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    series: list[ScanSeries] = []

    for idx, scenario_dir in enumerate(scenario_dirs(root)):
        color = colors[idx % len(colors)]
        scenario_label = pretty_label(scenario_dir.name)
        for element, linestyle, linewidth in ((1, "-", 2.8), (2, "--", 1.5)):
            path = scenario_dir / f"option_{option}_element_{element}.csv"
            if not path.exists():
                continue
            x, y, best = load_scan(path)
            series.append(
                ScanSeries(
                    label=f"{scenario_label} e{element}",
                    x=x,
                    y=y,
                    best=best,
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
            path = scenario_dir / f"option_{option}_element_{element}.csv"
            if not path.exists():
                continue
            color, linestyle, linewidth, label = styles[(option, element)]
            x, y, best = load_scan(path)
            series.append(
                ScanSeries(
                    label=label,
                    x=x,
                    y=y,
                    best=best,
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                )
            )

    self_path = scenario_dir / "self_element.csv"
    if self_path.exists():
        x, y, best = load_scan(self_path)
        series.append(
            ScanSeries(
                label="self element",
                x=x,
                y=y,
                best=best,
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Directory containing scenario CSV folders.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG output DPI.")
    args = parser.parse_args()

    root = args.root
    if not root.exists():
        raise SystemExit(f"Input folder does not exist: {root}")

    outputs = [
        plot_all_scenarios(root, "a", args.dpi),
        plot_all_scenarios(root, "b", args.dpi),
    ]
    outputs.extend(plot_single_scenario(p, args.dpi) for p in scenario_dirs(root))

    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
