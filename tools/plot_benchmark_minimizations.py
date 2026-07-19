#!/usr/bin/env python3
"""Create comparative minimization plots as QA PNGs and vector PDFs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Point:
    x: float
    mean: float
    stdev: float
    samples: int


@dataclass(frozen=True)
class Preparation:
    history: str
    size: int
    threads: int
    elapsed_seconds: float
    function_evaluations: float
    edge_flips: float


LABELS = {
    "small_none": "Small event - no reconnect",
    "small_edge": "Small event - edge-flip history",
    "large_none": "Large event - no reconnect",
    "first_none": "First minimization - no reconnect",
    "first_edge": "First minimization - edge flipping",
}


def configure_matplotlib() -> bool:
    """Use native Matplotlib styling with LaTeX/Computer Modern typography."""
    plt.rcdefaults()
    latex_available = all(shutil.which(program) for program in ("latex", "dvipng"))
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [
                "Computer Modern Roman",
                "CMU Serif",
                "cmr10",
                "DejaVu Serif",
            ],
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "text.usetex": latex_available,
            "savefig.bbox": "tight",
        }
    )
    return latex_available


def read_json(path: Path):
    with path.open(encoding="utf-8") as source:
        return json.load(source)


def read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as source:
        return [json.loads(line) for line in source if line.strip()]


def event_key(row: dict) -> str | None:
    load = float(row["fixture_load"])
    history = str(row["history"])
    if math.isclose(load, 0.15) and history == "noReconnecting":
        return "small_none"
    if math.isclose(load, 0.15) and history == "edgeFlipping":
        return "small_edge"
    if math.isclose(load, 0.7) and history == "noReconnecting":
        return "large_none"
    return None


def load_preparations(paths: Iterable[Path]) -> dict[tuple[str, int], Preparation]:
    preparations: dict[tuple[str, int], Preparation] = {}
    for path in paths:
        for record in read_json(path):
            if record.get("status") not in {"generated", "existing"}:
                continue
            native = record.get("native") or {}
            elapsed = native.get("seconds_per_call", native.get("elapsed_seconds"))
            if elapsed is None:
                continue
            preparation = Preparation(
                history=str(record["history"]),
                size=int(record["rows"]),
                threads=int(record.get("preparation_threads", 1)),
                elapsed_seconds=float(elapsed),
                function_evaluations=float(native.get("function_evaluations", 0.0)),
                edge_flips=float(native.get("edge_flips", 0.0)),
            )
            preparations[(preparation.history, preparation.size)] = preparation
    return preparations


def errorbar(axis, points: list[Point], key: str, *, label: bool = True) -> None:
    if not points:
        return
    axis.errorbar(
        [point.x for point in points],
        [point.mean for point in points],
        yerr=[point.stdev for point in points],
        fmt="o-",
        label=LABELS[key] if label else None,
    )


def first_point(axis, preparation: Preparation, key: str, *, label: bool = True) -> None:
    axis.plot(
        preparation.threads,
        preparation.elapsed_seconds,
        label=LABELS[key] if label else None,
        marker="o",
        linestyle="none",
    )


def finish_axes(axis, threads: list[int], *, log_y: bool = True) -> None:
    axis.set_xscale("log", base=2)
    axis.set_xticks(threads)
    axis.set_xticklabels([str(thread) for thread in threads])
    if log_y:
        axis.set_yscale("log")
    axis.minorticks_off()


def add_figure_legend(figure, axes, *, columns: int = 3, y: float = 0.92) -> None:
    handles: list = []
    labels: list[str] = []
    for axis in axes:
        for handle, label in zip(*axis.get_legend_handles_labels()):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y), ncol=columns)


def save_figure(figure, output_dir: Path, stem: str) -> None:
    figure.savefig(output_dir / f"{stem}.png", dpi=220)
    figure.savefig(output_dir / f"{stem}.pdf")
    plt.close(figure)


def workload_points(rows: list[dict], size: int, key: str) -> list[Point]:
    return [
        Point(
            x=float(row["threads"]),
            mean=float(row["mean_seconds_per_call"]),
            stdev=float(row["stdev_seconds_per_call"]),
            samples=int(row["samples"]),
        )
        for row in sorted(rows, key=lambda item: item["threads"])
        if int(row["rows"]) == size and event_key(row) == key
    ]


def thread_scaling_plot(
    output_dir: Path,
    workload: list[dict],
    preparations: dict[tuple[str, int], Preparation],
    sizes: list[int],
    threads: list[int],
) -> None:
    figure, axes = plt.subplots(1, len(sizes), figsize=(10.8, 5.0), squeeze=False)
    axes = list(axes[0])
    for axis, size in zip(axes, sizes):
        for key in ("small_none", "small_edge", "large_none"):
            errorbar(axis, workload_points(workload, size, key), key)
        for history, key in (
            ("noReconnecting", "first_none"),
            ("edgeFlipping", "first_edge"),
        ):
            preparation = preparations.get((history, size))
            if preparation:
                first_point(axis, preparation, key)
        finish_axes(axis, threads)
        axis.set_title(rf"${size}\times {size}$")
    figure.suptitle("Minimization thread scaling\nCauchy, close core placement", y=0.995)
    add_figure_legend(figure, axes, columns=3, y=0.90)
    figure.supxlabel("OpenMP threads", y=0.065)
    figure.supylabel("Seconds per minimization", x=0.015)
    figure.text(
        0.5,
        0.015,
        r"Replay error bars show sample SD ($n=3$). First minimizations were measured once at one thread.",
        ha="center",
    )
    figure.tight_layout(rect=(0.02, 0.09, 1, 0.82))
    save_figure(figure, output_dir, "minimization-thread-scaling")


def size_scaling_plot(
    output_dir: Path,
    workload: list[dict],
    preparations: dict[tuple[str, int], Preparation],
    sizes: list[int],
) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 5.0))
    for key in ("small_none", "small_edge", "large_none"):
        points = [
            Point(
                x=float(row["rows"] * row["cols"]),
                mean=float(row["mean_seconds_per_call"]),
                stdev=float(row["stdev_seconds_per_call"]),
                samples=int(row["samples"]),
            )
            for row in workload
            if int(row["threads"]) == 1 and event_key(row) == key
        ]
        errorbar(axis, sorted(points, key=lambda point: point.x), key)
    for history, key in (
        ("noReconnecting", "first_none"),
        ("edgeFlipping", "first_edge"),
    ):
        selected = [preparation for preparation in preparations.values() if preparation.history == history]
        axis.plot(
            [preparation.size**2 for preparation in sorted(selected, key=lambda item: item.size)],
            [preparation.elapsed_seconds for preparation in sorted(selected, key=lambda item: item.size)],
            label=LABELS[key],
            marker="o",
        )
    nodes = [size**2 for size in sizes]
    axis.set_xscale("log", base=2)
    axis.set_yscale("log")
    axis.set_xticks(nodes)
    axis.set_xticklabels([rf"${size}\times {size}$" for size in sizes])
    axis.set_xlabel("System size")
    axis.set_ylabel("Seconds per minimization")
    axis.set_title("Minimization system-size scaling at one OpenMP thread")
    axis.minorticks_off()
    axis.legend(ncol=2, loc="upper left")
    figure.text(
        0.5,
        0.015,
        r"Replay error bars show sample SD ($n=3$); first-minimization points are single preparations ($n=1$).",
        ha="center",
    )
    figure.tight_layout(rect=(0, 0.055, 1, 1))
    save_figure(figure, output_dir, "minimization-size-scaling")


def speedup_plot(
    output_dir: Path, workload: list[dict], sizes: list[int], threads: list[int]
) -> None:
    figure, axes = plt.subplots(1, len(sizes), figsize=(10.8, 4.7), squeeze=False)
    axes = list(axes[0])
    for axis, size in zip(axes, sizes):
        maximum = 1.0
        for key in ("small_none", "small_edge", "large_none"):
            raw = workload_points(workload, size, key)
            if not raw:
                continue
            base = next(point for point in raw if point.x == 1)
            points: list[Point] = []
            for point in raw:
                speedup = base.mean / point.mean
                uncertainty = 0.0 if point.x == 1 else speedup * math.sqrt(
                    (base.stdev / base.mean) ** 2 + (point.stdev / point.mean) ** 2
                )
                maximum = max(maximum, speedup + uncertainty)
                points.append(Point(point.x, speedup, uncertainty, point.samples))
            errorbar(axis, points, key)
        finish_axes(axis, threads, log_y=False)
        axis.axhline(1.0, linestyle=":")
        axis.set_ylim(0.8, max(3.6, maximum * 1.12))
        axis.set_title(rf"${size}\times {size}$")
    figure.suptitle("Measured minimization speedup - close core placement", y=0.985)
    add_figure_legend(figure, axes, columns=3, y=0.90)
    figure.supxlabel("OpenMP threads", y=0.065)
    figure.supylabel("Speedup relative to one thread", x=0.015)
    figure.text(
        0.5,
        0.015,
        "Error bars propagate sample SD from the one-thread and multithread timing groups.",
        ha="center",
    )
    figure.tight_layout(rect=(0.02, 0.09, 1, 0.82))
    save_figure(figure, output_dir, "minimization-speedup")


def raw_group_points(
    samples: list[dict], size: int, threads: list[int], key: str, metric: str
) -> list[Point]:
    result: list[Point] = []
    for thread in threads:
        values = [
            float(row[metric])
            for row in samples
            if row.get("experiment") == "workload"
            and row.get("affinity_policy") == "close"
            and int(row["rows"]) == size
            and int(row["threads"]) == thread
            and event_key(row) == key
        ]
        if not values:
            continue
        result.append(
            Point(
                x=float(thread),
                mean=statistics.fmean(values),
                stdev=statistics.stdev(values) if len(values) > 1 else 0.0,
                samples=len(values),
            )
        )
    return result


def metric_plot(
    output_dir: Path,
    samples: list[dict],
    sizes: list[int],
    threads: list[int],
    *,
    metric: str,
    title: str,
    ylabel: str,
    stem: str,
) -> None:
    figure, axes = plt.subplots(1, len(sizes), figsize=(10.8, 4.7), squeeze=False)
    axes = list(axes[0])
    for axis, size in zip(axes, sizes):
        for key in ("small_none", "small_edge", "large_none"):
            errorbar(axis, raw_group_points(samples, size, threads, key, metric), key)
        finish_axes(axis, threads)
        axis.set_title(rf"${size}\times {size}$")
    figure.suptitle(title, y=0.985)
    add_figure_legend(figure, axes, columns=3, y=0.90)
    figure.supxlabel("OpenMP threads", y=0.065)
    figure.supylabel(ylabel, x=0.015)
    figure.text(0.5, 0.015, r"Error bars show sample standard deviation ($n=3$).", ha="center")
    figure.tight_layout(rect=(0.02, 0.09, 1, 0.82))
    save_figure(figure, output_dir, stem)


def affinity_plot(
    output_dir: Path, summaries: list[dict], sizes: list[int], threads: list[int]
) -> None:
    figure, axes = plt.subplots(1, len(sizes), figsize=(10.8, 4.7), squeeze=False)
    axes = list(axes[0])
    for axis, size in zip(axes, sizes):
        maximum = 1.0
        for key in ("small_none", "small_edge", "large_none"):
            points: list[Point] = []
            for thread in threads:
                matches = [
                    row
                    for row in summaries
                    if row["experiment"] == "workload"
                    and int(row["rows"]) == size
                    and int(row["threads"]) == thread
                    and event_key(row) == key
                ]
                close = next((row for row in matches if row["affinity_policy"] == "close"), None)
                spread = next((row for row in matches if row["affinity_policy"] == "spread"), None)
                if not close or not spread:
                    continue
                ratio = float(spread["mean_seconds_per_call"]) / float(close["mean_seconds_per_call"])
                uncertainty = ratio * math.sqrt(
                    (float(spread["stdev_seconds_per_call"]) / float(spread["mean_seconds_per_call"])) ** 2
                    + (float(close["stdev_seconds_per_call"]) / float(close["mean_seconds_per_call"])) ** 2
                )
                maximum = max(maximum, ratio + uncertainty)
                points.append(Point(float(thread), ratio, uncertainty, min(int(close["samples"]), int(spread["samples"]))))
            errorbar(axis, points, key)
        plotted_threads = [thread for thread in threads if thread > 1]
        finish_axes(axis, plotted_threads, log_y=False)
        axis.axhline(1.0, linestyle=":")
        axis.set_ylim(0.85, max(1.6, maximum * 1.08))
        axis.set_title(rf"${size}\times {size}$")
    figure.suptitle("Core-placement effect - values above one favor close placement", y=0.985)
    add_figure_legend(figure, axes, columns=3, y=0.90)
    figure.supxlabel("OpenMP threads", y=0.065)
    figure.supylabel("Spread time / close time", x=0.015)
    figure.text(0.5, 0.015, "Error bars propagate sample SD from close and spread timing groups.", ha="center")
    figure.tight_layout(rect=(0.02, 0.09, 1, 0.82))
    save_figure(figure, output_dir, "affinity-placement-penalty")


def first_minimization_plot(
    output_dir: Path,
    preparations: dict[tuple[str, int], Preparation],
    sizes: list[int],
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))
    histories = (
        ("noReconnecting", "first_none"),
        ("edgeFlipping", "first_edge"),
    )
    for history, key in histories:
        selected = sorted(
            (preparation for preparation in preparations.values() if preparation.history == history),
            key=lambda preparation: preparation.size,
        )
        x = [preparation.size**2 for preparation in selected]
        axes[0].plot(
            x,
            [preparation.elapsed_seconds for preparation in selected],
            marker="o",
            label=LABELS[key],
        )
        axes[1].plot(
            x,
            [preparation.function_evaluations for preparation in selected],
            marker="o",
            label=LABELS[key],
        )
    nodes = [size**2 for size in sizes]
    for axis, ylabel, title in zip(
        axes,
        ("Seconds", "Function evaluations"),
        ("Elapsed minimization time", "Minimizer path length"),
    ):
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks(nodes)
        axis.set_xticklabels([rf"${size}\times {size}$" for size in sizes])
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.minorticks_off()
    figure.suptitle("First minimization after random noise - one OpenMP thread", y=0.985)
    add_figure_legend(figure, list(axes), columns=2, y=0.90)
    figure.supxlabel("System size", y=0.065)
    figure.text(0.5, 0.015, r"Each point is one generated initial condition ($n=1$); no standard-deviation estimate is available yet.", ha="center")
    figure.tight_layout(rect=(0, 0.09, 1, 0.82))
    save_figure(figure, output_dir, "first-minimization-cost")


def write_plot_data(
    path: Path,
    summaries: list[dict],
    preparations: dict[tuple[str, int], Preparation],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target)
        writer.writerow(
            [
                "event",
                "history",
                "load",
                "size",
                "threads",
                "affinity",
                "mean_seconds",
                "stdev_seconds",
                "samples",
                "mean_function_evaluations",
                "stdev_function_evaluations",
            ]
        )
        for row in summaries:
            key = event_key(row)
            if row["experiment"] != "workload" or key is None:
                continue
            writer.writerow(
                [
                    LABELS[key],
                    row["history"],
                    row["fixture_load"],
                    f"{row['rows']}x{row['cols']}",
                    row["threads"],
                    row["affinity_policy"],
                    row["mean_seconds_per_call"],
                    row["stdev_seconds_per_call"],
                    row["samples"],
                    row["mean_function_evaluations"],
                    row["stdev_function_evaluations"],
                ]
            )
        for preparation in sorted(preparations.values(), key=lambda item: (item.history, item.size)):
            key = "first_edge" if preparation.history == "edgeFlipping" else "first_none"
            writer.writerow(
                [
                    LABELS[key],
                    preparation.history,
                    0.15,
                    f"{preparation.size}x{preparation.size}",
                    preparation.threads,
                    "close",
                    preparation.elapsed_seconds,
                    "",
                    1,
                    preparation.function_evaluations,
                    "",
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path)
    parser.add_argument(
        "--first-minimization",
        type=Path,
        action="append",
        default=[],
        help="initial_condition_preparation.json; may be passed more than once",
    )
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    latex_enabled = configure_matplotlib()
    result_dir = args.result_dir.resolve()
    output_dir = (args.output_dir or result_dir / "analysis_plots").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = read_json(result_dir / "summary.json")
    samples = read_jsonl(result_dir / "samples.jsonl")
    workload = [
        row
        for row in summaries
        if row["experiment"] == "workload"
        and row["affinity_policy"] == "close"
        and event_key(row) is not None
    ]
    preparations = load_preparations(path.resolve() for path in args.first_minimization)
    sizes = sorted({int(row["rows"]) for row in workload})
    threads = sorted({int(row["threads"]) for row in workload})
    if not sizes or not threads:
        raise SystemExit("No minimization workload summaries were found.")

    thread_scaling_plot(output_dir, workload, preparations, sizes, threads)
    size_scaling_plot(output_dir, workload, preparations, sizes)
    speedup_plot(output_dir, workload, sizes, threads)
    metric_plot(
        output_dir,
        samples,
        sizes,
        threads,
        metric="function_evaluations",
        title="Minimizer path length",
        ylabel="Function evaluations per minimization",
        stem="minimizer-function-evaluations",
    )
    metric_plot(
        output_dir,
        samples,
        sizes,
        threads,
        metric="seconds_per_function_evaluation",
        title="Cost per minimizer function evaluation",
        ylabel="Seconds per function evaluation",
        stem="time-per-function-evaluation",
    )
    affinity_plot(output_dir, summaries, sizes, threads)
    first_minimization_plot(output_dir, preparations, sizes)
    write_plot_data(output_dir / "plot-data.csv", summaries, preparations)

    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "latex_enabled": latex_enabled,
                "png_plots": sorted(path.name for path in output_dir.glob("*.png")),
                "pdf_plots": sorted(path.name for path in output_dir.glob("*.pdf")),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
