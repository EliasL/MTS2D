#!/usr/bin/env python3
"""Select and animate the longest element histories in the two shear runs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import meshio
import numpy as np


HERE = Path(__file__).resolve().parent
SIMULATION_SCRIPTS = HERE.parents[2] / "SimulationScripts"
sys.path.insert(0, str(SIMULATION_SCRIPTS))

from MTMath.poincareEnergy import C2PoincareDisk  # noqa: E402
from Plotting.dataFunctions import infer_strain_from_vtu, resolve_vtu_files  # noqa: E402
from Plotting.element_tracking import ElementMatrixHistory  # noqa: E402
import Plotting.element_tracking_animation as tracking_animation  # noqa: E402


RUNS = {
    "current": next((HERE / "current").glob("simpleShear*")),
    "mathcalF": next((HERE / "mathcalF").glob("simpleShear*")),
}
OUTPUT_DIR = HERE / "poincare_longest_paths"
SUMMARY_PATH = OUTPUT_DIR / "longest_path_selection.json"
RANKING_PATH = OUTPUT_DIR / "longest_path_rankings.csv"


def matrix_from_fields(fields: dict[str, np.ndarray], prefix: str) -> np.ndarray:
    return np.stack(
        [
            fields[f"{prefix}11"],
            fields[f"{prefix}12"],
            fields[f"{prefix}21"],
            fields[f"{prefix}22"],
        ],
        axis=1,
    ).reshape(-1, 2, 2)


def read_fields(path: str | Path) -> dict[str, np.ndarray]:
    mesh = meshio.read(path)
    return {
        name: np.asarray(values[0], dtype=float).reshape(-1)
        for name, values in mesh.cell_data.items()
    }


def matrices_for_frame(path: str | Path, run_name: str) -> tuple[np.ndarray, np.ndarray]:
    fields = read_fields(path)
    branch = matrix_from_fields(fields, "T")
    if run_name == "mathcalF":
        total = matrix_from_fields(fields, "Ftotal")
    else:
        elastic = matrix_from_fields(fields, "F_E")
        total = elastic @ branch
    return branch, total


def poincare_coordinates(matrices: np.ndarray) -> np.ndarray:
    metric = np.swapaxes(matrices, -1, -2) @ matrices
    x, y = C2PoincareDisk(metric)
    coordinates = np.column_stack((np.asarray(x).reshape(-1), np.asarray(y).reshape(-1)))
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("Non-finite Poincare coordinates encountered.")
    return coordinates


def hyperbolic_increment(previous: np.ndarray, current: np.ndarray) -> np.ndarray:
    z0 = previous[:, 0] + 1j * previous[:, 1]
    z1 = current[:, 0] + 1j * current[:, 1]
    denominator = 1.0 - np.conjugate(z0) * z1
    relative_radius = np.abs((z1 - z0) / denominator)
    if np.any(relative_radius > 1 + 1e-10):
        bad = int(np.flatnonzero(relative_radius > 1 + 1e-10)[0])
        raise ValueError(
            f"Element {bad} produced an invalid relative Poincare radius "
            f"{relative_radius[bad]}."
        )
    relative_radius = np.clip(relative_radius, 0.0, 1.0 - 1e-14)
    return 2.0 * np.arctanh(relative_radius)


def scan_run(run_name: str, run_dir: Path) -> dict[str, object]:
    paths = tuple(Path(path) for path in resolve_vtu_files(run_dir / "collection.pvd"))
    if len(paths) < 2:
        raise ValueError(f"Not enough VTU frames in {run_dir}")

    branch0, total0 = matrices_for_frame(paths[0], run_name)
    branch_length = np.zeros(len(branch0), dtype=float)
    total_length = np.zeros(len(total0), dtype=float)
    previous_branch = poincare_coordinates(branch0)
    previous_total = poincare_coordinates(total0)

    for frame_index, path in enumerate(paths[1:], start=1):
        branch, total = matrices_for_frame(path, run_name)
        branch_coordinates = poincare_coordinates(branch)
        total_coordinates = poincare_coordinates(total)
        branch_length += hyperbolic_increment(previous_branch, branch_coordinates)
        total_length += hyperbolic_increment(previous_total, total_coordinates)
        previous_branch = branch_coordinates
        previous_total = total_coordinates
        if (frame_index + 1) % 10 == 0 or frame_index + 1 == len(paths):
            print(f"{run_name}: scanned {frame_index + 1}/{len(paths)} frames", flush=True)

    total_order = np.argsort(total_length)[::-1]
    branch_order = np.argsort(branch_length)[::-1]
    return {
        "run": run_name,
        "frame_count": len(paths),
        "element_count": len(total_length),
        "selected_element": int(total_order[0]),
        "selected_total_path_length": float(total_length[total_order[0]]),
        "selected_branch_path_length": float(branch_length[total_order[0]]),
        "longest_branch_element": int(branch_order[0]),
        "longest_branch_path_length": float(branch_length[branch_order[0]]),
        "top_total": [
            {
                "element": int(index),
                "total_path_length": float(total_length[index]),
                "branch_path_length": float(branch_length[index]),
            }
            for index in total_order[:20]
        ],
        "total_path_lengths": total_length,
        "branch_path_lengths": branch_length,
    }


def select_longest_paths() -> dict[str, dict[str, object]]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    raw = {name: scan_run(name, directory) for name, directory in RUNS.items()}
    summary: dict[str, dict[str, object]] = {}
    for name, result in raw.items():
        summary[name] = {
            key: value
            for key, value in result.items()
            if key not in {"total_path_lengths", "branch_path_lengths"}
        }

    with SUMMARY_PATH.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    with RANKING_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["run", "rank", "element", "total_path_length", "branch_path_length"]
        )
        for name, result in raw.items():
            total_lengths = result["total_path_lengths"]
            branch_lengths = result["branch_path_lengths"]
            order = np.argsort(total_lengths)[::-1]
            for rank, element in enumerate(order, start=1):
                writer.writerow(
                    [
                        name,
                        rank,
                        int(element),
                        float(total_lengths[element]),
                        float(branch_lengths[element]),
                    ]
                )
    return summary


def load_one_element(
    run_name: str,
    run_dir: Path,
    element_index: int,
) -> tuple[ElementMatrixHistory, ElementMatrixHistory, np.ndarray]:
    paths = tuple(Path(path) for path in resolve_vtu_files(run_dir / "collection.pvd"))
    branch_matrices = np.empty((len(paths), 2, 2), dtype=float)
    total_matrices = np.empty((len(paths), 2, 2), dtype=float)
    loads = np.empty(len(paths), dtype=float)
    for frame_index, path in enumerate(paths):
        branch, total = matrices_for_frame(path, run_name)
        branch_matrices[frame_index] = branch[element_index]
        total_matrices[frame_index] = total[element_index]
        load = infer_strain_from_vtu(path)
        if load is None or not np.isfinite(load):
            raise ValueError(f"Could not infer load from {path}")
        loads[frame_index] = load
    branch_history = ElementMatrixHistory(paths, branch_matrices, element_index, "T")
    total_name = "Ftotal" if run_name == "mathcalF" else "Tfull"
    total_history = ElementMatrixHistory(paths, total_matrices, element_index, total_name)
    return branch_history, total_history, loads


def render_with_label(
    branch_history: ElementMatrixHistory,
    total_history: ElementMatrixHistory,
    timeline: tracking_animation.GammaTimeline,
    output_path: Path,
    label: str,
    *,
    fps: int,
) -> None:
    original_line = tracking_animation.Line2D

    def labelled_line(*args, **kwargs):
        if kwargs.get("label") == r"$F_E T$":
            kwargs["label"] = label
        return original_line(*args, **kwargs)

    tracking_animation.Line2D = labelled_line
    try:
        tracking_animation.render_poincare_animation(
            branch_history,
            timeline,
            output_path,
            ffmpeg_executable="/opt/homebrew/bin/ffmpeg",
            reconstructed_history=total_history,
            fps=fps,
            dpi=120,
            grid_depth=5,
            grid_samples=50,
        )
    finally:
        tracking_animation.Line2D = original_line


def render_four_videos(
    summary: dict[str, dict[str, object]],
    *,
    gamma_per_frame: float,
    fps: int,
) -> list[Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    selected = {
        origin: int(result["selected_element"]) for origin, result in summary.items()
    }
    outputs: list[Path] = []
    for origin_run, element_index in selected.items():
        for target_run, target_dir in RUNS.items():
            print(
                f"Rendering element {element_index} (selected from {origin_run}) "
                f"in {target_run}",
                flush=True,
            )
            branch_history, total_history, loads = load_one_element(
                target_run, target_dir, element_index
            )
            timeline = tracking_animation.build_gamma_timeline(
                branch_history,
                loads,
                total_history,
                gamma_per_frame=gamma_per_frame,
                camera_smoothing_gamma=0.04,
            )
            output_path = OUTPUT_DIR / (
                f"{target_run}_element_{element_index}_selected_from_{origin_run}.mp4"
            )
            label = r"$\mathcal{F}$" if target_run == "mathcalF" else r"$F_E T$"
            render_with_label(
                branch_history,
                total_history,
                timeline,
                output_path,
                label,
                fps=fps,
            )
            outputs.append(output_path)
    return outputs


def read_summary() -> dict[str, dict[str, object]]:
    if not SUMMARY_PATH.is_file():
        raise FileNotFoundError(
            f"Selection summary not found: {SUMMARY_PATH}. Run --select-only first."
        )
    with SUMMARY_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--select-only", action="store_true")
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--gamma-per-frame", type=float, default=0.004)
    parser.add_argument("--fps", type=int, default=30)
    args = parser.parse_args()
    if args.select_only and args.render_only:
        parser.error("Choose at most one of --select-only and --render-only.")

    if args.render_only:
        summary = read_summary()
    else:
        summary = select_longest_paths()
    if not args.select_only:
        outputs = render_four_videos(
            summary,
            gamma_per_frame=args.gamma_per_frame,
            fps=args.fps,
        )
        print("Rendered:")
        for path in outputs:
            print(path)


if __name__ == "__main__":
    main()
