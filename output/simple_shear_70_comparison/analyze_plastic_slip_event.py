#!/usr/bin/env python3
"""Classify the two large avalanches by integer square-lattice slips."""

from __future__ import annotations

from collections import Counter, deque
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import meshio
import numpy as np


HERE = Path(__file__).resolve().parent
SIMULATION_SCRIPTS = HERE.parents[2] / "SimulationScripts"
sys.path.insert(0, str(SIMULATION_SCRIPTS))

from MTMath.reduction import elastic_reduction  # noqa: E402


@dataclass(frozen=True)
class Event:
    key: str
    label: str
    directory: Path
    before: float
    after: float
    total_field: str


EVENTS = (
    Event(
        key="current",
        label=r"Plastic reconstruction $T_{\mathrm{total}}$",
        directory=next((HERE / "current").glob("simpleShear*")),
        before=0.66,
        after=0.67,
        total_field="Ttotal",
    ),
    Event(
        key="mathcalF",
        label=r"Direct history $\mathcal{F}$",
        directory=next((HERE / "mathcalF").glob("simpleShear*")),
        before=0.72,
        after=0.73,
        total_field="Ftotal",
    ),
)

GRID_MIN = -2
GRID_MAX = 2
SYSTEM_SIZE = 70.0


def horizontal_shear(amount: int) -> np.ndarray:
    return np.array([[1, amount], [0, 1]], dtype=int)


def vertical_shear(amount: int) -> np.ndarray:
    return np.array([[1, 0], [amount, 1]], dtype=int)


GENERATORS = (
    ("H", +1, horizontal_shear(+1)),
    ("H", -1, horizontal_shear(-1)),
    ("V", +1, vertical_shear(+1)),
    ("V", -1, vertical_shear(-1)),
)


def vtu_at(event: Event, load: float) -> Path:
    matches = list((event.directory / "data").glob(f"*load={load:g}_*.vtu"))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one VTU for {event.key} at gamma={load:g}; got {matches}"
        )
    return matches[0]


def read_vtu(event: Event, load: float) -> meshio.Mesh:
    return meshio.read(vtu_at(event, load))


def cell_field(mesh: meshio.Mesh, name: str) -> np.ndarray:
    return np.asarray(mesh.cell_data[name][0], dtype=float).reshape(-1)


def matrix_field(mesh: meshio.Mesh, prefix: str) -> np.ndarray:
    return np.stack(
        [cell_field(mesh, prefix + suffix) for suffix in ("11", "12", "21", "22")],
        axis=1,
    ).reshape(-1, 2, 2)


def total_history(event: Event, mesh: meshio.Mesh) -> np.ndarray:
    if event.total_field == "Ftotal":
        return matrix_field(mesh, "Ftotal")
    if event.total_field == "Ttotal":
        return matrix_field(mesh, "F_E") @ matrix_field(mesh, "T")
    raise ValueError(f"Unknown total field mode: {event.total_field}")


def plastic_branch(total: np.ndarray) -> np.ndarray:
    """Return integer F_P after the square-lattice elastic reduction."""
    metric = np.swapaxes(total, -1, -2) @ total
    _, reduction = elastic_reduction(metric, compute_M=True)
    plastic = np.linalg.inv(reduction)
    rounded = np.rint(plastic).astype(int)
    error = float(np.max(np.abs(plastic - rounded)))
    if error > 1e-9:
        raise RuntimeError(f"Plastic branch is not integer; max error={error:g}")
    determinants = np.rint(np.linalg.det(rounded)).astype(int)
    if not np.all(determinants == 1):
        raise RuntimeError("Plastic reduction returned a non-SL(2,Z) branch")
    return rounded


def matrix_key(matrix: np.ndarray) -> tuple[int, int, int, int]:
    return tuple(int(value) for value in matrix.ravel())


def shortest_shear_words(
    targets: set[tuple[int, int, int, int]],
) -> dict[tuple[int, int, int, int], tuple[tuple[str, int], ...]]:
    """Find a deterministic shortest word in H(+/-1), V(+/-1) for each target."""
    identity = (1, 0, 0, 1)
    pending = set(targets)
    words: dict[tuple[int, int, int, int], tuple[tuple[str, int], ...]] = {}
    queue = deque([(identity, tuple())])
    visited = {identity}
    largest_target_entry = max(max(abs(value) for value in key) for key in targets)
    entry_bound = max(24, 4 * largest_target_entry)
    max_depth = 18

    while queue and pending:
        key, word = queue.popleft()
        if key in pending:
            words[key] = word
            pending.remove(key)
        if len(word) >= max_depth:
            continue
        matrix = np.asarray(key, dtype=int).reshape(2, 2)
        for direction, sign, generator in GENERATORS:
            candidate = matrix @ generator
            candidate_key = matrix_key(candidate)
            if candidate_key in visited:
                continue
            if max(abs(value) for value in candidate_key) > entry_bound:
                continue
            visited.add(candidate_key)
            queue.append((candidate_key, word + ((direction, sign),)))

    if pending:
        raise RuntimeError(f"Could not factor plastic increments: {sorted(pending)}")
    return words


@dataclass
class EventResult:
    event: Event
    after_mesh: meshio.Mesh
    plastic_increment: np.ndarray
    horizontal: np.ndarray
    vertical: np.ndarray
    word_length: np.ndarray
    word_text: list[str]
    two_block_order: list[str]
    heatmap: np.ndarray
    outside_grid: int


def analyze_event(event: Event) -> EventResult:
    before_mesh = read_vtu(event, event.before)
    after_mesh = read_vtu(event, event.after)
    before_plastic = plastic_branch(total_history(event, before_mesh))
    after_plastic = plastic_branch(total_history(event, after_mesh))

    before_inverse = np.rint(np.linalg.inv(before_plastic)).astype(int)
    # Left relative increment: maps the old intermediate configuration to the new.
    increment = after_plastic @ before_inverse
    if not np.all(np.rint(np.linalg.det(increment)).astype(int) == 1):
        raise RuntimeError("Relative plastic increment left SL(2,Z)")

    keys = [matrix_key(matrix) for matrix in increment]
    words = shortest_shear_words(set(keys))
    horizontal = np.empty(len(keys), dtype=int)
    vertical = np.empty(len(keys), dtype=int)
    word_length = np.empty(len(keys), dtype=int)
    word_text: list[str] = []
    two_block_order: list[str] = []

    for index, key in enumerate(keys):
        word = words[key]
        horizontal[index] = sum(sign for direction, sign in word if direction == "H")
        vertical[index] = sum(sign for direction, sign in word if direction == "V")
        word_length[index] = len(word)
        word_text.append(" ".join(f"{direction}{sign:+d}" for direction, sign in word) or "I")

        a, b, c, d = (int(value) for value in increment[index].ravel())
        if d == 1 and a == 1 + b * c:
            two_block_order.append("HV")
        elif a == 1 and d == 1 + b * c:
            two_block_order.append("VH")
        else:
            two_block_order.append("longer_word")

        reconstructed = np.eye(2, dtype=int)
        for direction, sign in word:
            reconstructed = reconstructed @ (
                horizontal_shear(sign) if direction == "H" else vertical_shear(sign)
            )
        if not np.array_equal(reconstructed, increment[index]):
            raise RuntimeError(f"Shear word failed for element {index}")

    heatmap = np.zeros((GRID_MAX - GRID_MIN + 1, GRID_MAX - GRID_MIN + 1), dtype=int)
    inside = (
        (horizontal >= GRID_MIN)
        & (horizontal <= GRID_MAX)
        & (vertical >= GRID_MIN)
        & (vertical <= GRID_MAX)
    )
    for h, v in zip(horizontal[inside], vertical[inside]):
        heatmap[v - GRID_MIN, h - GRID_MIN] += 1

    return EventResult(
        event=event,
        after_mesh=after_mesh,
        plastic_increment=increment,
        horizontal=horizontal,
        vertical=vertical,
        word_length=word_length,
        word_text=word_text,
        two_block_order=two_block_order,
        heatmap=heatmap,
        outside_grid=int(np.count_nonzero(~inside)),
    )


def heatmap_colormap() -> mpl.colors.Colormap:
    cmap = mpl.colormaps["magma"].copy()
    cmap.set_bad("#eeeeee")
    return cmap


def slip_colormap() -> mpl.colors.Colormap:
    return mpl.colors.ListedColormap(
        ["#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"]
    )


def plot_heatmap(
    ax: plt.Axes,
    result: EventResult,
    norm: mpl.colors.Normalize,
    cmap: mpl.colors.Colormap,
) -> mpl.image.AxesImage:
    masked = np.ma.masked_where(result.heatmap == 0, result.heatmap)
    image = ax.imshow(
        masked,
        origin="lower",
        extent=(GRID_MIN - 0.5, GRID_MAX + 0.5, GRID_MIN - 0.5, GRID_MAX + 0.5),
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )
    ax.set_xticks(range(GRID_MIN, GRID_MAX + 1))
    ax.set_yticks(range(GRID_MIN, GRID_MAX + 1))
    ax.set_xlabel(r"net horizontal slips $H$")
    ax.set_ylabel(r"net vertical slips $V$")
    ax.set_aspect("equal")
    threshold = np.sqrt(norm.vmin * norm.vmax)
    for v in range(GRID_MIN, GRID_MAX + 1):
        for h in range(GRID_MIN, GRID_MAX + 1):
            count = int(result.heatmap[v - GRID_MIN, h - GRID_MIN])
            # Magma is dark at small counts and light at large counts.
            color = "white" if 0 < count < threshold else "black"
            ax.text(h, v, str(count), ha="center", va="center", color=color, fontsize=8)
    changed = int(np.count_nonzero((result.horizontal != 0) | (result.vertical != 0)))
    ax.set_title(
        result.event.label
        + "\n"
        + rf"$\gamma={result.event.before:.2f}\rightarrow{result.event.after:.2f}$; "
        + rf"nonzero: {changed:,}; outside: {result.outside_grid}"
    )
    return image


def set_mesh_view(ax: plt.Axes, gamma: float) -> None:
    pad = 1.5
    ax.set_xlim(-pad, SYSTEM_SIZE * (1.0 + gamma) + pad)
    ax.set_ylim(-pad, SYSTEM_SIZE + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    boundary = np.array(
        [
            [0.0, 0.0],
            [SYSTEM_SIZE, 0.0],
            [SYSTEM_SIZE * (1.0 + gamma), SYSTEM_SIZE],
            [SYSTEM_SIZE * gamma, SYSTEM_SIZE],
            [0.0, 0.0],
        ]
    )
    ax.plot(
        boundary[:, 0],
        boundary[:, 1],
        color=mpl.rcParams["text.color"],
        linewidth=0.7,
    )


def plot_mesh_field(
    ax: plt.Axes,
    result: EventResult,
    values: np.ndarray,
    symbol: str,
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
) -> mpl.collections.Collection:
    points = result.after_mesh.points[:, :2]
    triangles = result.after_mesh.cells[0].data
    triangulation = mtri.Triangulation(points[:, 0], points[:, 1], triangles=triangles)
    clipped = np.clip(values, GRID_MIN, GRID_MAX)
    collection = ax.tripcolor(
        triangulation,
        facecolors=clipped,
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        rasterized=True,
    )
    ax.triplot(
        triangulation,
        color=mpl.rcParams["text.color"],
        linewidth=0.055,
        alpha=0.12,
        rasterized=True,
    )
    set_mesh_view(ax, result.event.after)
    outside = int(np.count_nonzero((values < GRID_MIN) | (values > GRID_MAX)))
    suffix = f"; {outside} values saturated" if outside else ""
    ax.set_title(
        result.event.label
        + "\n"
        + rf"net ${symbol}$ field at $\gamma={result.event.after:.2f}$"
        + suffix
    )
    return collection


def save_individual_plots(
    results: tuple[EventResult, ...],
    count_norm: mpl.colors.Normalize,
    count_cmap: mpl.colors.Colormap,
    field_norm: mpl.colors.Normalize,
    field_cmap: mpl.colors.Colormap,
) -> None:
    for result in results:
        fig, ax = plt.subplots(figsize=(6.2, 5.4), constrained_layout=True)
        image = plot_heatmap(ax, result, count_norm, count_cmap)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.047, pad=0.04)
        colorbar.set_label("number of elements (logarithmic color scale)")
        for extension in ("png", "pdf"):
            fig.savefig(
                HERE / f"plastic_slip_heatmap_{result.event.key}.{extension}",
                dpi=240,
                bbox_inches="tight",
            )
        plt.close(fig)

        for values, symbol in ((result.horizontal, "H"), (result.vertical, "V")):
            fig, ax = plt.subplots(figsize=(9.4, 5.2), constrained_layout=True)
            collection = plot_mesh_field(
                ax, result, values, symbol, field_cmap, field_norm
            )
            colorbar = fig.colorbar(
                collection,
                ax=ax,
                ticks=range(GRID_MIN, GRID_MAX + 1),
                fraction=0.035,
                pad=0.025,
                extend="both",
            )
            colorbar.set_label(f"net {symbol} slip count")
            for extension in ("png", "pdf"):
                fig.savefig(
                    HERE
                    / f"plastic_slip_mesh_{result.event.key}_{symbol}.{extension}",
                    dpi=240,
                    bbox_inches="tight",
                )
            plt.close(fig)


def save_composite(
    results: tuple[EventResult, ...],
    count_norm: mpl.colors.Normalize,
    count_cmap: mpl.colors.Colormap,
    field_norm: mpl.colors.Normalize,
    field_cmap: mpl.colors.Colormap,
    suffix: str = "",
    dark: bool = False,
    stem: str = "plastic_slip_event_comparison",
    title: str = "Integer plastic-slip content of the two large avalanches",
) -> None:
    style = "dark_background" if dark else "default"
    with plt.style.context(style):
        fig, axes = plt.subplots(2, 3, figsize=(14.2, 7.4), constrained_layout=True)
        for row, result in enumerate(results):
            plot_heatmap(axes[row, 0], result, count_norm, count_cmap)
            plot_mesh_field(
                axes[row, 1], result, result.horizontal, "H", field_cmap, field_norm
            )
            plot_mesh_field(
                axes[row, 2], result, result.vertical, "V", field_cmap, field_norm
            )
        count_map = mpl.cm.ScalarMappable(norm=count_norm, cmap=count_cmap)
        count_map.set_array([])
        count_bar = fig.colorbar(
            count_map,
            ax=axes[:, 0],
            orientation="horizontal",
            fraction=0.045,
            pad=0.04,
        )
        count_bar.set_label("number of elements (logarithmic color scale)")
        field_map = mpl.cm.ScalarMappable(norm=field_norm, cmap=field_cmap)
        field_map.set_array([])
        field_bar = fig.colorbar(
            field_map,
            ax=axes[:, 1:],
            ticks=range(GRID_MIN, GRID_MAX + 1),
            orientation="horizontal",
            fraction=0.045,
            pad=0.04,
            extend="both",
        )
        field_bar.set_label("net integer slip count (values outside [-2,2] saturated)")
        fig.suptitle(title, fontsize=14)
        output = HERE / f"{stem}{suffix}.png"
        fig.savefig(output, dpi=165, bbox_inches="tight", facecolor=fig.get_facecolor())
        if not dark:
            fig.savefig(HERE / f"{stem}.pdf", bbox_inches="tight")
        plt.close(fig)


def write_tables(results: tuple[EventResult, ...]) -> None:
    with (HERE / "plastic_slip_fields.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "method",
                "before_gamma",
                "after_gamma",
                "element",
                "H",
                "V",
                "shortest_word_length",
                "shortest_word",
                "two_block_order",
                "dP11",
                "dP12",
                "dP21",
                "dP22",
            ]
        )
        for result in results:
            for index, (h, v, length, word, order, matrix) in enumerate(
                zip(
                    result.horizontal,
                    result.vertical,
                    result.word_length,
                    result.word_text,
                    result.two_block_order,
                    result.plastic_increment,
                )
            ):
                writer.writerow(
                    [
                        result.event.key,
                        result.event.before,
                        result.event.after,
                        index,
                        int(h),
                        int(v),
                        int(length),
                        word,
                        order,
                        *[int(value) for value in matrix.ravel()],
                    ]
                )

    with (HERE / "plastic_slip_heatmap_counts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "H", "V", "element_count"])
        for result in results:
            counter = Counter(zip(result.horizontal, result.vertical))
            for (h, v), count in sorted(counter.items()):
                writer.writerow([result.event.key, int(h), int(v), int(count)])

    with (HERE / "plastic_slip_event_classification_summary.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "method",
                "before_gamma",
                "after_gamma",
                "elements",
                "nonzero_net_slip_elements",
                "fraction_nonzero_net_slip",
                "elements_outside_minus2_plus2",
                "elements_not_factorable_as_HV_or_VH",
                "elements_with_shortest_word_longer_than_2",
                "net_H_sum",
                "net_V_sum",
                "absolute_net_slip_units",
            ]
        )
        for result in results:
            changed = (result.horizontal != 0) | (result.vertical != 0)
            writer.writerow(
                [
                    result.event.key,
                    result.event.before,
                    result.event.after,
                    len(changed),
                    int(np.count_nonzero(changed)),
                    float(np.mean(changed)),
                    result.outside_grid,
                    int(sum(order == "longer_word" for order in result.two_block_order)),
                    int(np.count_nonzero(result.word_length > 2)),
                    int(np.sum(result.horizontal)),
                    int(np.sum(result.vertical)),
                    int(np.sum(np.abs(result.horizontal) + np.abs(result.vertical))),
                ]
            )


def main() -> None:
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "figure.dpi": 160,
            "savefig.dpi": 240,
        }
    )
    results = tuple(analyze_event(event) for event in EVENTS)
    largest_count = max(int(np.max(result.heatmap)) for result in results)
    count_norm = mpl.colors.LogNorm(vmin=1, vmax=largest_count)
    count_cmap = heatmap_colormap()
    field_cmap = slip_colormap()
    field_norm = mpl.colors.BoundaryNorm(
        np.arange(GRID_MIN - 0.5, GRID_MAX + 1.5, 1.0), field_cmap.N
    )

    save_individual_plots(results, count_norm, count_cmap, field_norm, field_cmap)
    save_composite(results, count_norm, count_cmap, field_norm, field_cmap)
    save_composite(
        results,
        count_norm,
        count_cmap,
        field_norm,
        field_cmap,
        suffix="_dark",
        dark=True,
    )
    write_tables(results)

    for result in results:
        changed = int(
            np.count_nonzero((result.horizontal != 0) | (result.vertical != 0))
        )
        print(
            f"{result.event.key}: changed={changed}, "
            f"outside_grid={result.outside_grid}, "
            f"word_length_gt_2={np.count_nonzero(result.word_length > 2)}"
        )


if __name__ == "__main__":
    main()
