#!/usr/bin/env python3
"""Compare FF and T_total on the same direct-history avalanche trajectory."""

from __future__ import annotations

from collections import Counter
import csv
from pathlib import Path

import matplotlib as mpl
import numpy as np

import analyze_plastic_slip_event as base


HERE = Path(__file__).resolve().parent
RUN_DIRECTORY = next((HERE / "mathcalF").glob("simpleShear*"))
BEFORE = 0.72
AFTER = 0.73

MEASURES = (
    base.Event(
        key="FF_same_simulation",
        label=r"Direct full history $FF$",
        directory=RUN_DIRECTORY,
        before=BEFORE,
        after=AFTER,
        total_field="Ftotal",  # Legacy VTU name for FF.
    ),
    base.Event(
        key="Ttotal_same_simulation",
        label=r"Reconstructed history $T_{\mathrm{total}}$",
        directory=RUN_DIRECTORY,
        before=BEFORE,
        after=AFTER,
        total_field="Ttotal",
    ),
)


def write_tables(results: tuple[base.EventResult, ...]) -> None:
    with (HERE / "plastic_slip_same_simulation_fields.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "measure",
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
            for index, values in enumerate(
                zip(
                    result.horizontal,
                    result.vertical,
                    result.word_length,
                    result.word_text,
                    result.two_block_order,
                    result.plastic_increment,
                )
            ):
                h, v, length, word, order, matrix = values
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

    with (HERE / "plastic_slip_same_simulation_heatmap_counts.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["measure", "H", "V", "element_count"])
        for result in results:
            counter = Counter(zip(result.horizontal, result.vertical))
            for (h, v), count in sorted(counter.items()):
                writer.writerow([result.event.key, int(h), int(v), int(count)])

    ff, ttotal = results
    ff_changed = (ff.horizontal != 0) | (ff.vertical != 0)
    ttotal_changed = (ttotal.horizontal != 0) | (ttotal.vertical != 0)
    same_h = ff.horizontal == ttotal.horizontal
    same_v = ff.vertical == ttotal.vertical
    same_pair = same_h & same_v
    same_increment = np.all(
        ff.plastic_increment == ttotal.plastic_increment, axis=(1, 2)
    )

    with (HERE / "plastic_slip_same_simulation_comparison_summary.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value"])
        rows = (
            ("elements", len(ff.horizontal)),
            ("FF_nonzero_net_slip_elements", int(np.count_nonzero(ff_changed))),
            (
                "Ttotal_nonzero_net_slip_elements",
                int(np.count_nonzero(ttotal_changed)),
            ),
            ("same_H_elements", int(np.count_nonzero(same_h))),
            ("same_V_elements", int(np.count_nonzero(same_v))),
            ("same_HV_pair_elements", int(np.count_nonzero(same_pair))),
            (
                "same_plastic_increment_matrix_elements",
                int(np.count_nonzero(same_increment)),
            ),
            (
                "FF_only_nonzero_elements",
                int(np.count_nonzero(ff_changed & ~ttotal_changed)),
            ),
            (
                "Ttotal_only_nonzero_elements",
                int(np.count_nonzero(ttotal_changed & ~ff_changed)),
            ),
            (
                "both_nonzero_elements",
                int(np.count_nonzero(ff_changed & ttotal_changed)),
            ),
            ("FF_outside_minus2_plus2", ff.outside_grid),
            ("Ttotal_outside_minus2_plus2", ttotal.outside_grid),
        )
        writer.writerows(rows)


def main() -> None:
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "figure.dpi": 160,
            "savefig.dpi": 240,
        }
    )
    results = tuple(base.analyze_event(event) for event in MEASURES)
    largest_count = max(int(np.max(result.heatmap)) for result in results)
    count_norm = mpl.colors.LogNorm(vmin=1, vmax=largest_count)
    count_cmap = base.heatmap_colormap()
    field_cmap = base.slip_colormap()
    field_norm = mpl.colors.BoundaryNorm(
        np.arange(base.GRID_MIN - 0.5, base.GRID_MAX + 1.5, 1.0),
        field_cmap.N,
    )

    base.save_individual_plots(
        results, count_norm, count_cmap, field_norm, field_cmap
    )
    title = (
        r"Same avalanche: integer plastic slips inferred from $FF$ and "
        r"$T_{\mathrm{total}}$"
    )
    base.save_composite(
        results,
        count_norm,
        count_cmap,
        field_norm,
        field_cmap,
        stem="plastic_slip_same_simulation_comparison",
        title=title,
    )
    base.save_composite(
        results,
        count_norm,
        count_cmap,
        field_norm,
        field_cmap,
        suffix="_dark",
        dark=True,
        stem="plastic_slip_same_simulation_comparison",
        title=title,
    )
    write_tables(results)

    ff, ttotal = results
    same_pair = (ff.horizontal == ttotal.horizontal) & (
        ff.vertical == ttotal.vertical
    )
    print(
        f"FF changed: {np.count_nonzero((ff.horizontal != 0) | (ff.vertical != 0))}"
    )
    print(
        "T_total changed: "
        f"{np.count_nonzero((ttotal.horizontal != 0) | (ttotal.vertical != 0))}"
    )
    print(f"Same H,V pair: {np.count_nonzero(same_pair)}/{len(same_pair)}")


if __name__ == "__main__":
    main()
