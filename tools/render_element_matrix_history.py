#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path

DEFAULT_INPUT_JSON = (
    Path(__file__).resolve().parents[1]
    / "test_data"
    / "doubleDislocation8x8Inspection"
    / "element_48_matrix_history.json"
)


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "_": r"\_",
        "&": r"\&",
        "%": r"\%",
        "#": r"\#",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def latex_scalar(value: float) -> str:
    if abs(value) < 1e-4:
        value = 0.0
    return format(value, ".3g")


def matrix_to_latex(matrix) -> str:
    return (
        r"$\left[\begin{smallmatrix}"
        + latex_scalar(matrix[0][0])
        + " & "
        + latex_scalar(matrix[0][1])
        + r"\\"
        + latex_scalar(matrix[1][0])
        + " & "
        + latex_scalar(matrix[1][1])
        + r"\end{smallmatrix}\right]$"
    )


def normalize_triangle(nodes):
    min_x = min(node[0] for node in nodes)
    max_x = max(node[0] for node in nodes)
    min_y = min(node[1] for node in nodes)
    max_y = max(node[1] for node in nodes)
    center_x = 0.5 * (min_x + max_x)
    center_y = 0.5 * (min_y + max_y)
    scale = max(max_x - min_x, max_y - min_y, 1e-12)
    normalized = []
    for x, y in nodes:
        normalized.append(
            (
                1.0 + 1.4 * ((x - center_x) / scale),
                1.0 + 1.4 * ((y - center_y) / scale),
            )
        )
    return normalized


def triangle_to_tikz(nodes, draw_style: str, node_color: str, labels) -> str:
    normalized = normalize_triangle(nodes)
    label_offsets = [(0.10, 0.12), (0.10, -0.12), (-0.14, 0.12)]
    parts = [
        r"\begin{tikzpicture}[x=1.25em,y=1.25em,baseline=(current bounding box.center)]",
        r"\path[use as bounding box] (0,0) rectangle (2,2);",
        r"\filldraw["
        + draw_style
        + "] "
        + f"({latex_scalar(normalized[0][0])},{latex_scalar(normalized[0][1])}) -- "
        + f"({latex_scalar(normalized[1][0])},{latex_scalar(normalized[1][1])}) -- "
        + f"({latex_scalar(normalized[2][0])},{latex_scalar(normalized[2][1])}) -- cycle;",
    ]
    for x, y in normalized:
        parts.append(
            r"\fill["
            + node_color
            + "] "
            + f"({latex_scalar(x)},{latex_scalar(y)}) circle[radius=0.06];"
        )
    for (x, y), label, (dx, dy) in zip(normalized, labels, label_offsets):
        parts.append(
            r"\node["
            + node_color
            + r",font=\tiny] at "
            + f"({latex_scalar(x + dx)},{latex_scalar(y + dy)})"
            + "{"
            + latex_escape(label)
            + "};"
        )
    parts.append(r"\end{tikzpicture}")
    return "".join(parts)


def row_cells(row) -> str:
    return " & ".join(
        [
            triangle_to_tikz(
                row["reference_nodes"],
                "draw=blue!70!black,fill=blue!12,line width=0.7pt",
                "blue!70!black",
                ("A", "B", "C"),
            ),
            triangle_to_tikz(
                row["current_nodes"],
                "draw=orange!80!black,fill=orange!18,line width=0.7pt",
                "orange!80!black",
                ("a", "b", "c"),
            ),
            matrix_to_latex(row["F_P"]),
            matrix_to_latex(row["H"]),
            matrix_to_latex(row["T"]),
        ]
    )


def build_latex(document) -> str:
    config = document["config"]
    element_index = document["element_index"]
    lines = [
        r"\documentclass{article}",
        r"\usepackage[margin=0.45in]{geometry}",
        r"\usepackage{amsmath}",
        r"\usepackage{array}",
        r"\usepackage{booktabs}",
        r"\usepackage{longtable}",
        r"\usepackage{pdflscape}",
        r"\usepackage{tikz}",
        r"\begin{document}",
        r"\begin{landscape}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{center}",
        rf"\textbf{{Double-dislocation inspection in a {config['rows']}$\times${config['cols']} mesh.}}\\",
        rf"Load step $\Delta\gamma = {latex_scalar(config['load_increment'])}$. Each gamma block is separated by a horizontal rule. Static steps use one row, while reconnecting steps use a Before row followed by an After row. Here $T = F_p H$.",
        r"\end{center}",
        r"\begin{longtable}{cc*{5}{>{\centering\arraybackslash}p{0.095\linewidth}}}",
        r"\caption{" + latex_escape(document["label"]) + r"}\\",
        r"\toprule",
        rf"$\gamma$ & State & \multicolumn{{5}}{{c}}{{Element {element_index}}} \\",
        r"\cmidrule(lr){3-7}",
        r" &  & Reference & Current & $F_p$ & $H$ & $T$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        rf"$\gamma$ & State & \multicolumn{{5}}{{c}}{{Element {element_index}}} \\",
        r"\cmidrule(lr){3-7}",
        r" &  & Reference & Current & $F_p$ & $H$ & $T$ \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]

    for row in document["rows"]:
        before_cells = row_cells(row["before"])
        after_cells = row_cells(row["after"])
        if row["element_changed"]:
            lines.append(
                f"${latex_scalar(row['gamma'])}$ & Before & {before_cells} \\\\"
            )
            lines.append(f" & After & {after_cells} \\\\")
        else:
            lines.append(
                f"${latex_scalar(row['gamma'])}$ & Static & {before_cells} \\\\"
            )
        lines.append(r"\midrule")

    lines.extend([r"\end{longtable}", r"\end{landscape}", r"\end{document}"])
    return "\n".join(lines) + "\n"


def render_pdf(tex_path: Path) -> None:
    pdflatex = shutil.which("pdflatex")
    if pdflatex is None:
        raise RuntimeError("pdflatex was not found on PATH")
    output_dir = tex_path.parent
    for _ in range(2):
        subprocess.run(
            [
                pdflatex,
                "-interaction=nonstopmode",
                "-halt-on-error",
                "-output-directory",
                str(output_dir),
                str(tex_path),
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render the element matrix history report from JSON."
    )
    parser.add_argument(
        "input_json",
        type=Path,
        nargs="?",
        default=DEFAULT_INPUT_JSON,
        help=f"Path to the inspection JSON. Defaults to {DEFAULT_INPUT_JSON}.",
    )
    parser.add_argument("--output-tex", type=Path, default=None)
    parser.add_argument(
        "--skip-pdf",
        action="store_true",
        help="Write the .tex file but do not run pdflatex.",
    )
    args = parser.parse_args()

    with args.input_json.open("r", encoding="utf-8") as handle:
        document = json.load(handle)

    tex_path = args.output_tex or args.input_json.with_suffix(".tex")
    tex_path.parent.mkdir(parents=True, exist_ok=True)
    tex_path.write_text(build_latex(document), encoding="utf-8")

    if not args.skip_pdf:
        render_pdf(tex_path)


if __name__ == "__main__":
    main()
