#!/usr/bin/env python3
"""Regenerate the 70x70 mean-history comparison with final notation."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DATA = np.genfromtxt(
    HERE / "Ftotal_vs_Tfull_70x70_simple_shear.csv",
    delimiter=",",
    names=True,
)

fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.1), sharex=True, constrained_layout=True)
for ax, component in zip(axes.flat, ("11", "12", "21", "22")):
    ax.plot(
        DATA["gamma"], DATA[f"Ftotal{component}_mean"],
        "-o", color="#1769aa", linewidth=1.7, markersize=2.7,
        label=rf"$\langle\mathcal{{F}}_{{{component}}}\rangle$",
    )
    ax.plot(
        DATA["gamma"], DATA[f"Tfull{component}_mean"],
        "--s", color="#d95f02", linewidth=1.5, markersize=2.4,
        label=rf"$\langle(T_{{\mathrm{{total}}}})_{{{component}}}\rangle$",
    )
    ax.set_title(f"Component {component}")
    ax.set_ylabel("element-wise mean")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=8)

for ax in axes[-1]:
    ax.set_xlabel(r"imposed shear $\gamma$")

fig.suptitle(
    r"Transported history $\mathcal{F}$ vs. reconstructed "
    r"$T_{\mathrm{total}}=F_E T$",
    fontsize=13,
)
fig.savefig(HERE / "Ftotal_vs_Ttotal_70x70_simple_shear.png", dpi=300)
fig.savefig(HERE / "Ftotal_vs_Ttotal_70x70_simple_shear.pdf")
plt.close(fig)
