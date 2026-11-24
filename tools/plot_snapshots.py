#!/usr/bin/env python3
"""Scatter-plot galaxy_batch --snapshot files side by side.

  python tools/plot_snapshots.py out.png label1=snap1.csv label2=snap2.csv ...
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SURFACE, INK, MUTED, DOT = "#0b0b0b", "#f0efec", "#898781", "#f0efec"


def main() -> None:
    out = Path(sys.argv[1])
    items = [a.split("=", 1) for a in sys.argv[2:]]
    n = len(items)
    cols = min(n, 3)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 4.8 * rows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    lim = 0.0
    data = []
    for label, path in items:
        d = np.genfromtxt(path, delimiter=",", names=True)
        data.append((label, d["x"], d["y"]))
        lim = max(lim, float(np.percentile(np.hypot(d["x"], d["y"]), 99.5)))
    lim = min(lim * 1.05, 1.5e5)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (label, x, y) in zip(axes.ravel(), data):
        ax.axis("on")
        ax.set_facecolor(SURFACE)
        ax.scatter(x, y, s=0.6, c=DOT, alpha=0.55, linewidths=0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(label, color=INK, fontsize=10, loc="left")
        ax.tick_params(colors=MUTED, labelsize=7)
        for s in ax.spines.values():
            s.set_color("#383835")
    fig.tight_layout()
    fig.savefig(out, dpi=140, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
