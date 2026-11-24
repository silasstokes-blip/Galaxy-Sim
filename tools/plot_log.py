#!/usr/bin/env python3
"""Plot a Galaxy Sim CSV log.

Columns: time_myr, separation_ly, bound_system_pct, bound_galaxy_pct

Two stacked panels share the time axis:
  top    - distance between the two central bodies (two-galaxy runs only)
  bottom - percent of stars bound to the system and bound to a galaxy

Usage:
  python tools/plot_log.py                      # newest CSV in logs/
  python tools/plot_log.py logs/galaxy_X.csv
  python tools/plot_log.py logs/galaxy_X.csv -o fig.png --dark
  python tools/plot_log.py --show               # open an interactive window
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

import numpy as np

# Validated categorical slots 1-3 (see the dataviz palette). Color follows the
# series, not its position, so each series keeps its color in every chart.
THEMES = {
    "light": dict(
        surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
        grid="#e1e0d9", axis="#c3c2b7",
        separation="#2a78d6", system="#eb6834", galaxy="#1baf7a",
    ),
    "dark": dict(
        surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781",
        grid="#2c2c2a", axis="#383835",
        separation="#3987e5", system="#d95926", galaxy="#199e70",
    ),
}


def load(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    data = np.atleast_1d(data)  # a one-row file comes back 0-d
    return {name: np.asarray(data[name], dtype=float) for name in data.dtype.names}


def newest_log() -> Path:
    logs = sorted(Path("logs").glob("*.csv"), key=lambda p: p.stat().st_mtime)
    if not logs:
        sys.exit("No CSV files in logs/. Run the sim with logging enabled first.")
    return logs[-1]


def end_labels(ax, items, theme, x_last):
    """Name + value at each line's right end, in text ink. Labels are spread
    apart vertically if the lines end close together."""
    ymin, ymax = ax.get_ylim()
    min_gap = 0.09 * (ymax - ymin)
    placed: list[float] = []
    for y, text in sorted(items, key=lambda it: it[0]):
        ty = y if not placed or y - placed[-1] >= min_gap else placed[-1] + min_gap
        placed.append(ty)
        ax.text(x_last, ty, "  " + text, va="center", ha="left", fontsize=9,
                color=theme["ink2"], clip_on=False)


def style_axes(ax, theme):
    ax.set_facecolor(theme["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["axis"])
    ax.tick_params(colors=theme["muted"], labelsize=9, length=3)
    ax.grid(axis="y", color=theme["grid"], linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.set_axisbelow(True)
    ax.yaxis.label.set_color(theme["ink2"])
    ax.xaxis.label.set_color(theme["ink2"])


def make_figure(d: dict[str, np.ndarray], title: str, mode: str):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    theme = THEMES[mode]
    t = d["time_myr"]
    has_sep = "separation_ly" in d and np.isfinite(d["separation_ly"]).any()
    # The distributed-mass mode logs only the system measure (galaxy column empty).
    has_gal = "bound_galaxy_pct" in d and np.isfinite(d["bound_galaxy_pct"]).any()

    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans", "sans-serif"]
    n_panels = 2 if has_sep else 1
    fig, axes = plt.subplots(
        n_panels, 1, sharex=True, figsize=(10, 7 if has_sep else 4.2),
        gridspec_kw=dict(height_ratios=[1, 1] if has_sep else [1], hspace=0.12),
    )
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(theme["surface"])
    for ax in axes:
        style_axes(ax, theme)

    thousands = FuncFormatter(lambda v, _p: f"{v:,.0f}")
    x_last = t[-1]

    # --- top: separation ---
    if has_sep:
        ax = axes[0]
        sep = d["separation_ly"]
        ax.plot(t, sep, color=theme["separation"], linewidth=2.0, solid_capstyle="round")
        ax.set_ylabel("Distance between centers (ly)")
        ax.yaxis.set_major_formatter(thousands)
        ax.set_ylim(0, np.nanmax(sep) * 1.08)

        # Summary line (not a marker on one pass: repeated passes are near-identical).
        # Minimum is sampled at the log interval, so it overestimates the true pericenter.
        subtitle = (f"closest sampled separation {np.nanmin(sep):,.0f} ly, "
                    f"farthest {np.nanmax(sep):,.0f} ly")
        ax.set_title(title, loc="left", fontsize=12, color=theme["ink"], pad=26)
        ax.text(0, 1.04, subtitle, transform=ax.transAxes, fontsize=9,
                color=theme["ink2"], va="bottom", ha="left")

    # --- bottom: bound percentages ---
    ax = axes[-1]
    sysp = d["bound_system_pct"]
    galp = d["bound_galaxy_pct"] if has_gal else None
    ax.plot(t, sysp, color=theme["system"], linewidth=2.0, label="Bound to system",
            solid_capstyle="round")
    if has_sep and has_gal:
        # In a single-galaxy run the two definitions are identical; draw one line.
        ax.plot(t, galp, color=theme["galaxy"], linewidth=2.0, label="Bound to a galaxy",
                solid_capstyle="round")
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Stars bound (%)")
    ax.set_xlabel("Simulated time (Myr)")
    ax.xaxis.set_major_formatter(thousands)
    ax.set_xlim(t[0], x_last)
    if not has_sep:
        ax.set_title(title, loc="left", fontsize=12, color=theme["ink"], pad=14)

    # Legend (>= 2 series) plus direct labels at the line ends.
    if has_sep and has_gal:
        leg = ax.legend(loc="upper right", frameon=False, fontsize=9, handlelength=1.6)
        for txt in leg.get_texts():
            txt.set_color(theme["ink2"])
        end_labels(ax, [(sysp[-1], f"System {sysp[-1]:.1f}%"),
                        (galp[-1], f"Galaxy {galp[-1]:.1f}%")], theme, x_last)
    else:
        end_labels(ax, [(sysp[-1], f"{sysp[-1]:.1f}%")], theme, x_last)

    fig.subplots_adjust(left=0.1, right=0.88, top=0.89 if has_sep else 0.88, bottom=0.11)
    return fig


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("csv", nargs="?", type=Path, help="log file (default: newest in logs/)")
    ap.add_argument("-o", "--output", type=Path, help="output PNG (default: next to the CSV)")
    ap.add_argument("--dark", action="store_true", help="dark theme")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--show", action="store_true", help="open an interactive window")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")

    path = args.csv or newest_log()
    d = load(path)
    if len(d.get("time_myr", [])) < 2:
        sys.exit(f"{path}: need at least 2 rows to plot (got {len(d.get('time_myr', []))}).")

    fig = make_figure(d, f"Galaxy Sim log: {path.name}", "dark" if args.dark else "light")

    out = args.output or path.with_suffix(".png")
    fig.savefig(out, dpi=args.dpi, facecolor=fig.get_facecolor())
    print(f"wrote {out}")
    if args.show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
