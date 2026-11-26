#!/usr/bin/env python3
"""Memory diagrams for the Barnes-Hut quadtree and the galaxy sim's structure.

A small Python mirror of src/sim/quadtree.cpp (same quadrant bits, same
insertion, same opening criterion including the "target inside the cell" rule),
used only to draw pictures. Writes PNGs to docs/diagrams/.

  python tools/plot_quadtree.py
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parent.parent / "docs" / "diagrams"

# Palette: validated categorical slots 1-3 and neutrals (light surface).
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_LIGHT = "#cde2fb"

QUAD_NAMES = ["0: SW", "1: SE", "2: NW", "3: NE"]  # bit0 = x > centre, bit1 = y > centre


# --------------------------------------------------------------------------
# Mirror of the C++ tree
# --------------------------------------------------------------------------
class Node:
    __slots__ = ("px", "py", "half", "mass", "cx", "cy", "child", "body", "depth")

    def __init__(self, px, py, half, depth=0):
        self.px, self.py, self.half = px, py, half
        self.mass = 0.0
        self.cx = self.cy = 0.0
        self.child = -1  # index of first of 4 children
        self.body = -1  # -1 empty, -2 merged, >=0 index
        self.depth = depth


class Tree:
    MAX_DEPTH = 48

    def __init__(self, xs, ys, ms=None, root=None):
        self.x, self.y = xs, ys
        self.m = ms or [1.0] * len(xs)
        if root is None:
            lo_x, hi_x, lo_y, hi_y = min(xs), max(xs), min(ys), max(ys)
            half = 0.5 * max(hi_x - lo_x, hi_y - lo_y) * (1 + 1e-9)
            root = (0.5 * (lo_x + hi_x), 0.5 * (lo_y + hi_y), half)
        self.nodes = [Node(*root)]
        for i in range(len(xs)):
            self.insert(i)
        for nd in self.nodes:
            if nd.mass > 0:
                nd.cx /= nd.mass
                nd.cy /= nd.mass

    @staticmethod
    def quadrant(px, py, x, y):
        return (1 if x > px else 0) | (2 if y > py else 0)

    def insert(self, i):
        xi, yi, mi = self.x[i], self.y[i], self.m[i]
        node, add, depth = 0, True, 0
        while True:
            n = self.nodes[node]
            if add:
                n.mass += mi
                n.cx += mi * xi
                n.cy += mi * yi
            if n.child >= 0:
                node = n.child + self.quadrant(n.px, n.py, xi, yi)
                add = True
                depth += 1
                continue
            if n.body == -1:
                n.body = i
                return
            if n.body == -2 or depth >= self.MAX_DEPTH:
                n.body = -2
                return
            j = n.body
            h = 0.5 * n.half
            first = len(self.nodes)
            for q in range(4):
                self.nodes.append(
                    Node(n.px + (h if q & 1 else -h), n.py + (h if q & 2 else -h), h, n.depth + 1)
                )
            n.child, n.body = first, -1
            cj = self.nodes[first + self.quadrant(n.px, n.py, self.x[j], self.y[j])]
            cj.mass, cj.cx, cj.cy, cj.body = self.m[j], self.m[j] * self.x[j], self.m[j] * self.y[j], j
            add = False

    def walk(self, self_idx, x, y, theta):
        """Return (accepted, opened): lists of node indices, like the C++ walk."""
        accepted, opened = [], []
        stack = [0]
        t2 = theta * theta
        while stack:
            k = stack.pop()
            n = self.nodes[k]
            if n.mass <= 0:
                continue
            dx, dy = n.cx - x, n.cy - y
            r2 = dx * dx + dy * dy
            if n.child < 0:
                if n.body == self_idx:
                    continue
                accepted.append(k)
                continue
            s = 2 * n.half
            inside = abs(x - n.px) <= n.half and abs(y - n.py) <= n.half
            if s * s < t2 * r2 and not inside:
                accepted.append(k)
            else:
                opened.append(k)
                stack.extend(n.child + q for q in range(4))
        return accepted, opened


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def style(ax, title=None):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)


def draw_cells(ax, tree, lw=0.8, color=AXIS):
    for n in tree.nodes:
        ax.add_patch(
            Rectangle((n.px - n.half, n.py - n.half), 2 * n.half, 2 * n.half,
                      fill=False, ec=color, lw=lw)
        )


def draw_bodies(ax, tree, labels=None, size=26, color=BLUE):
    ax.scatter(tree.x, tree.y, s=size, c=color, zorder=3, edgecolors=SURFACE, linewidths=0.6)
    if labels:
        for i, name in enumerate(labels):
            ax.annotate(name, (tree.x[i], tree.y[i]), xytext=(5, 5), textcoords="offset points",
                        color=INK, fontsize=9, fontweight="bold")


def fit_square(ax, tree, pad=0.03):
    r = tree.nodes[0]
    ax.set_xlim(r.px - r.half - pad, r.px + r.half + pad)
    ax.set_ylim(r.py - r.half - pad, r.py + r.half + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def layout_tree(tree):
    """x positions: leaves get consecutive slots, parents sit above their children."""
    pos = {}
    counter = [0]

    def rec(k, depth):
        n = tree.nodes[k]
        if n.child < 0:
            pos[k] = (counter[0], -depth)
            counter[0] += 1
        else:
            xs = [rec(n.child + q, depth + 1) for q in range(4)]
            pos[k] = (sum(xs) / 4.0, -depth)
        return pos[k][0]

    rec(0, 0)
    return pos


def draw_tree(ax, tree, names=None, show_mass=True, fs=8, ms=430):
    """Nodes are scatter markers (fixed size in points), so they stay round
    whatever the axes' aspect ratio."""
    pos = layout_tree(tree)
    for k, n in enumerate(tree.nodes):
        if n.child >= 0:
            for q in range(4):
                c = n.child + q
                ax.plot([pos[k][0], pos[c][0]], [pos[k][1], pos[c][1]], color=AXIS, lw=1, zorder=1)
    for k, n in enumerate(tree.nodes):
        x, y = pos[k]
        if n.child >= 0:
            ax.scatter([x], [y], s=ms, fc=BLUE_LIGHT, ec=BLUE, lw=1.3, zorder=2)
            if show_mass:
                ax.text(x, y, f"{int(round(n.mass))}", ha="center", va="center", color=INK,
                        fontsize=fs, fontweight="bold", zorder=3)
        elif n.body >= 0:
            ax.scatter([x], [y], s=ms, fc=BLUE, ec=BLUE, zorder=2)
            nm = names[n.body] if names else str(n.body)
            ax.text(x, y, nm, ha="center", va="center", color="#ffffff", fontsize=fs,
                    fontweight="bold", zorder=3)
        else:
            ax.scatter([x], [y], s=ms * 0.32, marker="s", fc=SURFACE, ec=MUTED, lw=0.9, zorder=2)
    maxd = max(-p[1] for p in pos.values())
    ax.set_xlim(-0.8, max(p[0] for p in pos.values()) + 0.8)
    ax.set_ylim(-maxd - 0.6, 0.6)
    ax.axis("off")


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def fig_structure():
    pts = {
        "A": (0.12, 0.10), "B": (0.30, 0.20), "C": (0.72, 0.15), "D": (0.88, 0.40),
        "E": (0.60, 0.62), "F": (0.70, 0.80), "G": (0.22, 0.74), "H": (0.40, 0.88),
        "I": (0.04, 0.55), "J": (0.96, 0.96),
    }
    names = list(pts)
    xs = [pts[k][0] for k in names]
    ys = [pts[k][1] for k in names]
    tree = Tree(xs, ys, root=(0.5, 0.5, 0.5))

    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.4), gridspec_kw={"width_ratios": [1, 1.7]})
    fig.patch.set_facecolor(SURFACE)
    style(a, "Space: each cell splits into 4 until it holds at most 1 body")
    draw_cells(a, tree)
    draw_bodies(a, tree, names)
    fit_square(a, tree)
    # quadrant labels on the root, hugging the centre so they never hit a body
    for q, (dx, dy) in enumerate([(-1, -1), (1, -1), (-1, 1), (1, 1)]):
        a.text(0.5 + 0.02 * dx, 0.5 + 0.02 * dy, QUAD_NAMES[q], color=MUTED, fontsize=8,
               ha="right" if dx < 0 else "left", va="top" if dy < 0 else "bottom")

    style(b, "Tree: every split creates 4 children (empty ones are dashed squares)")
    draw_tree(b, tree, names)
    b.text(0.0, -0.02, "Circles with a number = internal cells (number = total mass inside, "
           "the COM is stored too).\nFilled circles = leaves holding one body. "
           "Children order is the quadrant number: 0 SW, 1 SE, 2 NW, 3 NE.",
           transform=b.transAxes, color=INK2, fontsize=9, va="top")
    fig.tight_layout()
    fig.savefig(OUT / "1_tree_structure.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_insertion():
    pts = [("A", (0.12, 0.10)), ("J", (0.90, 0.90)), ("B", (0.30, 0.20)), ("C", (0.72, 0.15))]
    notes = [
        "Insert A.\nRoot is an empty leaf,\nso it just takes A.",
        "Insert J.\nRoot is an occupied leaf:\nsplit into 4, push A down\ninto its quadrant (SW),\nthen put J in NE.",
        "Insert B.\nB walks down to SW,\nwhich holds A: split SW,\nplace A and B in its\nchildren (SW again here).",
        "Insert C.\nC walks down to SE,\nan empty leaf: takes it.\nMasses and COMs along\nthe path are updated.",
    ]
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.4), gridspec_kw={"height_ratios": [1.5, 1]})
    fig.patch.set_facecolor(SURFACE)
    for k in range(4):
        sub = pts[: k + 1]
        names = [s[0] for s in sub]
        t = Tree([s[1][0] for s in sub], [s[1][1] for s in sub], root=(0.5, 0.5, 0.5))
        top, bot = axes[0][k], axes[1][k]
        style(top, f"Step {k + 1}")
        draw_cells(top, t)
        draw_bodies(top, t, names)
        fit_square(top, t)
        top.text(0.02, -0.04, notes[k], transform=top.transAxes, color=INK2, fontsize=9, va="top")
        bot.set_facecolor(SURFACE)
        draw_tree(bot, t, names)
    fig.suptitle("Building the tree: insert one body at a time, splitting an occupied leaf",
                 color=INK, fontsize=13, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.subplots_adjust(hspace=0.45)
    fig.savefig(OUT / "2_insertion_steps.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def blob(n, seed):
    rnd = random.Random(seed)
    xs, ys = [], []
    for i in range(n):
        if i < 0.7 * n:
            r = 0.30 * math.sqrt(rnd.random())
            th = 6.2832 * rnd.random()
            xs.append(0.5 + r * math.cos(th)); ys.append(0.5 + r * math.sin(th))
        else:
            xs.append(0.5 + rnd.gauss(0, 0.28)); ys.append(0.5 + rnd.gauss(0, 0.28))
    return xs, ys


def fig_opening():
    n = 400
    xs, ys = blob(n, 4)
    tree = Tree(xs, ys)
    target = min(range(n), key=lambda i: (xs[i] - 0.68) ** 2 + (ys[i] - 0.40) ** 2)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.9))
    fig.patch.set_facecolor(SURFACE)
    for ax, theta in zip(axes, (1.0, 0.5, 0.2)):
        acc, _ = tree.walk(target, xs[target], ys[target], theta)
        style(ax, f"theta = {theta}: {len(acc)} interactions instead of {n - 1}")
        for k in acc:
            nd = tree.nodes[k]
            if nd.child < 0:
                continue
            ax.add_patch(Rectangle((nd.px - nd.half, nd.py - nd.half), 2 * nd.half, 2 * nd.half,
                                   fc=BLUE_LIGHT, ec=BLUE, lw=0.7, alpha=0.75))
        singles = [tree.nodes[k] for k in acc if tree.nodes[k].child < 0]
        ax.scatter([s.cx for s in singles], [s.cy for s in singles], s=7, c=BLUE, zorder=3)
        coms = [tree.nodes[k] for k in acc if tree.nodes[k].child >= 0]
        ax.scatter([c.cx for c in coms], [c.cy for c in coms], s=14, c=INK, zorder=4, marker="+")
        ax.scatter([xs[target]], [ys[target]], s=120, c=ORANGE, zorder=5, edgecolors=INK, linewidths=1)
        fit_square(ax, tree, pad=0.02)
    fig.text(0.01, 0.015,
             "Orange = the body whose force is being computed. Blue boxes = cells treated as ONE point mass at "
             "their centre of mass (black +). Single dots = bodies used individually.\n"
             "A cell is accepted when (cell width / distance to its COM) < theta, and never when the target is "
             "inside it. Smaller theta opens more cells: more accurate, slower.",
             color=INK2, fontsize=9, va="bottom")
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(OUT / "3_opening_angle.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_scaling():
    sizes = [250, 500, 1000, 2000, 4000, 8000, 16000]
    series = {1.0: [], 0.5: [], 0.2: []}
    for n in sizes:
        xs, ys = blob(n, 11)
        t = Tree(xs, ys)
        rnd = random.Random(2)
        targets = rnd.sample(range(n), 25)
        for th in series:
            series[th].append(sum(len(t.walk(i, xs[i], ys[i], th)[0]) for i in targets) / len(targets))
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    fig.patch.set_facecolor(SURFACE)
    style(ax, "Interactions per body: direct sum grows with N, the tree barely does")
    ax.plot(sizes, [n - 1 for n in sizes], color=ORANGE, lw=2.2, marker="o", ms=5, label="direct sum (N - 1)")
    cols = {1.0: AQUA, 0.5: BLUE, 0.2: "#6250d6"}
    for th, vals in series.items():
        ax.plot(sizes, vals, color=cols[th], lw=2.2, marker="o", ms=5, label=f"Barnes-Hut, theta = {th}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("number of bodies N", color=INK2); ax.set_ylabel("interactions per body", color=INK2)
    ax.grid(color=GRID, lw=0.7)
    ax.legend(frameon=False, labelcolor=INK, fontsize=9)
    ax.text(0.0, -0.17, "Measured with a Python mirror of the C++ walk on a clumpy galaxy-like cloud (25 sample "
            "bodies per N).\nTotal work = N bodies x interactions per body: O(N^2) direct vs roughly O(N log N) tree.",
            transform=ax.transAxes, color=INK2, fontsize=8.5, va="top")
    fig.tight_layout()
    fig.savefig(OUT / "4_cost_scaling.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def box(ax, x, y, w, h, text, fc=BLUE_LIGHT, ec=BLUE, fs=9.5, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=1.4))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=INK, fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.4)


def arrow(ax, p, q, color=INK2, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=14, color=color, lw=1.4,
                                 connectionstyle=f"arc3,rad={rad}"))


def fig_pipeline():
    fig, ax = plt.subplots(figsize=(16, 7.2))
    fig.patch.set_facecolor(SURFACE)
    ax.set_xlim(0, 16); ax.set_ylim(0, 7.2); ax.axis("off")
    ax.text(0.2, 6.85, "One timestep of the N-body mode (kick-drift-kick leapfrog)", color=INK, fontsize=14,
            fontweight="bold")
    ax.text(0.2, 6.45, "Start of the step: positions x, velocities v, and accelerations a cached from the "
            "previous step.", color=INK2, fontsize=10)
    y1, h = 4.35, 1.2
    box(ax, 0.3, y1, 2.9, h, "1. Half kick\nv += a * dt / 2")
    box(ax, 3.9, y1, 2.9, h, "2. Drift\nx += v * dt")
    box(ax, 7.5, y1, 3.1, h, "3. Rebuild the quadtree\nat the new positions\nO(N log N)", fc="#fbe0d4", ec=ORANGE)
    box(ax, 11.3, y1, 4.3, h, "4. Walk the tree for every body\n(parallel, OpenMP) gives new a\nO(log N) each, theta = 0.5",
        fc="#fbe0d4", ec=ORANGE)
    for a, b in [((3.2, y1 + h / 2), (3.9, y1 + h / 2)), ((6.8, y1 + h / 2), (7.5, y1 + h / 2)),
                 ((10.6, y1 + h / 2), (11.3, y1 + h / 2))]:
        arrow(ax, a, b)
    y2 = 2.0
    box(ax, 11.3, y2, 4.3, h, "5. Half kick with the NEW a\nv += a * dt / 2")
    box(ax, 6.9, y2, 3.6, h, "6. t += dt\nIf peak |a| > 4 x its reference:\ndt = dt / 2 (never grows)", fs=9)
    box(ax, 2.0, y2, 3.9, h, "Loop: the new a is cached,\nso only ONE force pass\nper step", fc="#d8f2e8", ec=AQUA, fs=9)
    arrow(ax, (13.45, y1), (13.45, y2 + h))
    arrow(ax, (11.3, y2 + h / 2), (10.5, y2 + h / 2))
    arrow(ax, (6.9, y2 + h / 2), (5.9, y2 + h / 2))
    arrow(ax, (3.9, y2 + h), (1.7, y1), rad=-0.35)
    notes = [
        "Why leapfrog: second order (error / 4 when dt halves, measured), time-reversible and symplectic,",
        "so energy error stays bounded and oscillates instead of drifting steadily like Euler.",
        "What the tree changes: only step 4. Everything else is the same as a plain integrator.",
        "Known imperfections: the monopole-only tree adds an energy drift that scales with theta^2;",
        "dt is chosen from the start state (sqrt(eps / max |a|)) and the ratchet only ever makes it smaller.",
    ]
    for i, s in enumerate(notes):
        ax.text(0.3, 1.25 - 0.27 * i, s, color=INK2, fontsize=9.5)
    fig.savefig(OUT / "5_timestep_pipeline.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_galaxy_model():
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.4), gridspec_kw={"width_ratios": [1.2, 1]})
    fig.patch.set_facecolor(SURFACE)
    # schematic, R = 1
    R, hole, core, rh = 1.0, 0.2, 0.25, 2.0
    a.set_facecolor(SURFACE)
    rnd = random.Random(1)
    hx, hy = [], []
    for _ in range(700):  # halo: surface density ~ 1/sqrt(R^2 + core^2)
        s = core + rnd.random() * (math.sqrt(rh * rh + core * core) - core)
        r = math.sqrt(max(0.0, s * s - core * core))
        th = 6.2832 * rnd.random()
        hx.append(r * math.cos(th)); hy.append(r * math.sin(th))
    a.scatter(hx, hy, s=3, c=BLUE, alpha=0.35, zorder=1)
    dx, dy = [], []
    for _ in range(500):
        r = math.sqrt(hole * hole + rnd.random() * (R * R - hole * hole))
        th = 6.2832 * rnd.random()
        dx.append(r * math.cos(th)); dy.append(r * math.sin(th))
    a.scatter(dx, dy, s=5, c=INK2, zorder=2)
    a.add_patch(Circle((0, 0), R, fill=False, ec=INK, lw=1.4, zorder=3))
    a.add_patch(Circle((0, 0), hole, fill=False, ec=ORANGE, lw=1.4, ls="--", zorder=3))
    a.add_patch(Circle((0, 0), core, fill=False, ec=AQUA, lw=1.4, ls=":", zorder=3))
    a.add_patch(Circle((0, 0), rh, fill=False, ec=BLUE, lw=1.4, ls="--", zorder=3))
    a.scatter([0], [0], s=70, c="#000000", zorder=5)
    lab = dict(fontsize=9.5, fontweight="bold", va="center")
    arr = lambda c: dict(arrowstyle="-", color=c, lw=1.0, shrinkA=2, shrinkB=2)
    for txt, tip, at, col in [
        ("halo truncation 2 R", (-1.41, 1.41), (-2.15, 2.0), BLUE),
        ("disk edge R", (0.71, 0.71), (1.15, 2.05), INK),
        ("black hole", (0.0, 0.0), (-2.15, -1.7), INK),
        ("halo core a = 0.25 R", (-0.18, -0.18), (-2.15, -2.05), AQUA),
        ("hole (optional, 20% shown)", (0.14, -0.14), (0.35, -2.05), ORANGE),
    ]:
        a.annotate(txt, tip, xytext=at, color=col, arrowprops=arr(col), **lab)
    a.set_xlim(-2.3, 2.3); a.set_ylim(-2.3, 2.3); a.set_aspect("equal"); a.axis("off")
    a.set_title("One galaxy (not to scale): grey = visible disk particles, blue = dark-matter halo",
                color=INK, fontsize=11, loc="left")

    style(b, "Mass and particle budget at the defaults")
    b.set_xlim(0, 100); b.set_ylim(0, 4.2); b.set_yticks([]); b.spines[["left", "top", "right"]].set_visible(False)
    b.set_xlabel("percent of galaxy mass", color=INK2)
    segs = [(0, 1, INK, "BH 1%"), (1, 9, INK2, "disk 9%"), (10, 90, BLUE, "dark-matter halo 90%")]
    for x0, w, col, lab in segs:
        b.barh(3.2, w, left=x0, height=0.7, color=col, edgecolor=SURFACE, linewidth=2)
    b.annotate("black hole 1%", (0.5, 2.85), xytext=(0, 2.1), color=INK, fontsize=9.5,
               arrowprops=dict(arrowstyle="-", color=INK2, lw=1), ha="left")
    b.annotate("visible disk 9%", (5.5, 2.85), xytext=(22, 2.1), color=INK, fontsize=9.5,
               arrowprops=dict(arrowstyle="-", color=INK2, lw=1), ha="left")
    b.text(55, 3.2, "dark-matter halo 90%", color="#ffffff", fontsize=11, ha="center", va="center", fontweight="bold")
    b.text(0, 1.5, "Particles (per galaxy, default): 5,000 disk + 15,000 halo (3 per disk particle) + 1 black hole.\n"
           "Halo particles are heavier, so 90% of the mass sits in 75% of the particles.\n\n"
           "Every particle attracts every other one (Barnes-Hut tree).\n"
           "Disk speeds and halo pressure come from the measured in-plane force;\n"
           "each galaxy is then relaxed alone (800 Myr) before t = 0.",
           color=INK2, fontsize=9.5, va="top")
    fig.tight_layout()
    fig.savefig(OUT / "6_galaxy_model.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for f in (fig_structure, fig_insertion, fig_opening, fig_scaling, fig_pipeline, fig_galaxy_model):
        f()
        print("wrote", f.__name__)


if __name__ == "__main__":
    main()
