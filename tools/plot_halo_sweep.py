import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = os.path.join(os.path.dirname(__file__), "..", "logs", "halo_sweep")
runs = [("90", "#0b3d7a"), ("50", "#2a78d6"), ("20", "#7fb0ea"), ("00", "#c3d9f5")]
fig, ax = plt.subplots(figsize=(8, 4.6), facecolor="#fcfcfb")
ax.set_facecolor("#fcfcfb")
for tag, col in runs:
    t, b = [], []
    for r in csv.DictReader(open(os.path.join(d, f"halo_{tag}.csv"))):
        t.append(float(r["time_myr"])); b.append(float(r["bound_system_pct"]))
    ax.plot(t, b, color=col, lw=2, label=f"{tag}% halo")
    ax.text(t[-1] + 20, b[-1], f"{tag}%", va="center", fontsize=9, color="#333")
ax.set_xlabel("time (Myr)"); ax.set_ylabel("disk particles bound to system (%)")
ax.set_title("Bound fraction vs halo mass fraction (2 Gyr, 2 galaxies)", loc="left", fontsize=11)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e6e3", lw=0.8)
ax.legend(frameon=False, loc="lower left", ncol=4)
ax.set_xlim(0, 2200)
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(__file__), "..", "docs", "diagrams", "7_halo_sweep.png"), dpi=150)
