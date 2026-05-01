#!/usr/bin/env python3
"""
Plot HyCal W (PbWO4) module statistics map from single-cluster events.
Each module is colored by event count and annotated with its module ID.
"""

import json
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ── Parameters ──────────────────────────────────────────────────────────────
root_file  = 'data/prad_024004_recon.root'
geo_file   = 'hycal_modules.json'
output     = 'hycal_w_stats.png'
cmap_name  = 'YlOrRd'
use_log    = True          # set True to use log color scale
title      = 'HyCal W (PbWO4) Single-Cluster Event Statistics'
ann_fontsize = 7.5          # font size for module ID labels
# ────────────────────────────────────────────────────────────────────────────

# ── 1. Count single-cluster events per W module ──────────────────────────────
# cl_center encodes module ID: mod_id = cl_center - 1000
# W modules: mod_id in [0, 1155]  →  module name "W{mod_id+1}"
mod_counts = {}   # key: "W{n}", value: int count

print(f"Reading {root_file} ...")
with uproot.open(root_file) as f:
    tree = f["recon"]
    total = tree.num_entries
    print(f"  Total entries: {total}")

    for batch in tree.iterate(
        ["n_clusters", "cl_center"],
        step_size=50000,
        library="np"
    ):
        n_cl     = batch["n_clusters"]
        cl_cen   = batch["cl_center"]    # 1-D array of objects (jagged)
        for nc, centers in zip(n_cl, cl_cen):
            if nc != 1:
                continue
            mid = int(centers[0]) - 1000
            if 0 < mid <= 1156:
                name = f"W{mid}"
                mod_counts[name] = mod_counts.get(name, 0) + 1

print(f"  W modules with hits: {len(mod_counts)}")

# ── 2. Load geometry ──────────────────────────────────────────────────────────
with open(geo_file) as f:
    all_modules = json.load(f)

w_modules = [m for m in all_modules if m["t"] == "PbWO4"]
print(f"  W modules in geometry: {len(w_modules)}")

# ── 2b. Load DAQ map ──────────────────────────────────────────────────────────
with open("hycal_daq_map.json") as f:
    daq_list = json.load(f)
daq_map = {entry["name"]: entry for entry in daq_list}  # name → {crate, slot, channel}

# ── 3. Color normalization ────────────────────────────────────────────────────
counts = [mod_counts.get(m["n"], 0) for m in w_modules]

# ── 3a. Report zero-count (dead/missing) modules ─────────────────────────────
dead = [m for m in w_modules if mod_counts.get(m["n"], 0) == 0]
print(f"\n{'─'*60}")
print(f"  Zero-count W modules: {len(dead)}")
print(f"  {'Module':<10} {'x(mm)':>8} {'y(mm)':>8}  {'DAQ info'}")
print(f"  {'─'*56}")
for m in sorted(dead, key=lambda m: int(m["n"][1:])):
    name = m["n"]
    daq  = daq_map.get(name)
    if daq:
        daq_str = f"crate={daq['crate']}  slot={daq['slot']:>2}  ch={daq['channel']:>2}"
    else:
        daq_str = "not in DAQ map"
    print(f"  {name:<10} {m['x']:>8.2f} {m['y']:>8.2f}  {daq_str}")
print(f"{'─'*60}\n")
lo = min(c for c in counts if c > 0) if any(c > 0 for c in counts) else 1
hi = max(counts) if max(counts) > 0 else 1

if use_log:
    norm = mcolors.LogNorm(vmin=max(lo, 1), vmax=hi)
else:
    norm = mcolors.Normalize(vmin=lo, vmax=hi)

cmap = plt.get_cmap(cmap_name)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array(counts)

# ── 4. Plot ───────────────────────────────────────────────────────────────────
xs = [m["x"] for m in w_modules]
ys = [m["y"] for m in w_modules]
margin = 30
x_lo, x_hi = min(xs) - margin, max(xs) + margin
y_lo, y_hi = min(ys) - margin, max(ys) + margin

fig, ax = plt.subplots(figsize=(16, 14), dpi=160,
                       gridspec_kw=dict(left=0.08, right=0.88, bottom=0.06, top=0.94))

for m in w_modules:
    name = m["n"]
    count = mod_counts.get(name, 0)
    fc = sm.to_rgba(count) if count > 0 else "white"
    rx = m["x"] - m["sx"] / 2
    ry = m["y"] - m["sy"] / 2
    ax.add_patch(Rectangle((rx, ry), m["sx"], m["sy"],
                            fc=fc, ec="black", lw=0.25, zorder=1))
    # annotate with numeric ID (strip leading "W")
    mod_id_label = name[1:]   # e.g. "W555" → "555"
    ax.text(m["x"], m["y"], mod_id_label,
            ha="center", va="center",
            fontsize=ann_fontsize, color="black", zorder=2)

ax.set_xlim(x_lo, x_hi)
ax.set_ylim(y_lo, y_hi)
ax.set_aspect("equal")
ax.tick_params(labelsize=16)
ax.set_xlabel("X (mm)", fontsize=18)
ax.set_ylabel("Y (mm)", fontsize=18)
ax.set_title(title, fontsize=20)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)
cbar = plt.colorbar(sm, cax=cax)
cax.tick_params(labelsize=14)
cbar.set_label("Single-cluster event count", fontsize=16)

fig.savefig(output)
print(f"Saved: {output}")
