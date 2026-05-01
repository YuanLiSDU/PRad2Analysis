#!/usr/bin/env python3
"""
Plot HyCal W (PbWO4) module map of calibration factor difference:
    delta = (fadc_1 - factor_1) / factor_1  [relative difference, %]
Each module is colored by the relative difference and annotated with its ID.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ── Parameters ──────────────────────────────────────────────────────────────
file_ref   = 'calibration_factor_1.json'
file_fadc  = 'calibration_factor_fadc_1.json'
geo_file   = '../hycal_modules.json'
output     = 'hycal_calib_diff.png'
cmap_name  = 'RdBu_r'          # red = fadc > ref, blue = fadc < ref
vmin       = -30.0              # color scale lower limit (%); None = auto
vmax       =  30.0              # color scale upper limit (%); None = auto
title      = 'Calibration Factor Relative Difference: (FADC − Physics) / Physics  [%]'
ann_fontsize = 7.5
# ────────────────────────────────────────────────────────────────────────────

# ── 1. Load calibration factors ──────────────────────────────────────────────
with open(file_ref) as f:
    ref_list = json.load(f)
with open(file_fadc) as f:
    fadc_list = json.load(f)

ref_map  = {e["name"]: e["factor"] for e in ref_list}
fadc_map = {e["name"]: e["factor"] for e in fadc_list}

# ── 2. Load geometry (W modules only) ────────────────────────────────────────
with open(geo_file) as f:
    all_modules = json.load(f)
w_modules = [m for m in all_modules if m["t"] == "PbWO4"]
print(f"W modules in geometry: {len(w_modules)}")

# ── 3. Compute relative difference [%] ───────────────────────────────────────
diffs = {}
skipped = []
for m in w_modules:
    name = m["n"]
    r = ref_map.get(name, 0.0)
    f = fadc_map.get(name, 0.0)
    if r == 0.0:
        skipped.append(name)
        diffs[name] = None          # no reference → no color
    else:
        diffs[name] = (f - r) / r * 100.0   # percent

valid_vals = [v for v in diffs.values() if v is not None]
print(f"Modules with valid reference factor: {len(valid_vals)}")
print(f"Skipped (zero reference factor):     {len(skipped)}")
print(f"Diff range: {min(valid_vals):.3f}% to {max(valid_vals):.3f}%")
print(f"Mean: {np.mean(valid_vals):.3f}%   Std: {np.std(valid_vals):.3f}%")

# ── 4. Color normalization ────────────────────────────────────────────────────
_vmin = vmin if vmin is not None else min(valid_vals)
_vmax = vmax if vmax is not None else max(valid_vals)
norm = mcolors.Normalize(vmin=_vmin, vmax=_vmax)

cmap = plt.get_cmap(cmap_name)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array(valid_vals)

# ── 5. Plot ───────────────────────────────────────────────────────────────────
xs = [m["x"] for m in w_modules]
ys = [m["y"] for m in w_modules]
margin = 30
x_lo, x_hi = min(xs) - margin, max(xs) + margin
y_lo, y_hi = min(ys) - margin, max(ys) + margin

fig, ax = plt.subplots(figsize=(16, 14), dpi=160,
                       gridspec_kw=dict(left=0.08, right=0.88, bottom=0.06, top=0.94))

for m in w_modules:
    name = m["n"]
    val  = diffs[name]
    if val is None:
        fc = "lightgrey"       # no reference factor
    else:
        fc = sm.to_rgba(val)
    rx = m["x"] - m["sx"] / 2
    ry = m["y"] - m["sy"] / 2
    ax.add_patch(Rectangle((rx, ry), m["sx"], m["sy"],
                            fc=fc, ec="black", lw=0.25, zorder=1))
    mod_id_label = name[1:]    # "W555" → "555"
    ax.text(m["x"], m["y"], mod_id_label,
            ha="center", va="center",
            fontsize=ann_fontsize, color="black", zorder=2)

ax.set_xlim(x_lo, x_hi)
ax.set_ylim(y_lo, y_hi)
ax.set_aspect("equal")
ax.tick_params(labelsize=16)
ax.set_xlabel("X (mm)", fontsize=18)
ax.set_ylabel("Y (mm)", fontsize=18)
ax.set_title(title, fontsize=16)

mean_val = np.mean(valid_vals)
std_val  = np.std(valid_vals)
ax.text(0.02, 0.97,
        f"Mean = {mean_val:.2f}%\nStd  = {std_val:.2f}%",
        transform=ax.transAxes, fontsize=14, va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)
cbar = plt.colorbar(sm, cax=cax)
cax.tick_params(labelsize=14)
cbar.set_label("Relative difference (%)", fontsize=16)

fig.savefig(output)
print(f"\nSaved: {output}")

# ── 6. Factor distribution histogram ─────────────────────────────────────────
ref_vals_hist  = [ref_map[m["n"]]  for m in w_modules if ref_map.get(m["n"], 0.0) != 0.0]
fadc_vals_hist = [fadc_map[m["n"]] for m in w_modules if fadc_map.get(m["n"], 0.0) != 0.0]

all_vals_hist = ref_vals_hist + fadc_vals_hist
bins = np.linspace(min(all_vals_hist), max(all_vals_hist), 60)

fig2, ax2 = plt.subplots(figsize=(9, 6), dpi=160)
ax2.hist(ref_vals_hist,  bins=bins, alpha=0.6, label='Physics', color='steelblue')
ax2.hist(fadc_vals_hist, bins=bins, alpha=0.6, label='FADC',    color='tomato')
ax2.axvline(0.122, color='black', linestyle='--', linewidth=1.8, label='Original factor (0.122)')
ax2.set_xlabel('Calibration Factor', fontsize=14)
ax2.set_ylabel('Counts', fontsize=14)
ax2.set_title('Calibration Factor Distribution', fontsize=15)
ax2.legend(fontsize=13)
ax2.tick_params(labelsize=12)
fig2.tight_layout()
output2 = 'hycal_factor_dist.png'
fig2.savefig(output2)
print(f"Saved: {output2}")

# ── 7. HyCal map: FADC factor vs. original (0.122), relative difference [%] ──
orig_factor = 0.122

diffs_fadc_orig = {}
for m in w_modules:
    name = m["n"]
    f = fadc_map.get(name, 0.0)
    if f == 0.0:
        diffs_fadc_orig[name] = None
    else:
        diffs_fadc_orig[name] = (f - orig_factor) / orig_factor * 100.0

valid_fo = [v for v in diffs_fadc_orig.values() if v is not None]
print(f"\nFADC vs Original diff range: {min(valid_fo):.3f}% to {max(valid_fo):.3f}%")
print(f"Mean: {np.mean(valid_fo):.3f}%   Std: {np.std(valid_fo):.3f}%")

_vmin3 = -30.0
_vmax3 =  30.0
norm3 = mcolors.Normalize(vmin=_vmin3, vmax=_vmax3)
sm3 = ScalarMappable(cmap=plt.get_cmap(cmap_name), norm=norm3)
sm3.set_array(valid_fo)

fig3, ax3 = plt.subplots(figsize=(16, 14), dpi=160,
                         gridspec_kw=dict(left=0.08, right=0.88, bottom=0.06, top=0.94))

for m in w_modules:
    name = m["n"]
    val  = diffs_fadc_orig[name]
    fc   = "lightgrey" if val is None else sm3.to_rgba(val)
    rx = m["x"] - m["sx"] / 2
    ry = m["y"] - m["sy"] / 2
    ax3.add_patch(Rectangle((rx, ry), m["sx"], m["sy"],
                             fc=fc, ec="black", lw=0.25, zorder=1))
    ax3.text(m["x"], m["y"], name[1:],
             ha="center", va="center",
             fontsize=ann_fontsize, color="black", zorder=2)

ax3.set_xlim(x_lo, x_hi)
ax3.set_ylim(y_lo, y_hi)
ax3.set_aspect("equal")
ax3.tick_params(labelsize=16)
ax3.set_xlabel("X (mm)", fontsize=18)
ax3.set_ylabel("Y (mm)", fontsize=18)
ax3.set_title(f'Calibration Factor Relative Difference: (FADC − Original) / Original  [%]  (Original = {orig_factor})', fontsize=15)

mean_fo = np.mean(valid_fo)
std_fo  = np.std(valid_fo)
ax3.text(0.02, 0.97,
         f"Mean = {mean_fo:.2f}%\nStd  = {std_fo:.2f}%",
         transform=ax3.transAxes, fontsize=14, va="top",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

divider3 = make_axes_locatable(ax3)
cax3 = divider3.append_axes("right", size="3%", pad=0.1)
cbar3 = plt.colorbar(sm3, cax=cax3)
cax3.tick_params(labelsize=14)
cbar3.set_label("Relative difference (%)", fontsize=16)

output3 = 'hycal_fadc_vs_orig.png'
fig3.savefig(output3)
print(f"Saved: {output3}")
