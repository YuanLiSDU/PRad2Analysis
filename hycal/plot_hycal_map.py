#!/usr/bin/env python3
"""
Plot HyCal detector module map with signal integral values.
"""

import os
import json
import codecs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ── Parameters ───────────────────────────────────────────────
value_file = 'test_23232IntegralADC.txt'
geo_file   = 'hycal_modules.json'
output     = 'hycal_map_23232.png'
cmap_name  = 'Spectral_r'       # reversed colormap
use_log    = False               # log scale for color mapping
vmin       = None                # None = auto
vmax       = None                # None = auto
annotate   = False               # annotate module names
ann_color  = 'k'
ann_fs_w   = 4                   # font size for PbWO4
ann_fs_g   = 7                   # font size for PbGlass
title      = 'HyCal ADC Signal Integrals'
# ─────────────────────────────────────────────────────────────


def read_adc_file(path):
    """Read ADC integral file, handling UTF-16 or UTF-8."""
    for enc in ['utf-16', 'utf-8']:
        try:
            with codecs.open(path, 'r', encoding=enc) as f:
                lines = f.readlines()
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    else:
        raise RuntimeError(f"Cannot decode {path}")

    values = {}
    for line in lines[1:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                values[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return values


# Load data
with open(geo_file) as f:
    modules = json.load(f)
values = read_adc_file(value_file) if value_file else {}

cmap = plt.get_cmap(cmap_name)

# Plot bounds
min_x = min(m['x'] - m['sx'] / 2 for m in modules)
max_x = max(m['x'] + m['sx'] / 2 for m in modules)
min_y = min(m['y'] - m['sy'] / 2 for m in modules)
max_y = max(m['y'] + m['sy'] / 2 for m in modules)

fig, ax = plt.subplots(figsize=(18, 16), dpi=160,
                       gridspec_kw=dict(left=0.10, right=0.90, bottom=0.06, top=0.95))

# Color normalization
sm = None
if values:
    matched = [values[m['n']] for m in modules if m['n'] in values]
    if matched:
        lo = vmin if vmin is not None else min(matched)
        hi = vmax if vmax is not None else max(matched)
        norm = LogNorm(vmin=max(lo, 1), vmax=hi) if use_log else plt.Normalize(vmin=lo, vmax=hi)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array(matched)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        cbar = plt.colorbar(sm, cax=cax)
        cax.tick_params(labelsize=20)
        cbar.set_label('ADC Integral', fontsize=22)

# Draw modules
for m in modules:
    name = m['n']
    ax0 = m['x'] - m['sx'] / 2
    ay0 = m['y'] - m['sy'] / 2
    fc = sm.to_rgba(values[name]) if (sm and name in values) else 'linen'
    ax.add_patch(Rectangle((ax0, ay0), m['sx'], m['sy'], fc=fc, ec='black', lw=0.4))

    if annotate:
        label = name.lstrip('WG')
        fs = ann_fs_w if name.startswith('W') else ann_fs_g
        ax.annotate(label, (m['x'], m['y']), color=ann_color, fontsize=fs,
                    ha='center', va='center')

ax.set_xlim(min_x, max_x)
ax.set_ylim(min_y, max_y)
ax.set_aspect('equal')
ax.tick_params(labelsize=20)
ax.set_xlabel('X (mm)', fontsize=22)
ax.set_ylabel('Y (mm)', fontsize=22)
if title:
    ax.set_title(title, fontsize=26)

fig.savefig(output)
print(f"Saved: {output}  ({len(modules)} modules, {len(values)} values, "
      f"{sum(1 for m in modules if m['n'] in values)} matched)")
