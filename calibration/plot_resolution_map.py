#!/usr/bin/env python3
"""
Re-fit module energy histograms from CalibResult_iter8.root for config0 and config2,
compute energy resolution = sigma / sqrt(ExpectedPeak), then:
  1) Plot 1-D resolution distributions for both configs on one figure
  2) Plot resolution maps for both configs with the same color range
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

# ── Paths ────────────────────────────────────────────────────────────────────
BASE   = os.path.dirname(os.path.abspath(__file__))
CALIB  = os.path.join(BASE, '..', 'data', 'calib', 'Physics_calib')
GEO    = os.path.join(BASE, '..', 'hycal_modules.json')
ITER   = 8

CONFIGS = {
    'config0': os.path.join(CALIB, 'config0_5by5'),
    'config2': os.path.join(CALIB, 'config2_5by5'),
}

# ── Plot control options ──────────────────────────────────────────────────────
# 1-D distribution
HIST1D_XMIN  = 2    # lower edge of x-axis; None = auto (2nd percentile)
HIST1D_XMAX  = 8    # upper edge of x-axis; None = auto (98th percentile)
HIST1D_BINS  = 40      # number of bins

# 2-D resolution map
MAP_VMIN     = 0    # colorbar min; None = auto (2nd percentile)
MAP_VMAX     = 6    # colorbar max; None = auto (98th percentile)
MAP_CMAP     = 'rainbow'   # matplotlib colormap name

# 2-D ratio map  (MeasuredPeak / ExpectedPeak)
RATIO_VMIN   = 0.99   # colorbar min; None = auto (2nd percentile)
RATIO_VMAX   = 1.01   # colorbar max; None = auto (98th percentile)
RATIO_CMAP   = 'RdBu_r'   # matplotlib colormap name
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Parse fitting_parameters_iterN.dat ────────────────────────────────────
def parse_dat(path):
    """Return dict: module_name -> {'expected': float, 'sigma': float}"""
    result = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('Module'):
                continue
            parts = line.split()
            # columns: Module  ExpectedPeak  MeasuredPeak  oldFactor  Ratio  Sigma  Chi2/ndf
            # W63 / W66 in config0 are missing oldFactor → fewer columns
            name = parts[0]
            try:
                expected = float(parts[1])
                # sigma is always the second-to-last column
                sigma    = float(parts[-2])
            except (ValueError, IndexError):
                # skip if ExpectedPeak is missing
                continue
            if expected > 0 and sigma > 0:
                result[name] = {'expected': expected, 'sigma': sigma}
    return result

# ── 2. Re-fit one histogram ───────────────────────────────────────────────────
def refit_gauss(hist, center, sigma_hint):
    """
    Fit a Gaussian to hist in [center - sigma_hint, center + sigma_hint].
    Returns (mean, sigma) or (None, None) on failure.
    """
    lo = center - sigma_hint
    hi = center + sigma_hint
    g  = ROOT.TF1('gfit', 'gaus', lo, hi)
    g.SetParameters(hist.GetMaximum(), center, sigma_hint)
    status = hist.Fit(g, 'RQNS')   # R=range, Q=quiet, N=no store, S=save result
    if status and status.IsValid():
        return g.GetParameter(1), abs(g.GetParameter(2))
    # Fallback: try with looser range (2 sigma)
    lo2 = center - 2 * sigma_hint
    hi2 = center + 2 * sigma_hint
    g2  = ROOT.TF1('gfit2', 'gaus', lo2, hi2)
    g2.SetParameters(hist.GetMaximum(), center, sigma_hint)
    status2 = hist.Fit(g2, 'RQNS')
    if status2 and status2.IsValid():
        return g2.GetParameter(1), abs(g2.GetParameter(2))
    return None, None

# ── 3. Process one config ─────────────────────────────────────────────────────
def process_config(config_dir):
    dat_path  = os.path.join(config_dir, f'fitting_parameters_iter{ITER}.dat')
    root_path = os.path.join(config_dir, f'CalibResult_iter{ITER}.root')

    params = parse_dat(dat_path)
    print(f"  Parsed {len(params)} modules from dat")

    rfile = ROOT.TFile(root_path, 'READ')
    if rfile.IsZombie():
        sys.exit(f"Cannot open {root_path}")
    d_energy = rfile.Get('module_energy')
    if not d_energy:
        sys.exit("TDirectory 'module_energy' not found")

    resolutions = {}   # module_name -> resolution value
    ratios      = {}   # module_name -> MeasuredPeak / ExpectedPeak
    for name, p in params.items():
        hist = d_energy.Get(f'h_{name}')
        if not hist:
            print(f"  WARNING: histogram h_{name} not found, skipping")
            continue
        mean, sigma = refit_gauss(hist, p['expected'], p['sigma'])
        if mean is None:
            print(f"  WARNING: fit failed for {name}, skipping")
            continue
        # relative resolution = sigma / sqrt(ExpectedPeak)
        resolutions[name] = 100 * sigma / p['expected'] * np.sqrt(p['expected']/1000)
        ratios[name]      = mean / p['expected']

    rfile.Close()
    print(f"  Computed resolution for {len(resolutions)} modules")
    return resolutions, ratios

# ── 4. Load geometry ──────────────────────────────────────────────────────────
def load_geometry(geo_path):
    with open(geo_path) as f:
        mods = json.load(f)
    return {m['n']: m for m in mods if m['t'] == 'PbWO4'}

# ── 5. Build outer-3-layer module set ───────────────────────────────────────
def build_outer_set(geo):
    """
    Return a set of module names that belong to the outermost 3 rings
    (i.e. whose x-column index or y-row index is within the first/last 3).
    """
    xs = sorted(set(round(m['x'], 1) for m in geo.values()))
    ys = sorted(set(round(m['y'], 1) for m in geo.values()))
    outer_x = set(xs[:3] + xs[-3:])
    outer_y = set(ys[:3] + ys[-3:])
    return {name for name, m in geo.items()
            if round(m['x'], 1) in outer_x or round(m['y'], 1) in outer_y}

def build_center_set(geo, n=6):
    """
    Return a set of module names in the central n×n region
    (the n/2 column-positions nearest x=0 AND n/2 row-positions nearest y=0).
    """
    xs = sorted(set(round(m['x'], 1) for m in geo.values()))
    ys = sorted(set(round(m['y'], 1) for m in geo.values()))
    half = n // 2
    mid  = len(xs) // 2
    center_x = set(xs[mid - half : mid + half])
    center_y = set(ys[mid - half : mid + half])
    return {name for name, m in geo.items()
            if round(m['x'], 1) in center_x and round(m['y'], 1) in center_y}

# ── 6. Main ───────────────────────────────────────────────────────────────────
def main():
    results  = {}   # cfg -> resolutions dict
    ratios   = {}   # cfg -> ratios dict
    for cfg, cdir in CONFIGS.items():
        print(f"Processing {cfg} ...")
        results[cfg], ratios[cfg] = process_config(cdir)

    geo = load_geometry(GEO)
    outer_set  = build_outer_set(geo)
    center_set = build_center_set(geo, n=6)
    exclude_1d = outer_set | center_set
    print(f"Outer 3-layer modules: {len(outer_set)},  center 6×6: {len(center_set)},  used in 1D: {len(geo) - len(exclude_1d)}")

    # ── collect all resolution values to set a common range ──
    all_res = [v for cfg in results.values() for v in cfg.values()]
    auto_lo = np.percentile(all_res, 2)
    auto_hi = np.percentile(all_res, 98)
    print(f"\nResolution range (2\u201398 pct): {auto_lo:.4f} \u2013 {auto_hi:.4f}")

    h1d_lo = HIST1D_XMIN if HIST1D_XMIN is not None else auto_lo
    h1d_hi = HIST1D_XMAX if HIST1D_XMAX is not None else auto_hi
    map_lo = MAP_VMIN    if MAP_VMIN    is not None else auto_lo
    map_hi = MAP_VMAX    if MAP_VMAX    is not None else auto_hi

    # ── Figure 1: 1-D distributions (inner modules only) ─────────────────────
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    colors = {'config0': 'steelblue', 'config2': 'tomato'}
    for cfg, res in results.items():
        vals = [v for name, v in res.items() if name not in exclude_1d]
        ax1.hist(vals, bins=HIST1D_BINS, range=(h1d_lo, h1d_hi),
                 histtype='step', linewidth=1.8,
                 color=colors[cfg], label=cfg)
    ax1.set_xlabel(r'$\sigma \,/\, \sqrt{E}$  (%)', fontsize=13)
    ax1.set_ylabel('Modules', fontsize=13)
    ax1.set_title('PbWO4 Energy Resolution Distribution (3.5 GeV)', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.xaxis.set_major_locator(ticker.AutoLocator())
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax1.tick_params(axis='x', which='major', labelsize=11, length=6)
    ax1.tick_params(axis='x', which='minor', length=3)
    ax1.tick_params(axis='y', labelsize=11)
    fig1.tight_layout()
    out1 = os.path.join(BASE, 'resolution_distribution.png')
    fig1.savefig(out1, dpi=150)
    print(f"Saved: {out1}")
    plt.close(fig1)

    # ── Figure 2: resolution maps ─────────────────────────────────────────────
    cmap = plt.get_cmap(MAP_CMAP)
    norm = mcolors.Normalize(vmin=map_lo, vmax=map_hi)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    fig2, axes = plt.subplots(1, 2, figsize=(22, 10),
                              gridspec_kw=dict(left=0.04, right=0.90,
                                               bottom=0.05, top=0.85,
                                               wspace=0.06))

    xs = [m['x'] for m in geo.values()]
    ys = [m['y'] for m in geo.values()]
    margin = 25
    x_lo, x_hi = min(xs) - margin, max(xs) + margin
    y_lo, y_hi = min(ys) - margin, max(ys) + margin

    cfg_list = list(CONFIGS.keys())
    for ax, cfg in zip(axes, cfg_list):
        res = results[cfg]
        # draw all W modules
        for name, m in geo.items():
            val = res.get(name, None)
            fc  = sm.to_rgba(val) if val is not None else 'lightgrey'
            rx  = m['x'] - m['sx'] / 2
            ry  = m['y'] - m['sy'] / 2
            ax.add_patch(Rectangle((rx, ry), m['sx'], m['sy'],
                                    fc=fc, ec='black', lw=0.15, zorder=1))
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect('equal')
        ax.set_title(cfg, fontsize=32)
        ax.tick_params(labelsize=22)

    # shared colorbar on the right
    cbar_ax = fig2.add_axes([0.91, 0.08, 0.018, 0.82])
    cb = fig2.colorbar(sm, cax=cbar_ax)
    cb.set_label(r'$\sigma \,/\, \sqrt{E}$  (%)', fontsize=26)
    cb.ax.tick_params(labelsize=22)

    fig2.suptitle('PbWO4 Energy Resolution Map  (3.5 GeV)', fontsize=34, y=0.97)
    out2 = os.path.join(BASE, 'resolution_map.png')
    fig2.savefig(out2, dpi=150)
    print(f"Saved: {out2}")
    plt.close(fig2)

    # ── Figure 3: E/E_expect ratio maps ──────────────────────────────────────
    all_ratio = [v for cfg in ratios.values() for v in cfg.values()]
    ratio_auto_lo = np.percentile(all_ratio, 2)
    ratio_auto_hi = np.percentile(all_ratio, 98)
    ratio_lo = RATIO_VMIN if RATIO_VMIN is not None else ratio_auto_lo
    ratio_hi = RATIO_VMAX if RATIO_VMAX is not None else ratio_auto_hi

    cmap3 = plt.get_cmap(RATIO_CMAP)
    norm3 = mcolors.Normalize(vmin=ratio_lo, vmax=ratio_hi)
    sm3   = ScalarMappable(cmap=cmap3, norm=norm3)
    sm3.set_array([])

    fig3, axes3 = plt.subplots(1, 2, figsize=(22, 10),
                               gridspec_kw=dict(left=0.04, right=0.90,
                                                bottom=0.05, top=0.85,
                                                wspace=0.06))
    for ax, cfg in zip(axes3, cfg_list):
        rat = ratios[cfg]
        for name, m in geo.items():
            val = rat.get(name, None)
            fc  = sm3.to_rgba(val) if val is not None else 'lightgrey'
            rx  = m['x'] - m['sx'] / 2
            ry  = m['y'] - m['sy'] / 2
            ax.add_patch(Rectangle((rx, ry), m['sx'], m['sy'],
                                    fc=fc, ec='black', lw=0.15, zorder=1))
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect('equal')
        ax.set_title(cfg, fontsize=32)
        ax.tick_params(labelsize=22)

    cbar_ax3 = fig3.add_axes([0.91, 0.08, 0.018, 0.82])
    cb3 = fig3.colorbar(sm3, cax=cbar_ax3)
    cb3.set_label(r'$E_\mathrm{meas}\,/\,E_\mathrm{expect}$', fontsize=26)
    cb3.ax.tick_params(labelsize=22)

    fig3.suptitle(r'PbWO4  $E_\mathrm{meas}\,/\,E_\mathrm{expect}$  Map  (3.5 GeV)', fontsize=34, y=0.97)
    out3 = os.path.join(BASE, 'ratio_map.png')
    fig3.savefig(out3, dpi=150)
    print(f"Saved: {out3}")
    plt.close(fig3)


if __name__ == '__main__':
    main()
