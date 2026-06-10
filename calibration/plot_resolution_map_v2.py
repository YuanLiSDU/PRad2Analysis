#!/usr/bin/env python3
"""
Plot PbWO4 energy resolution map from a recon ROOT file.

Algorithm:
  - Read 'recon' tree, single-cluster events (n_clusters == 1, cl_nblocks >= 3)
  - Assign each hit to a module by checking which module's bounding box contains
    the hit centre (cl_x, cl_y)
  - Fill a per-module energy histogram
  - Fit a Gaussian around the expected ep-elastic peak
  - Resolution = sigma / E_peak * sqrt(E_peak / 1000) * 100  [%]
  - Show only modules in the "middle ring":
      - at least 3 layers away from the centre 2×2 hole
      - at least 3 layers away from the outer edge

Usage (must source setup_env.sh first):
  python3 plot_resolution_map_v2.py [recon_file.root]

Defaults to ../data/calib/prad_024512_recon.root if no argument given.
"""

import os
import sys
import json
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

# ── Configuration ─────────────────────────────────────────────────────────────
BASE      = os.path.dirname(os.path.abspath(__file__))
GEO_FILE  = os.path.join(BASE, '..', 'hycal_modules.json')

DEFAULT_ROOT_FILE = os.path.join(BASE, '..', 'data', 'recon', '3.5GeV',
                                 'prad_024917_recon.root')

M_PROTON  = 938.272   # MeV
SHIFT_X   = 0.7       # mm  (position alignment, same as energy_plot.C)
SHIFT_Y   = 1.83      # mm
MIN_NBLOCKS = 3       # minimum blocks per cluster

# Histogram range and bins for each module's energy distribution
HIST_NBINS  = 420
HIST_EMIN   = 0.0   # MeV
HIST_EMAX   = 4200.0  # MeV

# Only use the central fraction of each module for the resolution fit
# (|xd| < CENTER_FRAC and |yd| < CENTER_FRAC, where xd = (x-cx)/sx)
CENTER_FRAC = 0.3

# Fit window: ± FIT_NSIGMA × sigma_hint around expected peak
FIT_NSIGMA  = 1

# 1D resolution plot range
HIST1D_XMIN = 2.0    # %
HIST1D_XMAX = 5.0    # %
HIST1D_BINS = 30

# 2D map colour range
MAP_VMIN    = 2.0    # %
MAP_VMAX    = 8.0    # %
MAP_CMAP    = 'rainbow'
# ─────────────────────────────────────────────────────────────────────────────


# ── 1. Geometry ───────────────────────────────────────────────────────────────

def load_geometry(geo_path):
    """Load PbWO4 module geometry. Returns dict name -> module dict."""
    with open(geo_path) as f:
        mods = json.load(f)
    return {m['n']: m for m in mods if m['t'] == 'PbWO4'}


def build_grid_layers(geo):
    """
    Assign each module a (col_idx, row_idx) in the 34×34 grid and compute:
      layer_from_hole  – Chebyshev distance from nearest hole cell
                         (hole sits at column/row indices 16 and 17)
      layer_from_edge  – 1-indexed distance from the outer border
                         (edge modules = 1)

    Returns dict name -> {'ci': int, 'ri': int,
                          'layer_hole': int, 'layer_edge': int}
    """
    xs = sorted(set(round(m['x'], 1) for m in geo.values()))
    ys = sorted(set(round(m['y'], 1) for m in geo.values()))
    x_idx = {x: i for i, x in enumerate(xs)}
    y_idx = {y: i for i, y in enumerate(ys)}
    N = len(xs)  # 34
    M = len(ys)  # 34

    # The centre hole occupies col indices 16 & 17 (and row indices 16 & 17).
    # col distance from the hole for column index ci:
    #   ci <= 15 → 16 - ci   (adjacent to hole = 1)
    #   ci >= 18 → ci - 17
    # (indices 16 & 17 have no modules since they are the hole)
    def col_dist(ci):
        if ci <= 15:
            return 16 - ci
        else:                  # ci >= 18
            return ci - 17

    def row_dist(ri):
        if ri <= 15:
            return 16 - ri
        else:
            return ri - 17

    result = {}
    for name, m in geo.items():
        ci = x_idx[round(m['x'], 1)]
        ri = y_idx[round(m['y'], 1)]
        layer_hole = max(col_dist(ci), row_dist(ri))
        layer_edge = min(ci, N - 1 - ci, ri, M - 1 - ri) + 1
        result[name] = {'ci': ci, 'ri': ri,
                        'layer_hole': layer_hole,
                        'layer_edge': layer_edge}
    return result


def build_module_lookup(geo):
    """
    Build a fast lookup: for a given (x, y) return the module name or None.
    Uses a dict keyed by (col_idx, row_idx) for O(1) lookup after
    snapping the hit to the nearest grid position.
    """
    xs = sorted(set(round(m['x'], 1) for m in geo.values()))
    ys = sorted(set(round(m['y'], 1) for m in geo.values()))
    xs_arr = np.array(xs)
    ys_arr = np.array(ys)

    # Map (x_val, y_val) -> name for quick lookup
    pos_to_name = {(round(m['x'], 1), round(m['y'], 1)): name
                   for name, m in geo.items()}

    def find_module(hx, hy):
        # Snap to nearest grid column and row
        ci = int(np.argmin(np.abs(xs_arr - hx)))
        ri = int(np.argmin(np.abs(ys_arr - hy)))
        key = (xs_arr[ci], ys_arr[ri])
        return pos_to_name.get(key, None)

    return find_module


# ── 2. Expected energy (ep elastic) ──────────────────────────────────────────

def expected_energy_ep(hx, hy, hz, ebeam):
    """
    Compute expected energy of ep-elastically-scattered electron.
    hx, hy, hz: cluster position in mm; ebeam: beam energy in MeV.
    Returns MeV, or 0 if angle is not physical.
    """
    r2 = hx * hx + hy * hy
    if hz <= 0 or r2 <= 0:
        return 0.0
    cos_t = hz / math.sqrt(r2 + hz * hz)
    denom = M_PROTON + ebeam * (1.0 - cos_t)
    if denom <= 0:
        return 0.0
    return ebeam * M_PROTON / denom


# ── 3. Fill histograms ────────────────────────────────────────────────────────

def fill_histograms(root_path, geo, find_module, ebeam_override=None, center_only=True):
    """
    Loop over the 'recon' tree in root_path.
    Single-cluster events (n_clusters == 1, cl_nblocks[0] >= MIN_NBLOCKS):
      - apply position shift
      - find module by hit centre
      - if center_only=True, only use hits with |xd|<CENTER_FRAC, |yd|<CENTER_FRAC
      - fill per-module TH1F and per-module E_expected accumulator
    Returns:
      hists      – dict name -> TH1F
      eexp_sum   – dict name -> (sum_Eexp, count)  for average expected energy
    """
    rfile = ROOT.TFile.Open(root_path, 'READ')
    if not rfile or rfile.IsZombie():
        sys.exit(f'ERROR: cannot open {root_path}')
    tree = rfile.Get('recon')
    if not tree:
        sys.exit("ERROR: 'recon' tree not found")

    # Pre-create histograms
    hists    = {}
    eexp_acc = {}   # (sum, count)
    for name in geo:
        hists[name]    = ROOT.TH1F(f'h_{name}', name,
                                   HIST_NBINS, HIST_EMIN, HIST_EMAX)
        hists[name].SetDirectory(0)   # detach from file
        eexp_acc[name] = [0.0, 0]

    # Check whether EBeam branch exists in this file
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    has_ebeam = 'EBeam' in branch_names
    if not has_ebeam and ebeam_override is None:
        ebeam_override = 3485.43   # MeV — default for 3.5 GeV runs
        print(f'EBeam branch not found; using default {ebeam_override} MeV')

    n_entries = tree.GetEntries()
    print(f'Tree entries: {n_entries:,}')
    for i, ev in enumerate(tree):
        if i % 100_000 == 0:
            print(f'  {i:,} / {n_entries:,}\r', end='', flush=True)
        if ev.n_clusters != 1:
            continue
        if ev.cl_nblocks[0] < MIN_NBLOCKS:
            continue

        hx = ev.cl_x[0] + SHIFT_X
        hy = ev.cl_y[0] + SHIFT_Y
        hz = ev.cl_z[0]
        en = ev.cl_energy[0]

        name = find_module(hx, hy)
        if name is None:
            continue

        # Centre-region filter: only use hits near the module centre
        m = geo[name]
        xd = (hx - m['x']) / m['sx']
        yd = (hy - m['y']) / m['sy']
        if center_only and (abs(xd) >= CENTER_FRAC or abs(yd) >= CENTER_FRAC):
            continue

        if ebeam_override is not None:
            ebeam = ebeam_override
        else:
            ebeam = float(ev.EBeam) if ev.EBeam > 100 else 3485.43
        eexp  = expected_energy_ep(hx, hy, hz, ebeam)

        hists[name].Fill(en)
        if eexp > 0:
            eexp_acc[name][0] += eexp
            eexp_acc[name][1] += 1

    print()   # newline after progress
    rfile.Close()
    return hists, eexp_acc


# ── 4. Fit and compute resolution ─────────────────────────────────────────────

def gauss_fit(hist, center, sigma_hint):
    """
    Fit a Gaussian to hist in [center - FIT_NSIGMA*sigma_hint,
                                center + FIT_NSIGMA*sigma_hint].
    Returns (mean, sigma) or (None, None) on failure.
    """
    lo = center - FIT_NSIGMA * sigma_hint
    hi = center + FIT_NSIGMA * sigma_hint
    if hi <= lo or hist.GetEntries() < 20:
        return None, None
    g = ROOT.TF1('_gfit', 'gaus', lo, hi)
    g.SetParameters(hist.GetMaximum(), center, sigma_hint)
    r = hist.Fit(g, 'RQNS')
    if r and r.IsValid() and g.GetParameter(2) > 0:
        return g.GetParameter(1), abs(g.GetParameter(2))
    # Fallback with looser window
    lo2, hi2 = center - 1 * sigma_hint, center + 1 * sigma_hint
    g2 = ROOT.TF1('_gfit2', 'gaus', lo2, hi2)
    g2.SetParameters(hist.GetMaximum(), center, sigma_hint)
    r2 = hist.Fit(g2, 'RQNS')
    if r2 and r2.IsValid() and g2.GetParameter(2) > 0:
        return g2.GetParameter(1), abs(g2.GetParameter(2))
    return None, None


def compute_resolutions(hists, eexp_acc):
    """
    For each module, fit histogram, compute resolution.
    Returns dict name -> resolution (%)
    """
    resolutions = {}
    for name, h in hists.items():
        acc = eexp_acc[name]
        if acc[1] < 20:
            continue
        eexp = acc[0] / acc[1]    # average expected energy (MeV)
        # sigma hint: assume ~3% * sqrt(E/GeV) resolution as starting point
        sigma_hint = 0.03 * eexp * math.sqrt(eexp / 1000.0)
        sigma_hint = max(sigma_hint, 20.0)   # at least 20 MeV

        mean, sigma = gauss_fit(h, eexp, sigma_hint)
        if mean is None or mean <= 0:
            continue
        # resolution = (sigma/E) * sqrt(E/GeV) * 100 %
        resolutions[name] = 100.0 * sigma / mean * math.sqrt(mean / 1000.0)
    return resolutions


# ── 5. Layer filtering ────────────────────────────────────────────────────────

def get_valid_modules(layers, min_hole_layer=2, min_edge_layer=6):
    """
    Return set of module names satisfying the layer constraints:
      layer_from_hole >= min_hole_layer  AND  layer_from_edge >= min_edge_layer
    """
    return {name for name, lyr in layers.items()
            if lyr['layer_hole'] >= min_hole_layer
            and lyr['layer_edge'] >= min_edge_layer}


# ── 6. Plotting ───────────────────────────────────────────────────────────────

def plot_1d(resolutions_center, resolutions_all, valid_modules, out_path):
    vals_c = [v for name, v in resolutions_center.items() if name in valid_modules]
    vals_a = [v for name, v in resolutions_all.items()    if name in valid_modules]
    if not vals_c and not vals_a:
        print('WARNING: no resolution values in valid region for 1D plot')
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    hist_kw = dict(bins=HIST1D_BINS, range=(HIST1D_XMIN, HIST1D_XMAX),
                   histtype='step', linewidth=2.0)
    if vals_c:
        mu_c  = np.mean(vals_c)
        med_c = np.median(vals_c)
        ax.hist(vals_c, label=f'center |xd|<{CENTER_FRAC}  mean={mu_c:.2f}%  median={med_c:.2f}%',
                color='steelblue', **hist_kw)
        ax.axvline(mu_c,  ls='--', color='steelblue', lw=1.2, alpha=0.7)
    if vals_a:
        mu_a  = np.mean(vals_a)
        med_a = np.median(vals_a)
        ax.hist(vals_a, label=f'all hits            mean={mu_a:.2f}%  median={med_a:.2f}%',
                color='tomato', **hist_kw)
        ax.axvline(mu_a,  ls='--', color='tomato', lw=1.2, alpha=0.7)
    ax.set_xlabel(r'$\sigma\,/\,\sqrt{E}\;\;(\%)$', fontsize=13)
    ax.set_ylabel('Modules', fontsize=13)
    ax.set_title('PbWO4 Energy Resolution (middle ring)', fontsize=14)
    ax.legend(fontsize=11)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.tick_params(axis='x', which='major', labelsize=11, length=6)
    ax.tick_params(axis='x', which='minor', length=3)
    ax.tick_params(axis='y', labelsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f'Saved: {out_path}')
    plt.close(fig)


def plot_2d_map(resolutions, geo, valid_modules, out_path):
    vals = [v for name, v in resolutions.items() if name in valid_modules]
    vmin = MAP_VMIN if MAP_VMIN is not None else (np.percentile(vals, 2) if vals else 0)
    vmax = MAP_VMAX if MAP_VMAX is not None else (np.percentile(vals, 98) if vals else 10)

    cmap = plt.get_cmap(MAP_CMAP)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    xs = [m['x'] for m in geo.values()]
    ys = [m['y'] for m in geo.values()]
    margin = 25
    x_lo, x_hi = min(xs) - margin, max(xs) + margin
    y_lo, y_hi = min(ys) - margin, max(ys) + margin

    fig, ax = plt.subplots(1, 1, figsize=(11, 10))
    for name, m in geo.items():
        val = resolutions.get(name, None) if name in valid_modules else None
        if val is not None:
            fc = sm.to_rgba(val)
        elif name in valid_modules:
            fc = 'lightgrey'   # fit failed
        else:
            fc = 'white'       # excluded (inner / outer layers)
        rx = m['x'] - m['sx'] / 2
        ry = m['y'] - m['sy'] / 2
        ax.add_patch(Rectangle((rx, ry), m['sx'], m['sy'],
                                fc=fc, ec='gray', lw=0.2, zorder=1))

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect('equal')
    ax.set_xlabel('X (mm)', fontsize=13)
    ax.set_ylabel('Y (mm)', fontsize=13)
    ax.tick_params(labelsize=11)

    divider_ax = fig.add_axes([0.88, 0.08, 0.025, 0.82])
    cb = fig.colorbar(sm, cax=divider_ax)
    cb.set_label(r'$\sigma\,/\,\sqrt{E}\;\;(\%)$', fontsize=14)
    cb.ax.tick_params(labelsize=12)

    fig.suptitle('PbWO4 Energy Resolution Map  (middle ring)', fontsize=16, y=0.98)
    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.07, top=0.93)
    fig.savefig(out_path, dpi=150)
    print(f'Saved: {out_path}')
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='PbWO4 resolution map from recon ROOT file')
    parser.add_argument('root_file', nargs='?', default=DEFAULT_ROOT_FILE,
                        help='Path to recon ROOT file (default: prad_024512_recon.root)')
    parser.add_argument('--ebeam', type=float, default=None,
                        help='Override beam energy in MeV (default: use EBeam branch)')
    parser.add_argument('--min-hole-layer', type=int, default=2,
                        help='Minimum layer from centre hole (default: 2)')
    parser.add_argument('--min-edge-layer', type=int, default=6,
                        help='Minimum layer from outer edge (default: 6)')
    args = parser.parse_args()

    root_path = os.path.abspath(args.root_file)
    if not os.path.isfile(root_path):
        sys.exit(f'ERROR: file not found: {root_path}')

    print(f'Input file : {root_path}')
    print(f'Geometry   : {GEO_FILE}')
    print(f'Layer filter: hole >= {args.min_hole_layer}, edge >= {args.min_edge_layer}')

    # 1. Geometry
    geo = load_geometry(GEO_FILE)
    print(f'Loaded {len(geo)} PbWO4 modules')

    layers       = build_grid_layers(geo)
    find_module  = build_module_lookup(geo)
    valid_mods   = get_valid_modules(layers,
                                     args.min_hole_layer, args.min_edge_layer)
    print(f'Valid modules (middle ring): {len(valid_mods)} / {len(geo)}')

    # 2. Fill histograms (centre region and all hits)
    print('--- Filling centre-region histograms ---')
    hists_c, eexp_acc_c = fill_histograms(root_path, geo, find_module, args.ebeam,
                                           center_only=True)
    print('--- Filling all-hits histograms ---')
    hists_a, eexp_acc_a = fill_histograms(root_path, geo, find_module, args.ebeam,
                                           center_only=False)

    # 3. Fit
    resolutions_c = compute_resolutions(hists_c, eexp_acc_c)
    resolutions_a = compute_resolutions(hists_a, eexp_acc_a)
    print(f'Fitted (centre): {len(resolutions_c)} modules')
    print(f'Fitted (all)   : {len(resolutions_a)} modules')
    for label, res in [('centre', resolutions_c), ('all', resolutions_a)]:
        valid_res = {n: v for n, v in res.items() if n in valid_mods}
        if valid_res:
            vals = list(valid_res.values())
            print(f'  [{label}] mean={np.mean(vals):.2f}%  '
                  f'median={np.median(vals):.2f}%  '
                  f'range=[{np.min(vals):.2f}, {np.max(vals):.2f}]%')

    # 4. Output paths (alongside this script)
    out_1d  = os.path.join(BASE, 'resolution_distribution_v2.png')
    out_map = os.path.join(BASE, 'resolution_map_v2.png')

    plot_1d(resolutions_c, resolutions_a, valid_mods, out_1d)
    plot_2d_map(resolutions_c, geo, valid_mods, out_map)


if __name__ == '__main__':
    main()
