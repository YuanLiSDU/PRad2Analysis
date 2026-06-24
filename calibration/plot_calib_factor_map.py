#!/usr/bin/env python3
"""
Plot PbWO4 calibration factors from one calibration JSON file.

Outputs:
  - calib_factor_map.png: HyCal W-module map colored by calibration factor
  - calib_factor_hist.png: 1-D factor distribution with mean marked
"""

import argparse
import json
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


BASE = Path(__file__).resolve().parent

DEFAULT_CALIB = BASE / "calib_iter11.json"
DEFAULT_GEO = BASE.parent / "hycal_modules.json"

MAP_OUTPUT = BASE / "calib_factor_map.png"
HIST_OUTPUT = BASE / "calib_factor_hist.png"
OUTLIER_OUTPUT = BASE / "calib_factor_outliers.txt"

HYCal_RAINBOW_STOPS = [
    (0.00, (30 / 255,  58 / 255,  95 / 255)),
    (0.25, (59 / 255, 130 / 255, 246 / 255)),
    (0.50, (45 / 255, 212 / 255, 160 / 255)),
    (0.75, (234 / 255, 179 / 255,   8 / 255)),
    (1.00, (245 / 255, 101 / 255, 101 / 255)),
]
ANN_FONTSIZE = 7.5
MAP_VMIN = None
MAP_VMAX = 0.2
LOW_FACTOR_CUT = 0.11
HIGH_FACTOR_CUT = 0.185


def load_calib_factor(calib_file):
    with open(calib_file) as f:
        entries = json.load(f)
    return {entry["name"]: entry["factor"] for entry in entries}


def load_w_modules(geo_file):
    with open(geo_file) as f:
        modules = json.load(f)
    return [m for m in modules if m["t"] == "PbWO4"]


def collect_w_factors(w_modules, factor_map):
    factors = {}
    missing = []
    invalid = []

    for m in w_modules:
        name = m["n"]
        factor = factor_map.get(name)
        if factor is None:
            missing.append(name)
            factors[name] = None
        elif factor <= 0.0:
            invalid.append(name)
            factors[name] = None
        else:
            factors[name] = factor

    valid_vals = np.array([v for v in factors.values() if v is not None], dtype=float)
    return factors, valid_vals, missing, invalid


def plot_factor_map(w_modules, factors, valid_vals, output):
    xs = [m["x"] for m in w_modules]
    ys = [m["y"] for m in w_modules]
    margin = 30

    vmin = np.min(valid_vals) if MAP_VMIN is None else MAP_VMIN
    vmax = np.max(valid_vals) if MAP_VMAX is None else MAP_VMAX
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap = mcolors.LinearSegmentedColormap.from_list("hycal_rainbow", HYCal_RAINBOW_STOPS)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(valid_vals)

    fig, ax = plt.subplots(
        figsize=(16, 14),
        dpi=160,
        gridspec_kw=dict(left=0.08, right=0.88, bottom=0.06, top=0.94),
    )

    for m in w_modules:
        name = m["n"]
        factor = factors[name]
        fc = "lightgrey" if factor is None else sm.to_rgba(factor)

        rx = m["x"] - m["sx"] / 2
        ry = m["y"] - m["sy"] / 2
        ax.add_patch(Rectangle((rx, ry), m["sx"], m["sy"],
                               fc=fc, ec="black", lw=0.25, zorder=1))
        ax.text(m["x"], m["y"], name[1:],
                ha="center", va="center",
                fontsize=ANN_FONTSIZE, color="black", zorder=2)

    mean = np.mean(valid_vals)
    std = np.std(valid_vals)
    ax.text(
        0.02,
        0.97,
        f"Mean = {mean:.6f}\nStd  = {std:.6f}\nN    = {len(valid_vals)}",
        transform=ax.transAxes,
        fontsize=14,
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.75),
    )

    ax.set_xlim(min(xs) - margin, max(xs) + margin)
    ax.set_ylim(min(ys) - margin, max(ys) + margin)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)", fontsize=18)
    ax.set_ylabel("Y (mm)", fontsize=18)
    ax.set_title("PbWO4 Calibration Factor Map", fontsize=18)
    ax.tick_params(labelsize=16)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label("Calibration factor", fontsize=16)
    cax.tick_params(labelsize=14)

    fig.savefig(output)
    plt.close(fig)


def plot_factor_hist(valid_vals, output):
    mean = np.mean(valid_vals)
    std = np.std(valid_vals)

    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
    ax.hist(valid_vals, bins=70, histtype="stepfilled",
            color="steelblue", alpha=0.65, edgecolor="black", linewidth=0.8)
    ax.axvline(mean, color="crimson", linestyle="--", linewidth=2.0,
               label=f"Mean = {mean:.6f}")

    ax.set_xlabel("Calibration factor", fontsize=14)
    ax.set_ylabel("Counts", fontsize=14)
    ax.set_title("PbWO4 Calibration Factor Distribution", fontsize=15)
    ax.text(
        0.97,
        0.95,
        f"Mean = {mean:.6f}\nStd  = {std:.6f}\nN    = {len(valid_vals)}",
        transform=ax.transAxes,
        fontsize=12,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.75),
    )
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=12)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def write_outliers(w_modules, factors, output):
    outliers = []
    for m in w_modules:
        name = m["n"]
        factor = factors.get(name)
        if factor is None:
            continue
        if factor < LOW_FACTOR_CUT or factor > HIGH_FACTOR_CUT:
            outliers.append((name, factor))

    outliers.sort(key=lambda item: item[1])

    with open(output, "w") as f:
        f.write("module factor\n")
        for name, factor in outliers:
            f.write(f"{name} {factor:.9f}\n")

    print(f"\nModules with factor < {LOW_FACTOR_CUT} or > {HIGH_FACTOR_CUT}: {len(outliers)}")
    print("module factor")
    for name, factor in outliers:
        print(f"{name} {factor:.9f}")


def main():
    parser = argparse.ArgumentParser(description="Plot PbWO4 calibration factor map and histogram.")
    parser.add_argument("calib_file", nargs="?", default=DEFAULT_CALIB,
                        help="Calibration JSON file. Default: calibration/calib_iter11.json")
    parser.add_argument("--geo", default=DEFAULT_GEO,
                        help="HyCal geometry JSON file. Default: hycal_modules.json")
    parser.add_argument("--map-output", default=MAP_OUTPUT,
                        help="Output PNG for factor map.")
    parser.add_argument("--hist-output", default=HIST_OUTPUT,
                        help="Output PNG for 1-D histogram.")
    parser.add_argument("--outlier-output", default=OUTLIER_OUTPUT,
                        help="Output text file for modules with factor < 0.11 or > 0.185.")
    args = parser.parse_args()

    calib_file = Path(args.calib_file)
    geo_file = Path(args.geo)
    map_output = Path(args.map_output)
    hist_output = Path(args.hist_output)
    outlier_output = Path(args.outlier_output)

    factor_map = load_calib_factor(calib_file)
    w_modules = load_w_modules(geo_file)
    factors, valid_vals, missing, invalid = collect_w_factors(w_modules, factor_map)

    if len(valid_vals) == 0:
        raise RuntimeError("No valid PbWO4 calibration factors found.")

    print(f"Calibration file: {calib_file}")
    print(f"W modules in geometry: {len(w_modules)}")
    print(f"Valid W factors:       {len(valid_vals)}")
    print(f"Missing in calib JSON: {len(missing)}")
    print(f"Invalid/zero factors:  {len(invalid)}")
    print(f"Mean: {np.mean(valid_vals):.6f}   Std: {np.std(valid_vals):.6f}")
    print(f"Range: {np.min(valid_vals):.6f} to {np.max(valid_vals):.6f}")

    plot_factor_map(w_modules, factors, valid_vals, map_output)
    plot_factor_hist(valid_vals, hist_output)
    write_outliers(w_modules, factors, outlier_output)

    print(f"Saved: {map_output}")
    print(f"Saved: {hist_output}")
    print(f"Saved: {outlier_output}")


if __name__ == "__main__":
    main()
