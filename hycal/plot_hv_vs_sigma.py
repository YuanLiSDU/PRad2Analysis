#!/usr/bin/env python3
"""
Plot HV (V0Set) vs fitting sigma for W modules with |ratio - 1| <= RATIO_TOL,
excluding the outer EXCLUDE_LAYERS layers (keeping inner modules only).
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import uproot

FIT_FILE       = "fitting_parameters_iter5.dat"
HV_FILE        = "hv_settings_20260509_220811.json"
ROOT_FILE      = "test.root"
MODULE_FILE    = "../hycal_modules.json"
CALIB_FILE     = "calib_iter4.json"
RATIO_TOL      = 0.1
EXCLUDE_LAYERS = 4   # remove outer 4 layers, keep inner modules
DV_RMS_MAX     = 20  # V — remove obvious dV RMS outliers (gap between ~13 V cluster and 152/300 V outliers)


def load_fit_params(path):
    """Return dict: module_name -> {"ratio": float, "sigma": float}"""
    data = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            name   = parts[0]
            ratio  = float(parts[4])
            sigma  = float(parts[5])
            data[name] = {"ratio": ratio, "sigma": sigma}
    return data


def load_hv(path):
    """Return dict: module_name -> V0Set (float)"""
    with open(path) as f:
        doc = json.load(f)
    hv = {}
    for ch in doc["channels"]:
        name = ch.get("name", "")
        v = ch.get("params", {}).get("V0Set")
        if name and v is not None:
            hv[name] = float(v)
    return hv


def get_inner_modules(module_file, exclude_layers):
    """Return set of W module names excluding the outer exclude_layers layers."""
    with open(module_file) as f:
        mods = json.load(f)
    w_mods = [m for m in mods if m["n"].startswith("W")]
    xs = sorted(set(round(m["x"], 2) for m in w_mods))
    ys = sorted(set(round(m["y"], 2) for m in w_mods))
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    N = len(xs)  # 34
    result = set()
    for m in w_mods:
        ci = xi[round(m["x"], 2)]
        ri = yi[round(m["y"], 2)]
        layer = min(ci, N - 1 - ci, ri, N - 1 - ri) + 1
        if layer > exclude_layers:
            result.add(m["n"])
    return result


def load_dv_rms(root_file):
    """
    Read hv and hv_channels trees from ROOT file.
    Return dict: channel_name -> RMS of dv across all entries.
    """
    with uproot.open(root_file) as f:
        hvc  = f["hv_channels"]
        names = hvc.arrays(["name"], library="np")["name"]  # shape (N_ch,)

        hv   = f["hv"]
        dv_all = hv.arrays(["dv"], library="np")["dv"]      # shape (N_entries, N_ch)

    dv_matrix = np.stack(dv_all)          # (N_entries, N_ch)
    rms = np.std(dv_matrix, axis=0)       # sqrt(<dv²> - <dv>²), numerically stable
    return {name: rms[i] for i, name in enumerate(names)}


def load_calib_factors(path):
    """Return dict: module_name -> calib factor (float)"""
    with open(path) as f:
        doc = json.load(f)
    return {entry["name"]: entry["factor"] for entry in doc
            if entry.get("factor", 0.0) != 0.0}


def load_module_positions(module_file):
    """Return dict: module_name -> (x, y) for W modules."""
    with open(module_file) as f:
        mods = json.load(f)
    return {m["n"]: (m["x"], m["y"]) for m in mods if m["n"].startswith("W")}


def plot_2d_map(ax, mod_positions, values, title, cmap="RdYlBu_r"):
    """
    Plot a 2D spatial map of W modules on a regular grid.
    values: dict {module_name: float}  — only modules present are coloured.
    """
    all_x = sorted(set(round(x, 2) for x, _ in mod_positions.values()))
    all_y = sorted(set(round(y, 2) for _, y in mod_positions.values()))
    xi = {v: i for i, v in enumerate(all_x)}
    yi = {v: i for i, v in enumerate(all_y)}

    grid = np.full((len(all_y), len(all_x)), np.nan)
    for name, (x, y) in mod_positions.items():
        if name in values:
            grid[yi[round(y, 2)], xi[round(x, 2)]] = values[name]

    # Build cell edges for pcolormesh
    def edges(centers):
        c = np.array(centers, dtype=float)
        step = np.diff(c)
        e = np.empty(len(c) + 1)
        e[1:-1] = (c[:-1] + c[1:]) / 2
        e[0]  = c[0]  - step[0]  / 2
        e[-1] = c[-1] + step[-1] / 2
        return e

    xe = edges(all_x)
    ye = edges(all_y)
    XX, YY = np.meshgrid(xe, ye)

    masked = np.ma.masked_invalid(grid)
    pcm = ax.pcolormesh(XX, YY, masked, cmap=cmap, shading="flat")
    plt.colorbar(pcm, ax=ax)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel("y (mm)", fontsize=11)
    ax.set_title(title, fontsize=12)


def main():
    fit         = load_fit_params(FIT_FILE)
    hv          = load_hv(HV_FILE)
    dv_rms      = load_dv_rms(ROOT_FILE)
    calib       = load_calib_factors(CALIB_FILE)
    inner_mods  = get_inner_modules(MODULE_FILE, EXCLUDE_LAYERS)
    mod_pos     = load_module_positions(MODULE_FILE)

    hv_vals    = []
    sigma_vals = []
    sigma_map  = {}
    hv_map     = {}
    dv_rms_vals   = []
    dv_rms_sigmas = []
    calib_vals  = []
    calib_sigmas = []

    for name, info in fit.items():
        if name not in inner_mods:
            continue
        if abs(info["ratio"] - 1.0) <= RATIO_TOL and name in hv:
            hv_vals.append(hv[name])
            sigma_vals.append(info["sigma"])
            sigma_map[name] = info["sigma"]
            hv_map[name]    = hv[name]
            if name in dv_rms:
                dv_rms_vals.append(dv_rms[name])
                dv_rms_sigmas.append(info["sigma"])
            if name in calib and calib[name] < 0.15:
                calib_vals.append(calib[name])
                calib_sigmas.append(info["sigma"])

    # Remove dV RMS outliers: hard cutoff at DV_RMS_MAX
    if dv_rms_vals:
        dv_arr_all = np.array(dv_rms_vals)
        mask = dv_arr_all <= DV_RMS_MAX
        n_removed = np.sum(~mask)
        if n_removed:
            print(f"Removed {n_removed} dV RMS outlier(s) > {DV_RMS_MAX} V")
        dv_rms_vals   = dv_arr_all[mask].tolist()
        dv_rms_sigmas = np.array(dv_rms_sigmas)[mask].tolist()

    print(f"Excluding outer {EXCLUDE_LAYERS} layers, |ratio-1| <= {RATIO_TOL}: {len(hv_vals)} modules")

    hv_arr    = np.array(hv_vals)
    sigma_arr = np.array(sigma_vals)
    coeffs    = np.polyfit(hv_arr, sigma_arr, 1)
    slope, intercept = coeffs
    x_fit = np.linspace(hv_arr.min(), hv_arr.max(), 200)
    y_fit = np.polyval(coeffs, x_fit)
    print(f"Linear fit: sigma = {slope:.4f} * HV + {intercept:.2f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(hv_vals, sigma_vals, s=20, alpha=0.7, color="steelblue", label="data")
    ax.plot(x_fit, y_fit, color="tomato", linewidth=1.8,
            label=f"fit: {slope:.4f}·HV + {intercept:.2f}")

    ax.set_xlabel("HV V0Set (V)", fontsize=13)
    ax.set_ylabel("Sigma (MeV)", fontsize=13)
    ax.set_title(f"HV vs Fit Sigma  (inner, excl. outer {EXCLUDE_LAYERS} layers, |ratio−1| ≤ {RATIO_TOL})", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    out = "hv_vs_sigma.png"
    plt.savefig(out, dpi=150)
    print(f"Saved: {out}")

    # --- dV RMS vs Sigma scatter ---
    if dv_rms_vals:
        dv_arr2   = np.array(dv_rms_vals)
        sig_arr2  = np.array(dv_rms_sigmas)
        coeffs2   = np.polyfit(dv_arr2, sig_arr2, 1)
        x_fit2    = np.linspace(dv_arr2.min(), dv_arr2.max(), 200)
        y_fit2    = np.polyval(coeffs2, x_fit2)
        slope2, intercept2 = coeffs2
        print(f"dV RMS linear fit: sigma = {slope2:.4f} * dV_RMS + {intercept2:.2f}")

        fig3, ax3 = plt.subplots(figsize=(8, 6))
        ax3.scatter(dv_rms_vals, dv_rms_sigmas, s=20, alpha=0.7, color="steelblue", label="data")
        ax3.plot(x_fit2, y_fit2, color="tomato", linewidth=1.8,
                 label=f"fit: {slope2:.4f}·dV_RMS + {intercept2:.2f}")
        ax3.set_xlabel("dV RMS (V)", fontsize=13)
        ax3.set_ylabel("Sigma (ADC counts)", fontsize=13)
        ax3.set_title(
            f"dV RMS vs Fit Sigma  (inner, excl. outer {EXCLUDE_LAYERS} layers, |ratio−1|≤{RATIO_TOL})",
            fontsize=13)
        ax3.legend(fontsize=11)
        ax3.grid(True, linestyle="--", alpha=0.4)
        fig3.tight_layout()
        out3 = "dvrms_vs_sigma.png"
        fig3.savefig(out3, dpi=150)
        print(f"Saved: {out3}")

    # --- Sigma vs Calib Factor scatter ---
    if calib_vals:
        cf_arr  = np.array(calib_vals)
        sig_cf  = np.array(calib_sigmas)
        coeffs_cf = np.polyfit(cf_arr, sig_cf, 1)
        x_cf = np.linspace(cf_arr.min(), cf_arr.max(), 200)
        y_cf = np.polyval(coeffs_cf, x_cf)
        slope_cf, intercept_cf = coeffs_cf
        print(f"Calib factor linear fit: sigma = {slope_cf:.4f} * factor + {intercept_cf:.2f}")

        fig4, ax4 = plt.subplots(figsize=(8, 6))
        ax4.scatter(calib_vals, calib_sigmas, s=20, alpha=0.7, color="steelblue", label="data")
        ax4.plot(x_cf, y_cf, color="tomato", linewidth=1.8,
                 label=f"fit: {slope_cf:.4f}·factor + {intercept_cf:.2f}")
        ax4.set_xlabel("Calib Factor", fontsize=13)
        ax4.set_ylabel("Sigma (ADC counts)", fontsize=13)
        ax4.set_title(
            f"Calib Factor vs Fit Sigma  (inner, excl. outer {EXCLUDE_LAYERS} layers, |ratio−1|≤{RATIO_TOL})",
            fontsize=13)
        ax4.legend(fontsize=11)
        ax4.grid(True, linestyle="--", alpha=0.4)
        fig4.tight_layout()
        out4 = "calib_vs_sigma.png"
        fig4.savefig(out4, dpi=150)
        print(f"Saved: {out4}")

    # --- 2D maps ---
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    plot_2d_map(ax1, mod_pos, sigma_map,
                f"Sigma map (ADC counts)\n(excl. outer {EXCLUDE_LAYERS} layers, |ratio−1|≤{RATIO_TOL})",
                cmap="RdYlBu_r")
    plot_2d_map(ax2, mod_pos, hv_map,
                f"HV V0Set map (V)\n(excl. outer {EXCLUDE_LAYERS} layers, |ratio−1|≤{RATIO_TOL})",
                cmap="RdYlBu_r")
    fig2.tight_layout()
    out2 = "hv_sigma_2dmap.png"
    fig2.savefig(out2, dpi=150)
    print(f"Saved: {out2}")

    plt.show()


if __name__ == "__main__":
    main()
