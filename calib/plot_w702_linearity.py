#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Dict, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import uproot

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


ROOT_FILE = Path("build/nonlinearity_results.root")
MAP_FILE = Path("database/hycal_map.json")
OUT_FILE = Path("W702_linearity.png")
MODULE = "W702"

Z_HYCAL = 6270.0
M_PROTON = 938.2720813
M_ELECTRON = 0.51099895
RESOLUTION = 0.035

PEAKS = [
    ("ep_3p5", "ep 3.5 GeV", 3485.41, "ep", "energy_3p5GeV", "3p5", "tomato"),
    ("ee_3p5", "ee 3.5 GeV", 3485.41, "ee", "energy_3p5GeV", "3p5", "cornflowerblue"),
    ("ep_2p2", "ep 2.2 GeV", 2239.51, "ep", "energy_2p2GeV", "2p2", "tomato"),
    ("ee_2p2", "ee 2.2 GeV", 2239.51, "ee", "energy_2p2GeV", "2p2", "cornflowerblue"),
    ("ep_0p7", "ep 0.7 GeV", 728.9, "ep", "energy_0p7GeV", "0p7", "tomato"),
    ("ee_0p7", "ee 0.7 GeV", 728.9, "ee", "energy_0p7GeV", "0p7", "cornflowerblue"),
]


def energy_loss(theta_deg: float) -> float:
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sec = 1.0 / cos_t if cos_t > 0.01 else 100.0
    eloss = 0.500 * 1.6 * sec
    eloss += 0.120 * 1.6 * sec
    eloss += 0.100 * 2.0 * sec
    eloss += 0.480 * 1.8 * sec
    return eloss


def expected_energy(theta_deg: float, ebeam: float, kind: str) -> float:
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    if kind == "ep":
        e = ebeam * M_PROTON / (M_PROTON + ebeam * (1.0 - cos_t))
    elif kind == "ee":
        gamma = ebeam / M_ELECTRON
        num = (gamma + 1.0) + (gamma - 1.0) * cos_t * cos_t
        den = (gamma + 1.0) - (gamma - 1.0) * cos_t * cos_t
        e = M_ELECTRON * num / den if den > 0.0 else 0.0
    else:
        e = 0.0
    return max(0.0, e - energy_loss(theta_deg))


def estimate_sigma(e: float) -> float:
    return e * RESOLUTION / math.sqrt(max(e, 1.0) / 1000.0)


def gauss(x: np.ndarray, amp: float, mu: float, sig: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def gauss_fit_window(
    centers: np.ndarray,
    counts: np.ndarray,
    mu0: float,
    half_width: float,
) -> Tuple[float, float, float] | None:
    lo, hi = mu0 - half_width, mu0 + half_width
    mask = (centers >= lo) & (centers <= hi) & (counts > 0)
    if mask.sum() < 4 or not HAS_SCIPY:
        return None
    xf, yf = centers[mask], counts[mask].astype(float)
    amp0 = float(yf.max())
    bw = float(centers[1] - centers[0]) if len(centers) > 1 else 1.0
    sig0 = max(half_width / 3.0, bw)
    try:
        popt, _ = curve_fit(
            gauss,
            xf,
            yf,
            p0=[amp0, mu0, sig0],
            bounds=([0.0, lo, 1e-6], [np.inf, hi, half_width * 2.0]),
            maxfev=3000,
        )
        return float(popt[1]), abs(float(popt[2])), float(popt[0])
    except Exception:
        return None


def fit_peak_full(counts: np.ndarray, edges: np.ndarray, center: float) -> Tuple[float, float, float]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    sigma = estimate_sigma(center)
    mask = (centers >= center - 6.0 * sigma) & (centers <= center + 6.0 * sigma)
    w = counts[mask]
    x = centers[mask]
    if w.sum() <= 0.0:
        return 0.0, 0.0, 0.0

    mean = float((x * w).sum() / w.sum())
    amp = float(w.max())
    best = None
    r = gauss_fit_window(centers, counts, mean, 2.0 * estimate_sigma(mean))
    if r is not None:
        mean = r[0]
        best = r
    r = gauss_fit_window(centers, counts, mean, estimate_sigma(mean))
    if r is not None:
        return r
    if best is not None:
        return best
    return mean, estimate_sigma(mean), amp


def nl_ratio_1st(e_rec: np.ndarray, nl1: float, e_base: float) -> np.ndarray:
    return 1.0 + nl1 * (e_rec - e_base) / 1000.0


def nl_ratio_2nd(e_rec: np.ndarray, nl1: float, nl2: float, e_base: float) -> np.ndarray:
    t = (e_rec - e_base) / 1000.0
    return 1.0 + nl1 * t + nl2 * t * t


def fit_linearity(e_exp: np.ndarray, e_rec: np.ndarray, e_base: float):
    ratio = e_rec / e_exp
    nl1, nl1_err, chi1, ndf1 = 0.0, 0.0, 0.0, 0
    nl21, nl21_err, nl22, nl22_err, chi2, ndf2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0
    if HAS_SCIPY and len(e_rec) >= 2:
        popt, pcov = curve_fit(lambda x, a: nl_ratio_1st(x, a, e_base), e_rec, ratio, p0=[0.01])
        nl1 = float(popt[0])
        nl1_err = float(np.sqrt(max(pcov[0, 0], 0.0)))
        chi1 = float(np.sum((ratio - nl_ratio_1st(e_rec, nl1, e_base)) ** 2))
        ndf1 = len(e_rec) - 1
    if HAS_SCIPY and len(e_rec) >= 3:
        popt, pcov = curve_fit(
            lambda x, a, b: nl_ratio_2nd(x, a, b, e_base),
            e_rec,
            ratio,
            p0=[nl1, 0.0],
            maxfev=5000,
        )
        nl21, nl22 = float(popt[0]), float(popt[1])
        err = np.sqrt(np.maximum(np.diag(pcov), 0.0))
        nl21_err, nl22_err = float(err[0]), float(err[1])
        chi2 = float(np.sum((ratio - nl_ratio_2nd(e_rec, nl21, nl22, e_base)) ** 2))
        ndf2 = len(e_rec) - 2
    return nl1, nl1_err, chi1, ndf1, nl21, nl21_err, nl22, nl22_err, chi2, ndf2


def w702_theta() -> float:
    with MAP_FILE.open() as f:
        modules = json.load(f)
    mod = next(m for m in modules if m.get("n") == MODULE)
    geo = mod["geo"]
    return math.degrees(math.atan2(math.hypot(float(geo["x"]), float(geo["y"])), Z_HYCAL))


def load_peaks(root_file: Path, theta: float) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, str]]:
    peaks: Dict[str, float] = {}
    expected: Dict[str, float] = {}
    labels: Dict[str, str] = {}
    with uproot.open(root_file) as f:
        for key, label, ebeam, kind, root_dir, suffix, _color in PEAKS:
            counts, edges = f[f"{root_dir}/h_energy_{MODULE}_{suffix}"].to_numpy()
            e_exp = expected_energy(theta, ebeam, kind)
            fit_center = 545.0 if key == "ee_0p7" else e_exp
            mu, _sig, _amp = fit_peak_full(counts, edges, fit_center)
            if mu > 0.0:
                peaks[key] = mu
                expected[key] = e_exp
                labels[key] = label
                note = "refit at 545 MeV" if key == "ee_0p7" else "auto"
                print(f"{key:7s} exp={e_exp:8.2f}  peak={mu:8.2f}  ({note})")
            else:
                print(f"{key:7s} exp={e_exp:8.2f}  peak=N/A")
    return peaks, expected, labels


def main() -> None:
    theta = w702_theta()
    peaks, expected, labels = load_peaks(ROOT_FILE, theta)
    e_exp = np.array([expected[k] for k in peaks], dtype=float)
    e_rec = np.array([peaks[k] for k in peaks], dtype=float)
    e_base = expected["ep_3p5"]

    nl1, nl1_err, chi1, ndf1, nl21, nl21_err, nl22, nl22_err, chi2, ndf2 = fit_linearity(e_exp, e_rec, e_base)

    fig, ax = plt.subplots(figsize=(10.8, 6.4), facecolor="white")
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color("#333333")
    ax.tick_params(colors="#333333", labelsize=12)
    ax.set_title(f"Non-linearity of {MODULE}", color="#111111", fontsize=18, pad=10)
    ax.set_xlabel(r"$E_{\mathrm{reconstruct}}$ (MeV)", color="#111111", fontsize=16)
    ax.set_ylabel(r"$E_{\mathrm{reconstruct}} / E_{\mathrm{expect}}$", color="#111111", fontsize=16)
    ax.set_ylim(0.92, 1.025)
    ax.grid(color="#999999", alpha=0.35, linewidth=0.8)

    margin = (e_rec.max() - e_rec.min()) * 0.15
    xplot = np.linspace(max(0.0, e_rec.min() - margin), e_rec.max() + margin, 400)
    ax.axhline(1.0, color="tomato", linestyle="--", linewidth=1.5, alpha=0.75, label="y = 1 (ideal)")
    ax.plot(xplot, nl_ratio_1st(xplot, nl1, e_base), color="cornflowerblue", linewidth=2.0, label="1st order fit")
    ax.plot(
        xplot,
        nl_ratio_2nd(xplot, nl21, nl22, e_base),
        color="violet",
        linewidth=2.0,
        linestyle="--",
        label="2nd order fit",
    )

    for key, label, _ebeam, kind, _root_dir, _suffix, color in PEAKS:
        if key not in peaks:
            continue
        ratio = peaks[key] / expected[key]
        ax.scatter(peaks[key], ratio, color=color, marker="o", s=95, zorder=5, edgecolors="white", linewidths=0.5)
        text = label
        ax.annotate(
            text,
            xy=(peaks[key], ratio),
            xytext=(0, -20),
            textcoords="offset points",
            color=color,
            fontsize=16,
            ha="center",
            va="top",
        )

    for key in peaks:
        denom1 = 1.0 + nl1 * (peaks[key] - e_base) / 1000.0
        ec1 = peaks[key] / denom1 if denom1 != 0.0 else peaks[key]
        ax.scatter(ec1, ec1 / expected[key], color="limegreen", marker="s", s=81, zorder=6, edgecolors="white", linewidths=0.5)

        t = (peaks[key] - e_base) / 1000.0
        denom2 = 1.0 + nl21 * t + nl22 * t * t
        ec2 = peaks[key] / denom2 if denom2 != 0.0 else peaks[key]
        ax.scatter(ec2, ec2 / expected[key], color="orange", marker="s", s=81, zorder=6, edgecolors="white", linewidths=0.5)

    ax.scatter([], [], color="tomato", marker="o", s=95, label="e-p (measured)")
    ax.scatter([], [], color="cornflowerblue", marker="o", s=95, label="e-e (measured)")
    ax.scatter([], [], color="limegreen", marker="s", s=81, label="corrected (1st)")
    ax.scatter([], [], color="orange", marker="s", s=81, label="corrected (2nd)")
    ax.legend(facecolor="white", edgecolor="#555555", labelcolor="#111111", fontsize=16, loc="lower right")

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=180, facecolor="white", transparent=False)
    print(f"saved {OUT_FILE}")


if __name__ == "__main__":
    main()
