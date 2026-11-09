"""Figures for the paper.

Usage:
  python src/09_figures.py <processed_dir> <results_dir> <figures_dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "axes.titlesize": 9, "axes.titleweight": "bold", "figure.dpi": 150,
    "savefig.dpi": 300, "savefig.bbox": "tight", "lines.linewidth": 1.5,
})
FAM_LABEL = {"salience": "Salience (“heatwave”)", "coping": "Coping index", "concern": "Concern index"}


def fig_heat(proc, out):
    d = pd.read_csv(proc / "heat_daily.csv", parse_dates=["date"])
    d = d[(d.series == "UK") & (d.percentile == "p95")]
    years = [y for y in [2022, 2023, 2024, 2025, 2026] if ((d.date.dt.year == y) & (d.date.dt.month == 7)).any()]
    fig, axes = plt.subplots(len(years), 1, figsize=(7, 7.5), sharey=True)
    for ax, y in zip(axes, years):
        s = d[(d.date >= f"{y}-04-01") & (d.date <= f"{y}-10-31")]
        ax.plot(s.date, s.tmax_degC, color=BLUE, lw=1.0)
        thr = s.threshold_degC.iloc[0]
        ax.axhline(thr, color=INK2, lw=0.8, ls="--")
        hw = s[s.heatwave_day]
        for day in hw.date:
            ax.axvspan(day - pd.Timedelta(hours=12), day + pd.Timedelta(hours=12), color=ORANGE, alpha=0.35, lw=0)
        ax.set_title(str(y), loc="left")
        ax.set_ylabel("°C")
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
    fig.suptitle("UK daily maximum temperature (population-weighted), April–October\n"
                 f"dashed: heatwave threshold {thr:.1f}°C (95th percentile, June–August 1991–2020); "
                 "shaded: heatwave days", x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig1_heat_timeline.png")
    plt.close(fig)


def fig_event(res, out, rule):
    c = pd.read_csv(res / "analysis1_event_coefficients.csv")
    bands = pd.read_csv(res / "placebo_beta_bands.csv")
    fams = [f for f in ["salience", "coping", "concern"] if f"trends_{rule}_{f}" in set(c.outcome)]
    if not fams:
        return
    fig, axes = plt.subplots(2, len(fams), figsize=(2.6 * len(fams) + 0.6, 5), sharex="row", squeeze=False)
    for j, f in enumerate(fams):
        n = f"trends_{rule}_{f}"
        b = bands[bands.outcome == n].set_index("k")
        for i, (model, lab) in enumerate([("1a", "weeks since onset"), ("1b", "weeks since end")]):
            ax = axes[i, j]
            s = c[(c.outcome == n) & (c.model == model)].sort_values("k")
            bb = b.reindex(s.k)
            ax.fill_between(s.k, bb.q025.fillna(0), bb.q975.fillna(0), color=GRID, alpha=0.9, lw=0,
                            label="95% placebo range")
            ax.axhline(0, color=INK2, lw=0.8)
            ax.plot(s.k, s.beta, color=BLUE, marker="o", ms=4, label="estimate")
            ax.set_xlabel(lab)
            if j == 0:
                ax.set_ylabel("log points vs k = −1")
            if i == 0:
                ax.set_title(FAM_LABEL[f])
    axes[0, 0].legend(frameon=False, fontsize=7, loc="upper left")
    tag = "pre-registered sparsity rule" if rule == "plan" else "all pre-registered terms, deviation 2"
    fig.suptitle(f"Event-study estimates, UK Google Trends ({tag})", x=0.01, ha="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"fig2_event_study_{rule}.png")
    plt.close(fig)


def fig_dose(res, out, rule):
    fams = [f for f in ["coping", "concern"] if (res / f"analysis2_lags_{rule}_{f}_main.csv").exists()]
    if not fams:
        return
    fig, axes = plt.subplots(1, len(fams), figsize=(3.2 * len(fams), 2.8), squeeze=False)
    for j, f in enumerate(fams):
        ax = axes[0, j]
        l = pd.read_csv(res / f"analysis2_lags_{rule}_{f}_main.csv")
        l = l[l.term.str.fullmatch(r"d\d+")].assign(k=lambda d: d.term.str[1:].astype(int))
        pb = pd.read_csv(res / f"analysis2_permbands_{rule}_{f}_main.csv")
        ax.fill_between(pb.k, pb.perm_min, pb.perm_max, color=GRID, alpha=0.9, lw=0, label="range over 23 permutations")
        ax.axhline(0, color=INK2, lw=0.8)
        ax.plot(l.k, l.beta, color=BLUE, marker="o", ms=4, label="estimate")
        ax.set_title(FAM_LABEL[f])
        ax.set_xlabel("weeks after heat dose")
        if j == 0:
            ax.set_ylabel("log points per 10 degree-days")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Four-nation dose-response (within-week differences between nations)", x=0.01, ha="left",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"fig3_dose_response_{rule}.png")
    plt.close(fig)


def fig_size(res, out):
    s = pd.read_csv(res / "size_check.csv")
    fig, ax = plt.subplots(figsize=(5, 2.2))
    ax.barh(s.test, s.rejection_rate_5pct, color=BLUE, height=0.5)
    ax.axvline(0.05, color=ORANGE, lw=1.5)
    ax.text(0.052, -0.6, "nominal 5%", color=INK2, fontsize=7)
    for y, v in enumerate(s.rejection_rate_5pct):
        ax.text(v + 0.004, y, f"{v:.0%}", va="center", fontsize=7, color=INK)
    ax.set_xlabel("rejection rate on white noise")
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_title("Model-based tests are oversized with few events", loc="left")
    fig.tight_layout()
    fig.savefig(out / "fig4_size_check.png")
    plt.close(fig)


def main(proc, res, out):
    proc, res, out = Path(proc), Path(res), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fig_heat(proc, out)
    for rule in ["plan", "dev2"]:
        fig_event(res, out, rule)
        fig_dose(res, out, rule)
    if (res / "size_check.csv").exists():
        fig_size(res, out)
    print("figures written to", out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
