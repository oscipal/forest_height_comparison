"""Pooled, multi-scene figures.

These take a long-form table of paired cells drawn from every scene and encode
the grouping variable (site, or scene) as *identity* colour: fixed categorical
slots, assigned in order and never cycled. Past the available slots the palette
stops being separable, so grouping is refused rather than silently degraded.

Everything else -- the style, the title block, the metric definitions -- is
shared with the per-scene figures in ``bgt.viz``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import config
from bgt.metrics import compute_metrics
from bgt.viz import (
    AXIS,
    INK,
    INK_MUTED,
    INK_SECONDARY,
    SERIES,
    SURFACE,
    _annotate,
    _save,
    _titleblock,
)


def group_colors(groups: list[str]) -> dict[str, str]:
    """Map group labels to categorical slots, in fixed order."""
    if len(groups) > len(SERIES):
        raise ValueError(
            f"{len(groups)} groups exceeds the {len(SERIES)} categorical slots; "
            "fold the smallest into 'Other' or facet instead of cycling hues."
        )
    return {g: SERIES[i] for i, g in enumerate(groups)}


def _groups_of(df: pd.DataFrame, group_col: str) -> tuple[list[str], dict[str, str]]:
    groups = sorted(df[group_col].dropna().unique().tolist())
    return groups, group_colors(groups)


def plot_scatter_grouped(
    df: pd.DataFrame,
    ref_col: str,
    out_path: Path,
    ref_label: str,
    group_col: str = "site",
    prod_col: str = "fh",
    prod_label: str = "BIOMASS L2A forest height (m)",
    subtitle: str = "",
) -> Path:
    """Pooled scatter coloured by group; the fits are on the pooled sample."""
    work = df[[ref_col, prod_col, group_col]].dropna()
    x, y = work[ref_col].to_numpy(), work[prod_col].to_numpy()
    m = compute_metrics(x, y)
    groups, colors = _groups_of(work, group_col)

    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo)
    lims = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(6.8, 6.6), constrained_layout=True)
    for g in groups:
        sub = work[work[group_col] == g]
        ax.scatter(sub[ref_col], sub[prod_col], s=9, color=colors[g], alpha=0.35,
                   linewidths=0, label=f"{g}  (n = {len(sub)})", zorder=3)

    grid = np.array(lims)
    ax.plot(grid, grid, color=INK_MUTED, lw=1.4, ls="--", label="1:1", zorder=5)
    ax.plot(grid, m["ols_slope"] * grid + m["ols_intercept"], color=INK, lw=1.8,
            label=f"OLS pooled  y = {m['ols_slope']:.2f}x {m['ols_intercept']:+.1f}",
            zorder=6)
    ax.plot(grid, m["rma_slope"] * grid + m["rma_intercept"], color=INK, lw=1.8,
            ls=":",
            label=f"RMA pooled  y = {m['rma_slope']:.2f}x {m['rma_intercept']:+.1f}",
            zorder=6)

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.set_xlabel(ref_label)
    ax.set_ylabel(prod_label)
    legend = ax.legend(loc="lower right", markerscale=2)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    _annotate(
        ax,
        "\n".join([
            f"pooled over {len(groups)} sites",
            f"n = {m['n']:.0f} cells",
            f"bias = {m['bias']:+.2f} m",
            f"MAE = {m['mae']:.2f} m",
            f"RMSE = {m['rmse']:.2f} m ({m['rmse_pct']:.0f} %)",
            f"r = {m['pearson_r']:.3f}   CCC = {m['ccc']:.3f}",
            f"R2 vs 1:1 = {m['r2_one_to_one']:.3f}",
        ]),
    )
    _titleblock(fig, "BIOMASS vs ALS canopy height, all scenes pooled", subtitle)
    return _save(fig, out_path)


def plot_residuals_grouped(
    df: pd.DataFrame,
    ref_col: str,
    out_path: Path,
    ref_label: str,
    group_col: str = "site",
    prod_col: str = "fh",
    bins: list[float] | None = None,
    subtitle: str = "",
) -> Path:
    """Residual against reference height, with one binned mean line per group."""
    bins = bins or config.HEIGHT_BINS
    work = df[[ref_col, prod_col, group_col]].dropna().copy()
    work["resid"] = work[prod_col] - work[ref_col]
    work["bin"] = pd.cut(work[ref_col], bins=bins, right=False)
    groups, colors = _groups_of(work, group_col)

    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    ax.axhline(0.0, color=INK_MUTED, lw=1.2, ls="--", zorder=2)
    for g in groups:
        sub = work[work[group_col] == g]
        ax.scatter(sub[ref_col], sub["resid"], s=7, color=colors[g], alpha=0.18,
                   linewidths=0, zorder=3)
    for g in groups:
        sub = work[work[group_col] == g]
        agg = sub.groupby("bin", observed=True)["resid"].agg(["mean", "std", "count"])
        agg = agg[agg["count"] >= 3]
        if agg.empty:
            continue
        centres = np.array([(iv.left + iv.right) / 2.0 for iv in agg.index])
        ax.errorbar(centres, agg["mean"], yerr=agg["std"], color=colors[g], lw=2.0,
                    marker="o", ms=6, capsize=4, elinewidth=1.4, zorder=5,
                    markeredgecolor=SURFACE, markeredgewidth=1.2, label=g)
    ax.set_xlabel(ref_label)
    ax.set_ylabel("residual, BIOMASS - ALS (m)")
    ax.legend(loc="upper right", title="binned mean +/- 1 SD", alignment="left")
    _annotate(ax, "positive = BIOMASS taller than ALS", loc="lower left")
    _titleblock(fig, "Residual structure by site, all scenes pooled", subtitle)
    return _save(fig, out_path)


def plot_distributions_grouped(
    df: pd.DataFrame,
    ref_col: str,
    out_path: Path,
    ref_label: str,
    group_col: str = "site",
    prod_col: str = "fh",
    subtitle: str = "",
) -> Path:
    """Pooled histograms, plus per-site cumulative distributions of both estimates."""
    work = df[[ref_col, prod_col, group_col]].dropna()
    x, y = work[ref_col].to_numpy(), work[prod_col].to_numpy()
    edges = np.histogram_bin_edges(np.concatenate([x, y]), bins=40)
    groups, colors = _groups_of(work, group_col)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.8), constrained_layout=True)

    ax1.hist(x, bins=edges, color=SERIES[0], alpha=0.55, label="ALS reference")
    ax1.hist(y, bins=edges, color=SERIES[1], alpha=0.55, label="BIOMASS L2A FH")
    ax1.set_xlabel("canopy height (m)")
    ax1.set_ylabel("number of BIOMASS cells")
    ax1.set_title("Pooled height distributions")
    ax1.legend(loc="upper right")

    for g in groups:
        sub = work[work[group_col] == g]
        for col, style in ((ref_col, "-"), (prod_col, "--")):
            s = np.sort(sub[col].to_numpy())
            ax2.plot(s, np.arange(1, s.size + 1) / s.size, color=colors[g], ls=style,
                     lw=1.8, label=g if col == ref_col else None)
    ax2.set_xlabel("canopy height (m)")
    ax2.set_ylabel("cumulative fraction of cells")
    ax2.set_title("Cumulative distributions per site")
    ax2.legend(loc="lower right", title="solid ALS / dashed BIOMASS", alignment="left")
    _titleblock(fig, "Height distributions, all scenes pooled", subtitle)
    return _save(fig, out_path)


def plot_bland_altman_grouped(
    df: pd.DataFrame,
    ref_col: str,
    out_path: Path,
    group_col: str = "site",
    prod_col: str = "fh",
    subtitle: str = "",
) -> Path:
    """Bland-Altman over the pooled sample, points coloured by group."""
    work = df[[ref_col, prod_col, group_col]].dropna().copy()
    work["mean_h"] = (work[ref_col] + work[prod_col]) / 2.0
    work["diff"] = work[prod_col] - work[ref_col]
    bias, sd = work["diff"].mean(), work["diff"].std(ddof=1)
    groups, colors = _groups_of(work, group_col)

    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    for g in groups:
        sub = work[work[group_col] == g]
        ax.scatter(sub["mean_h"], sub["diff"], s=8, color=colors[g], alpha=0.25,
                   linewidths=0, label=g)
    for value, style, name in (
        (bias, "-", f"pooled bias {bias:+.2f} m"),
        (bias + 1.96 * sd, "--", f"+1.96 SD  {bias + 1.96 * sd:+.2f} m"),
        (bias - 1.96 * sd, "--", f"-1.96 SD  {bias - 1.96 * sd:+.2f} m"),
    ):
        ax.axhline(value, color=INK, lw=1.6, ls=style, zorder=6)
        ax.annotate(name, xy=(0.995, value), xycoords=("axes fraction", "data"),
                    ha="right", va="bottom", fontsize=8, color=INK_SECONDARY)
    ax.set_xlabel("mean of ALS and BIOMASS height (m)")
    ax.set_ylabel("difference, BIOMASS - ALS (m)")
    legend = ax.legend(loc="lower left", markerscale=2)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    _titleblock(fig, "Bland-Altman agreement, all scenes pooled", subtitle)
    return _save(fig, out_path)


def plot_scene_metrics(
    table: pd.DataFrame,
    out_path: Path,
    label_col: str = "scene_label",
    group_col: str = "site",
    subtitle: str = "",
) -> Path:
    """Per-scene bias, RMSE and correlation as bars, coloured by site."""
    tbl = table.sort_values([group_col, "biomass_date"]).reset_index(drop=True)
    groups, colors = _groups_of(tbl, group_col)
    bar_colors = [colors[g] for g in tbl[group_col]]
    pos = np.arange(len(tbl))

    fig, axes = plt.subplots(3, 1, figsize=(9.6, 9.2), constrained_layout=True,
                             sharex=True)
    panels = [
        ("bias", "bias (m)", "zero is best"),
        ("rmse", "RMSE (m)", "lower is better"),
        ("pearson_r", "Pearson r", "higher is better"),
    ]
    for ax, (col, label, hint) in zip(axes, panels):
        values = tbl[col].to_numpy()
        ax.bar(pos, values, width=0.66, color=bar_colors, linewidth=0)
        for p, v in zip(pos, values):
            ax.annotate(f"{v:.2f}", xy=(p, v), xytext=(0, 3 if v >= 0 else -11),
                        textcoords="offset points", ha="center", fontsize=8,
                        color=INK_SECONDARY)
        ax.axhline(0, color=AXIS, lw=1.0)
        ax.set_ylabel(label)
        ax.set_title(f"{label} - {hint}")
        ax.margins(y=0.20)

    # A figure-level legend: every panel's bars can reach the top of its axes,
    # so there is no reliably free corner inside one.
    handles = [plt.Line2D([], [], marker="s", ls="", ms=9, color=colors[g], label=g)
               for g in groups]
    fig.legend(handles=handles, loc="outside upper right", ncol=len(groups),
               frameon=False)
    axes[-1].set_xticks(pos)
    axes[-1].set_xticklabels(tbl[label_col], rotation=30, ha="right", fontsize=8)
    axes[-1].set_xlabel("BIOMASS scene")
    _titleblock(fig, "Agreement per scene", subtitle)
    return _save(fig, out_path)


def plot_exclusions_by_scene(
    table: pd.DataFrame,
    labels: dict[str, str],
    colors: dict[str, str],
    out_path: Path,
    label_col: str = "scene_label",
    subtitle: str = "",
    grid_name: str = "BIOMASS",
) -> Path:
    """Where every scene's cells went, as a share of its analysis window.

    Stacked to 100 % rather than to absolute counts: the windows differ in size
    between sites, and the question here is what proportion of each scene is
    lost to what, not which site has the bigger rectangle. Absolute counts are
    in ``exclusion_counts_all.csv``.
    """
    keys = [k for k in labels if k in table.columns]
    tbl = table.sort_values([label_col]).reset_index(drop=True)
    totals = tbl[keys].sum(axis=1).replace(0, np.nan)
    pos = np.arange(len(tbl))

    fig, ax = plt.subplots(figsize=(10.4, 0.62 * len(tbl) + 3.0),
                           constrained_layout=True)
    left = np.zeros(len(tbl))
    for key in keys:
        share = 100.0 * tbl[key] / totals
        ax.barh(pos, share, left=left, height=0.68, color=colors.get(key, INK_MUTED),
                linewidth=0.8, edgecolor=SURFACE, label=labels[key])
        for p, (s, l) in enumerate(zip(share, left)):
            if s >= 6:  # only label a segment wide enough to hold the text
                ax.annotate(f"{s:.0f}%", xy=(l + s / 2, p), ha="center",
                            va="center", fontsize=7.5, color=INK_SECONDARY)
        left = left + share.fillna(0).to_numpy()

    ax.set_yticks(pos)
    ax.set_yticklabels(tbl[label_col], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel(f"share of the scene's analysis window (%)")
    ax.grid(axis="y", visible=False)
    for p, total in enumerate(tbl[keys].sum(axis=1)):
        ax.annotate(f"{int(total):,} cells", xy=(100.6, p), va="center",
                    fontsize=7.5, color=INK_MUTED)
    ax.margins(x=0.12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16 / (len(tbl) ** 0.5 + 1) - 0.06),
              ncol=3, frameon=False)
    _titleblock(fig, f"Where the {grid_name} cells go, scene by scene", subtitle)
    return _save(fig, out_path)


def plot_site_metrics(
    table: pd.DataFrame,
    out_path: Path,
    group_col: str = "site",
    subtitle: str = "",
) -> Path:
    """Per-site pooled bias, RMSE, correlation and sample size."""
    tbl = table.sort_values(group_col).reset_index(drop=True)
    groups, colors = _groups_of(tbl, group_col)
    bar_colors = [colors[g] for g in tbl[group_col]]
    pos = np.arange(len(tbl))

    fig, axes = plt.subplots(1, 4, figsize=(13.6, 4.2), constrained_layout=True)
    panels = [
        ("bias", "bias (m)", "zero is best"),
        ("rmse", "RMSE (m)", "lower is better"),
        ("pearson_r", "Pearson r", "higher is better"),
        ("n", "paired cells", "sample size"),
    ]
    for ax, (col, label, hint) in zip(axes, panels):
        values = tbl[col].to_numpy()
        ax.bar(pos, values, width=0.62, color=bar_colors, linewidth=0)
        fmt = "{:.0f}" if col == "n" else "{:.2f}"
        for p, v in zip(pos, values):
            ax.annotate(fmt.format(v), xy=(p, v), xytext=(0, 3 if v >= 0 else -11),
                        textcoords="offset points", ha="center", fontsize=8,
                        color=INK_SECONDARY)
        ax.axhline(0, color=AXIS, lw=1.0)
        ax.set_xticks(pos)
        ax.set_xticklabels(tbl[group_col], rotation=20, ha="right", fontsize=8)
        ax.set_ylabel(label)
        ax.set_title(f"{label} - {hint}")
        ax.margins(y=0.20)
    _titleblock(fig, "Agreement per site, scenes pooled within each site", subtitle)
    return _save(fig, out_path)
