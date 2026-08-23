"""Figures for the ALS / BIOMASS forest height comparison.

Colour is assigned by the job it does, never decoratively:

* **magnitude** (height maps, density) -- one hue, light to dark;
* **polarity** (residuals) -- a two-hue diverging ramp with a neutral grey
  midpoint, always symmetric about zero so the midpoint means "no difference";
* **identity** (two distributions, several candidate statistics) -- fixed
  categorical slots in a fixed order, always paired with a legend.

Marks are thin, gridlines recessive, and labels carry text ink rather than the
series colour.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

import config
from bgt.metrics import compute_metrics

# --------------------------------------------------------------------------- #
# Design tokens
# --------------------------------------------------------------------------- #

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
BLUE, ORANGE, AQUA = SERIES[0], SERIES[1], SERIES[2]
RED = "#e34948"

#: Sequential ramp for magnitude (canopy height, point density).
CMAP_HEIGHT = LinearSegmentedColormap.from_list(
    "gca_blue",
    ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)
#: Diverging ramp for residuals. Grey midpoint so "no difference" reads as nothing.
CMAP_RESIDUAL = LinearSegmentedColormap.from_list(
    "gca_div",
    ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f0efec",
     "#f7c0bf", "#ec8281", "#d03b3b", "#8b1f1f"],
)

_RC = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.titleweight": "600",
    "axes.titlelocation": "left",
    "axes.titlepad": 8,
    "axes.labelsize": 9,
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.7,
    "xtick.color": INK_MUTED,
    "ytick.color": INK_MUTED,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "text.color": INK,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "lines.linewidth": 2.0,
    "figure.dpi": config.FIG_DPI,
    "savefig.dpi": config.FIG_DPI,
    "savefig.bbox": "tight",
}


def use_style() -> None:
    """Apply the figure style. Call once per run."""
    mpl.rcParams.update(_RC)


def _save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"    wrote {path.name}")
    return path


def _wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width)) if text else ""


def _titleblock(fig: plt.Figure, title: str, subtitle: str = "") -> None:
    """Figure title top-left, provenance caption bottom-left.

    The caption sits just below the figure box; ``savefig.bbox = "tight"``
    expands the saved image to include it, so it never collides with the axes.
    """
    fig.suptitle(title, fontsize=11.5, fontweight="600", x=0.01, ha="left")
    if subtitle:
        fig.text(0.01, -0.004, subtitle, ha="left", va="top", fontsize=7.5,
                 color=INK_MUTED, linespacing=1.6)


def _annotate(ax: plt.Axes, text: str, loc: str = "upper left") -> None:
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(
        x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8,
        color=INK_SECONDARY, linespacing=1.45,
        bbox=dict(facecolor=SURFACE, edgecolor=GRID, linewidth=0.7,
                  boxstyle="round,pad=0.45", alpha=0.92),
    )


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #


def plot_maps(
    als_grid: np.ndarray,
    fh_grid: np.ndarray,
    extent: tuple[float, float, float, float],
    out_path: Path,
    als_label: str,
    subtitle: str = "",
    prod_name: str = "BIOMASS L2A forest height",
    prod_short: str = "BIOMASS",
    grid_name: str = "BIOMASS",
) -> Path:
    """Side-by-side height maps and their difference, on the product grid."""
    residual = fh_grid - als_grid
    both = np.concatenate(
        [als_grid[np.isfinite(als_grid)].ravel(), fh_grid[np.isfinite(fh_grid)].ravel()]
    )
    vmin, vmax = (np.percentile(both, [1, 99]) if both.size else (0, 50))
    rlim = float(np.nanpercentile(np.abs(residual), 98)) if np.isfinite(residual).any() else 1.0
    rlim = max(rlim, 1e-3)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), constrained_layout=True)
    panels = [
        (als_grid, f"ALS reference ({als_label})", CMAP_HEIGHT, dict(vmin=vmin, vmax=vmax), "height (m)"),
        (fh_grid, prod_name, CMAP_HEIGHT, dict(vmin=vmin, vmax=vmax), "height (m)"),
        (residual, f"Difference ({prod_short} - ALS)", CMAP_RESIDUAL,
         dict(norm=TwoSlopeNorm(vcenter=0.0, vmin=-rlim, vmax=rlim)), "difference (m)"),
    ]
    for ax, (data, name, cmap, kw, cbar_label) in zip(axes, panels):
        im = ax.imshow(data, extent=extent, origin="upper", cmap=cmap,
                       interpolation="nearest", **kw)
        ax.set_title(name)
        ax.set_xlabel("longitude (deg)" if abs(extent[1] - extent[0]) < 5 else "easting (m)")
        ax.grid(False)
        ax.tick_params(labelsize=7)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label(cbar_label, size=8, color=INK_SECONDARY)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7, color=AXIS)
    axes[0].set_ylabel("latitude (deg)" if abs(extent[1] - extent[0]) < 5 else "northing (m)")
    _titleblock(fig, f"Canopy height on the {grid_name} grid", subtitle)
    return _save(fig, out_path)


def plot_scatter(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    ref_label: str,
    prod_label: str = "BIOMASS L2A forest height (m)",
    subtitle: str = "",
    prod_short: str = "BIOMASS",
) -> Path:
    """Density scatter against the 1:1 line, with OLS and RMA fits."""
    ok = np.isfinite(reference) & np.isfinite(product)
    x, y = np.asarray(reference)[ok], np.asarray(product)[ok]
    m = compute_metrics(x, y)

    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.05 * (hi - lo)
    lims = (lo - pad, hi + pad)

    fig, ax = plt.subplots(figsize=(6.4, 6.2), constrained_layout=True)
    hb = ax.hexbin(x, y, gridsize=42, extent=(*lims, *lims), mincnt=1,
                   cmap=CMAP_HEIGHT, linewidths=0.0)
    cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("cells per bin", size=8, color=INK_SECONDARY)
    cb.outline.set_visible(False)

    grid = np.array(lims)
    ax.plot(grid, grid, color=INK_MUTED, lw=1.4, ls="--", label="1:1", zorder=3)
    ax.plot(grid, m["ols_slope"] * grid + m["ols_intercept"], color=ORANGE, lw=2.0,
            label=f"OLS  y = {m['ols_slope']:.2f}x {m['ols_intercept']:+.1f}", zorder=4)
    ax.plot(grid, m["rma_slope"] * grid + m["rma_intercept"], color=AQUA, lw=2.0,
            label=f"RMA  y = {m['rma_slope']:.2f}x {m['rma_intercept']:+.1f}", zorder=4)

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.set_xlabel(ref_label)
    ax.set_ylabel(prod_label)
    _titleblock(fig, f"{prod_short} forest height vs ALS canopy height", subtitle)
    ax.legend(loc="lower right")
    _annotate(
        ax,
        "\n".join([
            f"n = {m['n']:.0f} cells",
            f"bias = {m['bias']:+.2f} m",
            f"MAE = {m['mae']:.2f} m",
            f"RMSE = {m['rmse']:.2f} m ({m['rmse_pct']:.0f} %)",
            f"r = {m['pearson_r']:.3f}   CCC = {m['ccc']:.3f}",
            f"R2 vs 1:1 = {m['r2_one_to_one']:.3f}",
        ]),
    )
    return _save(fig, out_path)


def plot_residuals(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    ref_label: str,
    bins: list[float] | None = None,
    subtitle: str = "",
    prod_short: str = "BIOMASS",
) -> Path:
    """Residual against reference height, with binned mean and spread."""
    bins = bins or config.HEIGHT_BINS
    ok = np.isfinite(reference) & np.isfinite(product)
    x, y = np.asarray(reference)[ok], np.asarray(product)[ok]
    resid = y - x

    df = pd.DataFrame({"ref": x, "resid": resid})
    df["bin"] = pd.cut(df["ref"], bins=bins, right=False)
    agg = df.groupby("bin", observed=True)["resid"].agg(["mean", "std", "count"])
    agg = agg[agg["count"] >= 3]
    centres = np.array([(iv.left + iv.right) / 2.0 for iv in agg.index])

    fig, ax = plt.subplots(figsize=(7.6, 5.0), constrained_layout=True)
    ax.axhline(0.0, color=INK_MUTED, lw=1.2, ls="--", zorder=2)
    ax.scatter(x, resid, s=7, color=BLUE, alpha=0.22, linewidths=0, zorder=3,
               label=f"{prod_short} cell")
    if len(agg):
        ax.errorbar(centres, agg["mean"], yerr=agg["std"], color=ORANGE, lw=2.0,
                    marker="o", ms=6, capsize=4, elinewidth=1.6, zorder=5,
                    markeredgecolor=SURFACE, markeredgewidth=1.2,
                    label="binned mean +/- 1 SD")
    ax.set_xlabel(ref_label)
    ax.set_ylabel(f"residual, {prod_short} - ALS (m)")
    _titleblock(fig, "Residual structure across the canopy height range", subtitle)
    legend = ax.legend(loc="upper right", markerscale=3)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
    _annotate(
        ax,
        f"positive = {prod_short} taller than ALS\n"
        "a downward trend indicates saturation over tall canopy",
        loc="lower left",
    )
    return _save(fig, out_path)


def plot_distributions(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    ref_label: str,
    subtitle: str = "",
    prod_short: str = "BIOMASS",
    prod_name: str = "BIOMASS L2A FH",
) -> Path:
    """Marginal distributions and empirical CDFs of the two height estimates."""
    ok = np.isfinite(reference) & np.isfinite(product)
    x, y = np.asarray(reference)[ok], np.asarray(product)[ok]
    edges = np.histogram_bin_edges(np.concatenate([x, y]), bins=40)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6), constrained_layout=True)

    ax1.hist(x, bins=edges, color=BLUE, alpha=0.55, label=f"ALS ({ref_label})")
    ax1.hist(y, bins=edges, color=ORANGE, alpha=0.55, label=prod_name)
    ax1.set_xlabel("canopy height (m)")
    ax1.set_ylabel(f"number of {prod_short} cells")
    ax1.set_title("Height distributions")
    ax1.legend(loc="upper right")

    for data, colour, name in ((x, BLUE, "ALS"), (y, ORANGE, prod_short)):
        s = np.sort(data)
        ax2.plot(s, np.arange(1, s.size + 1) / s.size, color=colour, label=name)
    ax2.set_xlabel("canopy height (m)")
    ax2.set_ylabel("cumulative fraction of cells")
    ax2.set_title("Empirical cumulative distributions")
    ax2.legend(loc="lower right")
    _annotate(
        ax2,
        f"ALS      mean {x.mean():.1f} m, SD {x.std(ddof=1):.1f} m\n"
        f"{prod_short}  mean {y.mean():.1f} m, SD {y.std(ddof=1):.1f} m",
        loc="upper left",
    )
    _titleblock(fig, "Height distributions of the two estimates", subtitle)
    return _save(fig, out_path)


def plot_bland_altman(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    subtitle: str = "",
    prod_short: str = "BIOMASS",
) -> Path:
    """Bland-Altman agreement plot with bias and 95 % limits of agreement."""
    ok = np.isfinite(reference) & np.isfinite(product)
    x, y = np.asarray(reference)[ok], np.asarray(product)[ok]
    mean_h, diff = (x + y) / 2.0, y - x
    bias, sd = diff.mean(), diff.std(ddof=1)

    fig, ax = plt.subplots(figsize=(7.6, 5.0), constrained_layout=True)
    ax.scatter(mean_h, diff, s=8, color=BLUE, alpha=0.25, linewidths=0)
    for value, colour, name in (
        (bias, ORANGE, f"bias {bias:+.2f} m"),
        (bias + 1.96 * sd, INK_MUTED, f"+1.96 SD  {bias + 1.96 * sd:+.2f} m"),
        (bias - 1.96 * sd, INK_MUTED, f"-1.96 SD  {bias - 1.96 * sd:+.2f} m"),
    ):
        ax.axhline(value, color=colour, lw=1.8,
                   ls="-" if colour == ORANGE else "--", zorder=4)
        ax.annotate(name, xy=(0.995, value), xycoords=("axes fraction", "data"),
                    ha="right", va="bottom", fontsize=8, color=INK_SECONDARY)
    ax.set_xlabel(f"mean of ALS and {prod_short} height (m)")
    ax.set_ylabel(f"difference, {prod_short} - ALS (m)")
    _titleblock(fig, "Bland-Altman agreement", subtitle)
    return _save(fig, out_path)


def plot_shift_surface(shift, out_path: Path, subtitle: str = "") -> Path:
    """Agreement as a function of the tested planimetric offset."""
    fine_x_m, fine_y_m = shift.fine_m
    east = shift.offsets_x * fine_x_m
    north = -shift.offsets_y * fine_y_m

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    for ax, data, cmap, name, label in (
        (ax1, shift.pearson, CMAP_HEIGHT, "Pearson r (higher is better)", "r"),
        (ax2, shift.rmse, CMAP_HEIGHT.reversed(), "RMSE (lower is better)", "RMSE (m)"),
    ):
        im = ax.imshow(
            # Row 0 is the most negative y offset, i.e. the largest northing:
            # origin="upper" puts it at the top of the panel.
            data, origin="upper", cmap=cmap, interpolation="nearest", aspect="auto",
            extent=(east.min(), east.max(), north.min(), north.max()),
        )
        ax.plot(shift.shift_east_m, shift.shift_north_m, marker="+", ms=14, mew=2.4,
                color=RED if shift.accepted else INK_MUTED, zorder=5)
        ax.plot(0, 0, marker="o", ms=8, mfc="none", mew=1.8, color=INK_MUTED, zorder=5)
        ax.set_title(name)
        ax.set_xlabel("ALS sampling window offset, east (m)")
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label(label, size=8, color=INK_SECONDARY)
        cb.outline.set_visible(False)
    ax1.set_ylabel("offset, north (m)")
    verdict = "APPLIED" if shift.accepted else "NOT APPLIED"
    _annotate(
        ax1,
        f"+  best offset ({shift.shift_east_m:+.0f}, {shift.shift_north_m:+.0f}) m"
        f"  -- {verdict}\n"
        f"o  no shift\n"
        f"r {shift.pearson_at_zero:.3f} -> {shift.pearson_at_best:.3f}\n"
        + _wrap(shift.reason, 46),
        loc="lower left",
    )
    _titleblock(fig, "Co-registration shift search", subtitle)
    return _save(fig, out_path)


def plot_stat_comparison(table: pd.DataFrame, out_path: Path,
                         subtitle: str = "",
                         prod_short: str = "BIOMASS") -> Path:
    """Which ALS aggregate best matches the product estimate."""
    if table.empty:
        raise ValueError("empty metrics table")
    order = ["mean", "p50", "p90", "p95", "p99", "max"]
    tbl = table.set_index("als_stat").reindex(
        [s for s in order if s in set(table["als_stat"])]
    ).reset_index()
    pos = np.arange(len(tbl))

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), constrained_layout=True)
    panels = [
        ("rmse", "RMSE (m)", BLUE, "lower is better"),
        ("bias", "bias (m)", ORANGE, "zero is best"),
        ("pearson_r", "Pearson r", AQUA, "higher is better"),
    ]
    for ax, (col, label, colour, hint) in zip(axes, panels):
        values = tbl[col].to_numpy()
        ax.bar(pos, values, width=0.62, color=colour, linewidth=0)
        for p, v in zip(pos, values):
            ax.annotate(f"{v:.2f}", xy=(p, v), xytext=(0, 3 if v >= 0 else -11),
                        textcoords="offset points", ha="center", fontsize=8,
                        color=INK_SECONDARY)
        if (values < 0).any():
            ax.axhline(0, color=AXIS, lw=1.0)
        ax.set_xticks(pos)
        ax.set_xticklabels(tbl["als_stat"])
        ax.set_ylabel(label)
        ax.set_title(f"{label} - {hint}")
        ax.set_xlabel(f"ALS statistic per {prod_short} cell")
        ax.margins(y=0.18)
    _titleblock(fig, "Agreement by choice of ALS aggregate", subtitle)
    return _save(fig, out_path)


def plot_quality_strata(table: pd.DataFrame, out_path: Path,
                        subtitle: str = "") -> Path:
    """RMSE, bias and sample size per BIOMASS quality-layer class."""
    if table.empty:
        raise ValueError("empty quality table")
    labels = [str(v) for v in table["quality"]]
    pos = np.arange(len(table))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.4), constrained_layout=True)
    ax1.bar(pos - 0.2, table["rmse"], width=0.38, color=BLUE, linewidth=0, label="RMSE")
    ax1.bar(pos + 0.2, table["bias"], width=0.38, color=ORANGE, linewidth=0, label="bias")
    ax1.axhline(0, color=AXIS, lw=1.0)
    ax1.set_xticks(pos)
    ax1.set_xticklabels(labels, rotation=45, ha="right")
    ax1.set_ylabel("metres")
    ax1.set_title("Error by BIOMASS quality class")
    ax1.set_xlabel("quality layer value")
    ax1.legend(loc="upper right")

    ax2.bar(pos, table["n_cells"], width=0.62, color=AQUA, linewidth=0)
    ax2.set_xticks(pos)
    ax2.set_xticklabels(labels, rotation=45, ha="right")
    ax2.set_ylabel("number of cells")
    ax2.set_title("Sample size per quality class")
    ax2.set_xlabel("quality layer value")
    _titleblock(fig, "Error and sample size by BIOMASS quality class", subtitle)
    return _save(fig, out_path)
