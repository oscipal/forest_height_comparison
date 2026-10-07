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
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch

import config
from bgt import metrics
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

#: Fraction of the sample left outside each end of the default y range of
#: :func:`plot_relation`.
RELATION_Y_TAIL = 0.005

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
    # Both scales are fixed rather than fitted to the scene: a colour has to
    # mean the same height, and the same error, on every map in the repository.
    vmin, vmax = config.HEIGHT_LIMITS_M
    rlim = float(config.RESIDUAL_LIMIT_M)

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


def height_limits(*arrays: np.ndarray) -> tuple[float, float]:
    """The fixed height-axis range every scatter and residual plot uses.

    Fixed rather than fitted, so figures from different scenes, from the pooled
    run and from the ETH grid are directly comparable. Cells outside the range
    are drawn outside the axes and so are not visible; they remain in every
    metric, since the statistics are computed on the paired sample rather than
    on the plotted subset. The arguments are ignored and kept so that a
    data-driven range can be restored in one place.
    """
    lo, hi = config.HEIGHT_LIMITS_M
    return float(lo), float(hi)


def plot_scatter(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    ref_label: str,
    prod_label: str = "BIOMASS L2A forest height (m)",
    subtitle: str = "",
    prod_short: str = "BIOMASS",
    lims: tuple[float, float] | None = None,
) -> Path:
    """Density scatter against the 1:1 line, with OLS and RMA fits.

    ``lims`` is the shared height range for the comparison this figure belongs
    to; without it the figure falls back to its own data.
    """
    ok = np.isfinite(reference) & np.isfinite(product)
    x, y = np.asarray(reference)[ok], np.asarray(product)[ok]
    m = compute_metrics(x, y)

    lims = lims or height_limits(x, y)

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


#: Most points a point-cloud figure draws. Statistics always use every cell;
#: only the drawing is thinned, since tens of millions of translucent markers
#: render no differently from half a million and take minutes. Every run on the
#: BIOMASS and ETH grids stays below this, so their figures are unaffected.
MAX_DRAWN_POINTS = 500_000


def _drawn(n: int) -> np.ndarray:
    """Indices of the points to draw: all of them, or a fixed random sample."""
    if n <= MAX_DRAWN_POINTS:
        return np.arange(n)
    return np.random.default_rng(0).choice(n, MAX_DRAWN_POINTS, replace=False)


def plot_residuals(
    reference: np.ndarray,
    product: np.ndarray,
    out_path: Path,
    ref_label: str,
    bins: list[float] | None = None,
    subtitle: str = "",
    prod_short: str = "BIOMASS",
    lims: tuple[float, float] | None = None,
) -> Path:
    """Residual against reference height, with binned mean and spread.

    ``lims`` fixes the height axis to the comparison's shared range; the
    residual axis is left free, since its useful span differs by an order of
    magnitude between sites.
    """
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
    drawn = _drawn(x.size)
    ax.scatter(x[drawn], resid[drawn], s=7, color=BLUE, alpha=0.22, linewidths=0,
               zorder=3, label=f"{prod_short} cell")
    if len(agg):
        ax.errorbar(centres, agg["mean"], yerr=agg["std"], color=ORANGE, lw=2.0,
                    marker="o", ms=6, capsize=4, elinewidth=1.6, zorder=5,
                    markeredgecolor=SURFACE, markeredgewidth=1.2,
                    label="binned mean +/- 1 SD")
    ax.set_xlim(lims or height_limits(x))
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
        frac = np.arange(1, s.size + 1) / s.size
        # Evenly spaced vertices of the sorted sample trace the same curve.
        step = max(1, s.size // MAX_DRAWN_POINTS)
        ax2.plot(s[::step], frac[::step], color=colour, label=name)
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
    drawn = _drawn(diff.size)
    ax.scatter(mean_h[drawn], diff[drawn], s=8, color=BLUE, alpha=0.25, linewidths=0)
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


# --------------------------------------------------------------------------- #
# Where the cells go
# --------------------------------------------------------------------------- #
#
# Colour here encodes identity -- which test removed a cell. The four hues
# carrying the real exclusion reasons were validated as an all-pairs set
# (blue / orange / aqua / violet: worst CVD dE 9.2, worst normal-vision dE
# 16.3 on this surface). "Kept" and "outside footprint" deliberately take
# recessive neutrals rather than a fifth and sixth hue: they are the states the
# reader is not asked to discriminate between, and holding the palette to four
# hues is what keeps the four that matter separable. Aqua sits below 3:1
# against the surface, so the legend carries counts and the same numbers appear
# in the funnel chart and the CSV -- identity is never colour alone.

EXCLUSION_COLORS = {
    # "kept" and "outside footprint" are both neutral, but they sit adjacent in
    # the stacked scene chart and are its two largest segments, so they need
    # real separation from each other -- hence a mid grey against a near-white
    # rather than two similar tints.
    "kept": "#bfbeb5",
    "als_absent": "#f2f1ec",
    "fh_nodata": "#2a78d6",
    "fh_out_of_range": "#eb6834",
    "fh_quality_rejected": "#1baf7a",
    "als_coverage_low": "#4a3aa7",
    "als_stat_missing": "#898781",
}


def plot_exclusions(
    counts: dict[str, int],
    labels: dict[str, str],
    out_path: Path,
    subtitle: str = "",
    grid_name: str = "BIOMASS",
) -> Path:
    """Funnel of how the analysis window is reduced to the paired sample.

    Each bar counts the cells removed by one test, applied in order, so the
    bars plus the kept bar sum to the whole window. The right-hand axis states
    what fraction of the cells still standing at that point each test removed,
    which is the number that says whether a test is doing a little or a lot of
    work.
    """
    keys = [k for k in labels if k in counts]
    total = sum(counts.values())
    if total == 0:
        raise ValueError("empty exclusion counts")

    # Share of the cells still remaining when each test is applied.
    remaining, shares = total, []
    for key in keys:
        n = counts[key]
        shares.append(100.0 * n / remaining if remaining else 0.0)
        if key != "kept":
            remaining -= n

    pos = np.arange(len(keys))
    values = [counts[k] for k in keys]
    colors = [EXCLUSION_COLORS.get(k, INK_MUTED) for k in keys]

    fig, ax = plt.subplots(figsize=(9.6, 0.62 * len(keys) + 2.4),
                           constrained_layout=True)
    ax.barh(pos, values, height=0.66, color=colors, linewidth=0.8,
            edgecolor=SURFACE)
    for p, v, share, key in zip(pos, values, shares, keys):
        note = (f"{v:,}  ({share:.0f} % of those remaining)"
                if key != "kept" else f"{v:,}  ({100.0 * v / total:.0f} % of window)")
        ax.annotate(note, xy=(v, p), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK_SECONDARY)
    ax.set_yticks(pos)
    ax.set_yticklabels([labels[k] for k in keys])
    ax.invert_yaxis()
    ax.set_xlabel(f"{grid_name} cells   |   analysis window {total:,}")
    ax.margins(x=0.26)
    ax.grid(axis="y", visible=False)
    # The note goes in the axes title rather than a floating box: the bars fill
    # the plot area from the left, so no corner is reliably free.
    ax.set_title("tests applied top to bottom; each sees only the cells that "
                 "survived the tests above it",
                 fontsize=8.5, fontweight="normal", color=INK_SECONDARY)
    _titleblock(fig, f"Where the {grid_name} cells go", subtitle)
    return _save(fig, out_path)


def plot_mask_map(
    reason: np.ndarray,
    counts: dict[str, int],
    labels: dict[str, str],
    codes: dict[int, str],
    extent: tuple[float, float, float, float],
    out_path: Path,
    subtitle: str = "",
    grid_name: str = "BIOMASS",
) -> Path:
    """Map of the analysis window, each cell coloured by why it was excluded."""
    present = [c for c in codes if (reason == c).any()]
    key_of = {c: codes[c] for c in present}

    # Remap the sparse codes onto 0..n-1 so the colour list indexes directly.
    lookup = np.full(int(max(codes)) + 1, -1, dtype=np.int16)
    for i, c in enumerate(present):
        lookup[c] = i
    indexed = lookup[reason]

    colors = [EXCLUSION_COLORS.get(key_of[c], INK_MUTED) for c in present]
    cmap = ListedColormap(colors)

    fig, ax = plt.subplots(figsize=(8.2, 6.6), constrained_layout=True)
    ax.imshow(indexed, extent=extent, origin="upper", cmap=cmap,
              interpolation="nearest", vmin=-0.5, vmax=len(present) - 0.5)
    ax.set_xlabel("longitude (deg)" if abs(extent[1] - extent[0]) < 5 else "easting (m)")
    ax.set_ylabel("latitude (deg)" if abs(extent[1] - extent[0]) < 5 else "northing (m)")
    ax.grid(False)
    ax.tick_params(labelsize=7)

    total = sum(counts.values())
    handles = [
        Patch(facecolor=EXCLUSION_COLORS.get(key_of[c], INK_MUTED),
              edgecolor=AXIS, linewidth=0.6,
              label=f"{labels[key_of[c]]}  --  {counts.get(key_of[c], 0):,} "
                    f"({100.0 * counts.get(key_of[c], 0) / total:.0f} %)")
        for c in present
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0.0, title="cell status", alignment="left")
    _titleblock(fig, f"Why each {grid_name} cell was kept or dropped", subtitle)
    return _save(fig, out_path)


def _rounded_range(lo: float, hi: float) -> tuple[float, float]:
    """Widen ``(lo, hi)`` outward to round numbers, whatever the unit's scale.

    Rounding to whole numbers is right for a quantity spanning tens -- decibels,
    digital numbers -- and useless for one spanning a fraction of 1, where it
    collapses every range to (0, 1). So the step is taken from the span itself.
    """
    span = hi - lo
    if span >= 5.0:
        return float(np.floor(lo)), float(np.ceil(hi))
    if span <= 0.0:
        return float(lo), float(lo + 1.0)
    step = 10.0 ** np.floor(np.log10(span)) / 2.0
    return float(np.floor(lo / step) * step), float(np.ceil(hi / step) * step)


def plot_relation(
    x: np.ndarray,
    y: np.ndarray,
    out_path: Path,
    x_label: str,
    y_label: str,
    title: str,
    subtitle: str = "",
    x_bin_width: float | None = None,
    y_lims: tuple[float, float] | None = None,
) -> Path:
    """Density scatter of two quantities that do not share a unit.

    The counterpart to :func:`plot_scatter`, which compares two estimates of the
    same thing and so is built around the 1:1 line. Here the units differ, so
    there is no 1:1 line and no bias: what the figure has to show is the shape
    of the relation. The binned median with its interquartile band carries that
    -- it is where saturation becomes visible -- and the OLS line is drawn only
    as the linear reference the correlation refers to.
    """
    ok = np.isfinite(x) & np.isfinite(y)
    a, b = np.asarray(x)[ok], np.asarray(y)[ok]
    m = metrics.association(a, b)

    x_lims = (float(np.floor(a.min())), float(np.ceil(a.max())))
    # A handful of very low cells -- open water, a clearing -- can span more of
    # the y range than the whole population does, leaving the relation squashed
    # into a strip. The default range is robust, and the annotation says how
    # many cells fall outside it so nothing is dropped quietly.
    y_lims = y_lims or _rounded_range(
        float(np.percentile(b, 100.0 * RELATION_Y_TAIL)),
        float(np.percentile(b, 100.0 * (1.0 - RELATION_Y_TAIL))),
    )
    off_scale = int(np.count_nonzero((b < y_lims[0]) | (b > y_lims[1])))

    fig, ax = plt.subplots(figsize=(6.8, 6.0), constrained_layout=True)
    hb = ax.hexbin(a, b, gridsize=44, extent=(*x_lims, *y_lims), mincnt=1,
                   cmap=CMAP_HEIGHT, linewidths=0.0)
    cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("cells per bin", size=8, color=INK_SECONDARY)
    cb.outline.set_visible(False)

    grid = np.array(x_lims)
    ax.plot(grid, m["ols_slope"] * grid + m["ols_intercept"], color=ORANGE, lw=2.0,
            label=f"OLS  y = {m['ols_slope']:.3g}x {m['ols_intercept']:+.3g}", zorder=4)

    width = x_bin_width or (x_lims[1] - x_lims[0]) / 12.0
    edges = np.arange(x_lims[0], x_lims[1] + width, width)
    profile = metrics.binned_profile(a, b, edges)
    if not profile.empty:
        ax.fill_between(profile["x_centre"], profile["y_q25"], profile["y_q75"],
                        color=AQUA, alpha=0.18, linewidth=0, zorder=3)
        ax.plot(profile["x_centre"], profile["y_median"], color=AQUA, lw=2.0,
                marker="o", markersize=3.5,
                label=f"median per {width:g}-wide bin (IQR shaded)", zorder=5)

    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    _titleblock(fig, title, subtitle)
    ax.legend(loc="lower right")
    _annotate(
        ax,
        "\n".join([
            f"n = {m['n']:.0f} cells",
            f"Pearson r = {m['pearson_r']:.3f}   R2 = {m['r2_fit']:.3f}",
            f"Spearman rho = {m['spearman_rho']:.3f}",
            f"scatter about fit = {m['rmse_fit']:.3g}",
        ] + ([f"{off_scale} cells outside the y range"] if off_scale else [])),
    )
    return _save(fig, out_path)
