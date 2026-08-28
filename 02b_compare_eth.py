"""Step 2b -- compare the ETH 10 m canopy height product against the ALS reference.

This runs on the **ETH grid**, not the BIOMASS grid: the ALS CHM is aggregated up
to the ~9.3 m ETH cells, which is where that product should be judged. It is one
comparison per site, not per scene, because the ETH map is a single global layer
rather than a set of acquisitions.

The BIOMASS-grid comparison of the same product lives in ``02_compare.py``, which
aggregates ETH up to the ~93 m BIOMASS cells so that the two products can be put
side by side on one grid. The two sets of numbers answer different questions and
are not interchangeable.

Results go to ``outputs_eth/<site>/`` -- a separate tree from ``outputs/``, so
that ``03_combine.py`` does not pool these pairs with the BIOMASS ones.

Examples
--------
    python 02b_compare_eth.py
    python 02b_compare_eth.py --site Luki2025
    python 02b_compare_eth.py --primary-stat p95 --min-coverage 0.8
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

import config
from bgt import als as als_mod
from bgt import coreg, metrics, viz
from bgt import eth as eth_mod


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable)")
    p.add_argument("--eth-dir", type=Path, default=config.ETH_DIR,
                   help="where the ETH clips are (default: %(default)s)")
    p.add_argument("--out-dir", type=Path, default=config.ETH_OUTPUT_DIR,
                   help="destination for results (default: %(default)s)")
    p.add_argument("--chm", default=config.ALS_CHM)
    p.add_argument("--primary-stat", default=config.PRIMARY_ALS_STAT,
                   choices=config.ALS_STATS)
    p.add_argument("--fine-cell", type=float, default=config.FINE_CELL_M,
                   help="target fine cell size in metres (default: %(default)s)")
    p.add_argument("--min-coverage", type=float, default=config.MIN_ALS_COVERAGE)
    p.add_argument("--shift-search", type=float, default=0.0,
                   help="planimetric shift search radius in metres. Off by "
                        "default: Sentinel-2 geolocation is good to well under "
                        "one ETH cell, and the search is expensive on a grid "
                        "with this many cells (default: %(default)s)")
    return p.parse_args(argv)


def select_sites(requested: list[str] | None) -> list[Path]:
    dirs = config.site_dirs()
    if not requested:
        return dirs
    by_name = {config.site_name(d): d for d in dirs}
    missing = [s for s in requested if s not in by_name]
    if missing:
        raise SystemExit(f"unknown site(s): {', '.join(missing)}; "
                         f"available: {', '.join(by_name)}")
    return [by_name[s] for s in requested]


def compare_site(als_dir: Path, args: argparse.Namespace) -> dict | None:
    """Aggregate the ALS onto the ETH grid for one site and score the agreement."""
    site = config.site_name(als_dir)
    print(f"\n=== {site} ===")

    chm = als_mod.load_chm(als_dir, chm_name=args.chm)
    footprint = rasterio.transform.array_bounds(*chm.shape, chm.transform)

    eth = eth_mod.load_eth(
        site, footprint_bounds=footprint, footprint_crs=chm.crs,
        pad_m=args.shift_search + 2 * args.fine_cell, eth_dir=args.eth_dir,
    )
    if eth is None:
        return None

    grid = coreg.build_fine_grid(chm, eth, args.fine_cell, args.shift_search)
    aggregator = coreg.BlockAggregator(grid, eth.shape)

    shift = None
    offset = (0, 0)
    if args.shift_search > 0:
        try:
            shift = coreg.find_shift(
                aggregator, eth.height, search_m=args.shift_search,
                step_m=config.SHIFT_STEP_M, criterion=config.SHIFT_CRITERION,
                min_coverage=args.min_coverage,
            )
            offset = shift.applied_offset
        except ValueError as exc:
            print(f"  shift search not possible ({exc}); using nominal geolocation")

    stats = aggregator.full_stats(*offset, config.ALS_STATS)

    rows, cols = np.indices(eth.shape)
    xs, ys = rasterio.transform.xy(eth.transform, rows.ravel(), cols.ravel())
    df = pd.DataFrame({
        "site": site,
        "product": "ETH_GlobalCanopyHeight_10m_2020",
        "row": rows.ravel(),
        "col": cols.ravel(),
        "x": np.asarray(xs),
        "y": np.asarray(ys),
        "eth": eth.height.ravel(),
        "coverage": stats["coverage"].ravel(),
        "n_fine": stats["n_fine"].ravel(),
    })
    if eth.sd is not None:
        df["eth_sd"] = eth.sd.ravel()
    for stat in config.ALS_STATS:
        df[f"als_{stat}"] = stats[stat].ravel()

    ref_col = f"als_{args.primary_stat}"
    keep = (
        np.isfinite(df["eth"])
        & (df["coverage"] >= args.min_coverage)
        & np.isfinite(df[ref_col])
    )
    paired = df[keep].copy()
    print(f"  {len(paired):,} paired cells at coverage >= {args.min_coverage:.0%} "
          f"(of {int(np.isfinite(df['eth']).sum()):,} valid ETH cells)")
    if len(paired) < 20:
        print("  too few paired cells for meaningful statistics -- skipped")
        return None

    headline = metrics.compute_metrics(paired[ref_col].to_numpy(),
                                       paired["eth"].to_numpy())
    print()
    print(metrics.format_summary(
        headline,
        ref_label=f"ALS {args.primary_stat} ({args.chm})",
        prod_label="ETH canopy height 10 m",
    ))

    out_dir = args.out_dir / site
    out_dir.mkdir(parents=True, exist_ok=True)
    stat_table = metrics.metrics_table(paired, config.ALS_STATS, product_col="eth")
    bin_table = metrics.binned_metrics(paired, ref_col, product_col="eth")
    paired.to_csv(out_dir / "paired_cells_eth.csv", index=False)
    stat_table.to_csv(out_dir / "metrics_by_als_stat.csv", index=False)
    bin_table.to_csv(out_dir / "metrics_by_height_bin.csv", index=False)

    # ---------------------------------------------------------------- #
    # Figures
    # ---------------------------------------------------------------- #
    print("\n  figures:")
    als_grid = np.full(eth.shape, np.nan)
    eth_grid = np.full(eth.shape, np.nan)
    als_grid[paired["row"], paired["col"]] = paired[ref_col]
    eth_grid[paired["row"], paired["col"]] = paired["eth"]
    left, bottom, right, top = rasterio.transform.array_bounds(*eth.shape,
                                                               eth.transform)
    extent = (left, right, bottom, top)

    site_label, als_dates = als_mod.acquisition_label(als_dir)
    px_x, px_y = eth.cell_size_m()
    subtitle = (
        f"Site: {site_label}   |   ALS: {args.chm} ({args.primary_stat} per cell), "
        f"acquired {als_dates}\n"
        f"ETH global canopy height 10 m, representative of 2020   |   "
        f"compared on the ETH grid, ~{px_x:.1f} x {px_y:.1f} m cells"
    )
    ref_label = f"ALS canopy height, {args.primary_stat} per cell (m)"
    prod_label = "ETH canopy height (m)"

    site_lims = viz.height_limits(paired[ref_col].to_numpy(),
                                  paired["eth"].to_numpy())

    viz.plot_maps(als_grid, eth_grid, extent, out_dir / "fig01_maps.png",
                  als_label=args.primary_stat, subtitle=subtitle,
                  prod_name="ETH canopy height 10 m", prod_short="ETH",
                  grid_name="ETH 10 m")
    viz.plot_scatter(paired[ref_col].to_numpy(), paired["eth"].to_numpy(),
                     out_dir / "fig02_scatter.png", ref_label=ref_label,
                     prod_label=prod_label, subtitle=subtitle,
                     prod_short="ETH", lims=site_lims)
    viz.plot_residuals(paired[ref_col].to_numpy(), paired["eth"].to_numpy(),
                       out_dir / "fig03_residuals.png", ref_label=ref_label,
                       bins=config.HEIGHT_BINS, subtitle=subtitle,
                       prod_short="ETH", lims=site_lims)
    viz.plot_distributions(paired[ref_col].to_numpy(), paired["eth"].to_numpy(),
                           out_dir / "fig04_distributions.png",
                           ref_label=args.primary_stat, subtitle=subtitle,
                           prod_short="ETH", prod_name="ETH canopy height")
    viz.plot_bland_altman(paired[ref_col].to_numpy(), paired["eth"].to_numpy(),
                          out_dir / "fig05_bland_altman.png", subtitle=subtitle,
                          prod_short="ETH")
    viz.plot_stat_comparison(stat_table, out_dir / "fig07_als_statistic.png",
                             subtitle=subtitle, prod_short="ETH")

    return {
        "site": site_label,
        "product": "ETH_GlobalCanopyHeight_10m_2020",
        "als_chm": args.chm,
        "als_stat": args.primary_stat,
        "als_dates": als_dates,
        "eth_cell_m_x": px_x,
        "eth_cell_m_y": px_y,
        "n_paired_cells": len(paired),
        "min_coverage": args.min_coverage,
        "shift_applied": bool(shift and shift.accepted),
        **headline,
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    viz.use_style()
    site_dirs = select_sites(args.site)
    if not site_dirs:
        print("No 03_processed_* site folders found.")
        return 1

    rows = [r for r in (compare_site(d, args) for d in site_dirs) if r]
    if not rows:
        print("\nNothing compared. Run 01b_download_eth.py first.")
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary_path = args.out_dir / "summary_eth_by_site.csv"
    summary.to_csv(summary_path, index=False)

    print(f"\n{'=' * 70}\nETH 10 m vs ALS, by site\n{'=' * 70}")
    print(summary[["site", "n_paired_cells", "bias", "mae", "rmse",
                   "pearson_r"]].to_string(index=False))
    print(f"\nWritten to {args.out_dir}  ({datetime.now():%Y-%m-%d %H:%M})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
