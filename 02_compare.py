"""Step 2 -- co-register each ALS canopy height model with the BIOMASS L2A
forest height products covering it, and quantify their agreement.

For every site, and every BIOMASS product downloaded for that site, the script

1. reads the product over the ALS footprint,
2. aggregates the masked 1 m ALS CHM onto the BIOMASS cells (the coarse grid is
   never resampled -- see ``bgt/coreg.py``),
3. searches the planimetric shift that best aligns the two,
4. writes per-cell pairs, metric tables and figures to
   ``outputs/<site>/<product>/``.

Run ``03_combine.py`` afterwards for the pooled, all-scene results.

Examples
--------
    python 02_compare.py
    python 02_compare.py --site Luki2025
    python 02_compare.py --no-shift-search --primary-stat mean
    python 02_compare.py --quality-max none
"""

from __future__ import annotations

import argparse
import gc
import re
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


# --------------------------------------------------------------------------- #
# Product discovery
# --------------------------------------------------------------------------- #


def _suffixed(out_dir: Path, suffix: str):
    """Return a namer that inserts ``suffix`` before each output's extension.

    A suffixed run writes a complete parallel set of figures and tables into the
    same folders, so a variant of the comparison (a different quality filter,
    say) can sit next to the default one instead of overwriting it.
    """
    def name(filename: str) -> Path:
        stem, dot, ext = filename.partition(".")
        return out_dir / f"{stem}{suffix}{dot}{ext}"
    return name


def _sensing_dates(item_id: str) -> tuple[str, str]:
    """``(ISO start date, human-readable sensing window)`` from a product name.

    BIOMASS names carry the first and last acquisition of the stack as
    ``..._L2A_<start>_<end>_...``; an L2A forest height product is derived from
    several passes, so both dates matter.
    """
    match = re.search(r"_L2A_(\d{8})T\d{6}_(\d{8})T\d{6}_", item_id)
    if not match:
        return "", "sensing dates unknown"
    start, end = (datetime.strptime(g, "%Y%m%d") for g in match.groups())
    window = (
        f"sensed {start:%d %b %Y}"
        if start.date() == end.date()
        else f"sensed {start:%d %b %Y} to {end:%d %b %Y}"
    )
    return f"{start:%Y-%m-%d}", window


def find_products(biomass_dir: Path, site: str) -> list[dict]:
    """Locate a site's FH rasters and their companion quality rasters.

    Products live in ``<biomass_dir>/<site>/<product-id>/``.
    """
    products = []
    for fh_path in sorted((biomass_dir / site).rglob("*_i_fh.tiff")):
        quality = next(iter(fh_path.parent.glob("*_i_quality.tiff")), None)
        annotation = next(iter(fh_path.parent.glob("*_annot.xml")), None)
        # Prefer the product identifier embedded in the file name; fall back to
        # the containing folder when the naming convention does not apply.
        match = re.match(r"(bio_fp_fh__l2a_.*?_t\d{3}_f\d{3})", fh_path.stem.lower())
        products.append(
            {
                "site": site,
                "id": (match.group(1).upper() if match else fh_path.parent.name),
                "fh": fh_path,
                "quality": quality,
                "annotation": annotation,
            }
        )
    return products


# --------------------------------------------------------------------------- #
# One product
# --------------------------------------------------------------------------- #


def compare_product(
    product: dict,
    chm: als_mod.AlsChm,
    chm_secondary: als_mod.AlsChm | None,
    als_dir: Path,
    args: argparse.Namespace,
) -> dict | None:
    """Run the full comparison for a single BIOMASS product."""
    item_id, site = product["id"], product["site"]
    print(f"\n--- {item_id} ---")

    footprint = rasterio.transform.array_bounds(*chm.shape, chm.transform)

    fh = coreg.load_fh(
        fh_path=product["fh"],
        quality_path=product["quality"],
        footprint_bounds=footprint,
        footprint_crs=chm.crs,
        pad_m=args.shift_search + 2 * args.fine_cell,
        quality_max=args.quality_max,
    )
    if fh is None:
        print("  no valid BIOMASS forest height over the ALS footprint -- skipped")
        return None

    px_x, px_y = fh.cell_size_m()
    print(f"  BIOMASS window {fh.shape} cells, ~{px_x:.1f} x {px_y:.1f} m, "
          f"CRS {fh.crs.to_string()}, "
          f"{int(np.isfinite(fh.height).sum())} valid cells")

    grid = coreg.build_fine_grid(chm, fh, args.fine_cell, args.shift_search)
    aggregator = coreg.BlockAggregator(grid, fh.shape)

    if args.no_shift_search:
        shift, offset = None, (0, 0)
        print("  shift search disabled; using the nominal geolocation")
    else:
        try:
            shift = coreg.find_shift(
                aggregator,
                fh.height,
                search_m=args.shift_search,
                step_m=args.shift_step,
                criterion=args.shift_criterion,
                min_coverage=args.min_coverage,
                min_gain=args.shift_min_gain,
                min_cells=args.shift_min_cells,
                force=args.force_shift,
            )
            offset = shift.applied_offset
        except ValueError as exc:
            # Small footprints can leave too few fully covered cells to search
            # on; that is a reason to skip the search, not the comparison.
            print(f"  shift search not possible ({exc}); using nominal geolocation")
            shift, offset = None, (0, 0)

    stats = aggregator.full_stats(*offset, config.ALS_STATS)
    # The per-cell statistics are all that is needed from here on; the fine grid
    # and its summed-area tables are the largest objects in the run, so let them
    # go before the next product builds its own.
    aggregator.release()
    del aggregator, grid
    gc.collect()

    # ------------------------------------------------------------------ #
    # Per-cell table
    # ------------------------------------------------------------------ #
    rows, cols = np.indices(fh.shape)
    xs, ys = rasterio.transform.xy(fh.transform, rows.ravel(), cols.ravel())
    date, sensed = _sensing_dates(item_id)
    df = pd.DataFrame(
        {
            "site": site,
            "scene": item_id,
            "biomass_date": date,
            "row": rows.ravel(),
            "col": cols.ravel(),
            "x": np.asarray(xs),
            "y": np.asarray(ys),
            "fh": fh.height.ravel(),
            "coverage": stats["coverage"].ravel(),
            "n_fine": stats["n_fine"].ravel(),
        }
    )
    if fh.quality is not None:
        df["fh_quality"] = fh.quality.ravel()
    for stat in config.ALS_STATS:
        df[f"als_{stat}"] = stats[stat].ravel()

    if chm_secondary is not None:
        grid2 = coreg.build_fine_grid(chm_secondary, fh, args.fine_cell,
                                      args.shift_search)
        agg2 = coreg.BlockAggregator(grid2, fh.shape)
        stats2 = agg2.full_stats(*offset, [args.primary_stat])
        df[f"als2_{args.primary_stat}"] = stats2[args.primary_stat].ravel()
        del agg2, grid2, stats2
        gc.collect()

    # ------------------------------------------------------------------ #
    # ETH global canopy height, aggregated onto the same BIOMASS cells
    # ------------------------------------------------------------------ #
    #
    # A second spaceborne product on one grid with the first. It is aggregated
    # at the nominal geolocation, not at the offset found for the ALS: that
    # shift describes where the ALS sits relative to BIOMASS and says nothing
    # about where a Sentinel-2 derived product sits.
    eth = None
    if not args.no_eth:
        eth = eth_mod.load_eth(
            site, footprint_bounds=footprint, footprint_crs=chm.crs,
            pad_m=args.shift_search + 2 * args.fine_cell, eth_dir=args.eth_dir,
        )
    if eth is not None:
        eth_grid_fine = coreg.build_fine_grid(eth, fh, args.fine_cell,
                                              args.shift_search)
        eth_agg = coreg.BlockAggregator(eth_grid_fine, fh.shape)
        eth_stats = eth_agg.full_stats(0, 0, config.ALS_STATS)
        df["eth_coverage"] = eth_stats["coverage"].ravel()
        for stat in config.ALS_STATS:
            df[f"eth_{stat}"] = eth_stats[stat].ravel()
        del eth_agg, eth_grid_fine, eth_stats
        gc.collect()

    # Attribute every cell of the window to the first test it fails, so the
    # losses can be reported rather than just implied by the surviving count.
    reason_grid = coreg.exclusion_grid(
        fh, stats["coverage"], stats[args.primary_stat], args.min_coverage
    )
    drop_counts = coreg.exclusion_counts(reason_grid)
    df["reject_reason"] = reason_grid.ravel()

    keep = (
        np.isfinite(df["fh"])
        & (df["coverage"] >= args.min_coverage)
        & np.isfinite(df[f"als_{args.primary_stat}"])
    )
    paired = df[keep].copy()
    print(f"  {len(paired)} paired cells at coverage >= {args.min_coverage:.0%} "
          f"(of {int(np.isfinite(df['fh']).sum())} valid BIOMASS cells)")

    # The funnel must reconcile with the pairing mask, or the map is a fiction.
    if drop_counts["kept"] != len(paired):
        raise AssertionError(
            f"exclusion accounting disagrees with the pairing mask: "
            f"{drop_counts['kept']} vs {len(paired)}"
        )
    in_window = sum(drop_counts.values())
    print("  cells: " + ", ".join(
        f"{coreg.EXCLUSION_LABELS[c]} {drop_counts[coreg.EXCLUSION_KEYS[c]]:,}"
        for c in coreg.EXCLUSION_ORDER
        if drop_counts[coreg.EXCLUSION_KEYS[c]]
    ) + f"  (window {in_window:,})")

    if len(paired) < 20:
        print("  too few paired cells for meaningful statistics -- skipped")
        return None

    # ------------------------------------------------------------------ #
    # Metrics
    # ------------------------------------------------------------------ #
    out_dir = args.out_dir / site / item_id
    out = _suffixed(out_dir, args.suffix)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_col = f"als_{args.primary_stat}"
    headline = metrics.compute_metrics(paired[ref_col].to_numpy(),
                                       paired["fh"].to_numpy())
    print()
    print(metrics.format_summary(
        headline,
        ref_label=f"ALS {args.chm} {args.primary_stat}",
        prod_label="BIOMASS L2A FH",
    ))

    stat_table = metrics.metrics_table(paired, config.ALS_STATS)
    bin_table = metrics.binned_metrics(paired, ref_col)
    qual_table = metrics.quality_metrics(paired, ref_col)

    paired.to_csv(out("paired_cells.csv"), index=False)
    stat_table.to_csv(out("metrics_by_als_stat.csv"), index=False)
    bin_table.to_csv(out("metrics_by_height_bin.csv"), index=False)
    if not qual_table.empty:
        qual_table.to_csv(out("metrics_by_quality_class.csv"), index=False)

    if chm_secondary is not None:
        sec = metrics.compute_metrics(
            paired[f"als2_{args.primary_stat}"].to_numpy(), paired["fh"].to_numpy()
        )
        pd.DataFrame(
            [{"chm": args.chm, **headline}, {"chm": args.chm_secondary, **sec}]
        ).to_csv(out("metrics_by_chm_product.csv"), index=False)

    # ------------------------------------------------------------------ #
    # Three-way comparison on the BIOMASS grid
    # ------------------------------------------------------------------ #
    #
    # Reported as three pairings over one common set of cells, so that the two
    # products are judged on exactly the same sample: neither is credited for
    # covering ground the other one misses.
    eth_ref_col = f"eth_{args.primary_stat}"
    eth_summary: dict[str, float] = {}
    common = None
    if eth is not None and eth_ref_col in paired:
        common = paired[
            np.isfinite(paired[eth_ref_col])
            & (paired["eth_coverage"] >= args.min_coverage)
        ].copy()
        print()
        print(f"  {len(common)} cells carry ALS, BIOMASS and ETH "
              f"(of {len(paired)} ALS/BIOMASS pairs)")
        if len(common) >= 20:
            pairings = [
                ("BIOMASS vs ALS", ref_col, "fh"),
                ("ETH vs ALS", ref_col, eth_ref_col),
                ("BIOMASS vs ETH", eth_ref_col, "fh"),
            ]
            rows = []
            for name, ref, prod in pairings:
                m = metrics.compute_metrics(common[ref].to_numpy(),
                                            common[prod].to_numpy())
                rows.append({"comparison": name, "reference": ref,
                             "product": prod, **m})
            three_way = pd.DataFrame(rows)
            three_way.to_csv(out("metrics_three_way.csv"), index=False)
            print(three_way[["comparison", "n", "bias", "mae", "rmse",
                             "pearson_r"]]
                  .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
            eth_summary = {
                "n_three_way_cells": len(common),
                "eth_vs_als_bias": three_way.loc[1, "bias"],
                "eth_vs_als_rmse": three_way.loc[1, "rmse"],
                "eth_vs_als_pearson_r": three_way.loc[1, "pearson_r"],
                "biomass_vs_eth_bias": three_way.loc[2, "bias"],
                "biomass_vs_eth_rmse": three_way.loc[2, "rmse"],
                "biomass_vs_eth_pearson_r": three_way.loc[2, "pearson_r"],
            }
        else:
            print("  too few three-way cells for meaningful statistics")
            common = None

    # ------------------------------------------------------------------ #
    # Figures
    # ------------------------------------------------------------------ #
    print("\n  figures:")
    als_grid = np.full(fh.shape, np.nan)
    fh_grid = np.full(fh.shape, np.nan)
    als_grid[paired["row"], paired["col"]] = paired[ref_col]
    fh_grid[paired["row"], paired["col"]] = paired["fh"]
    left, bottom, right, top = rasterio.transform.array_bounds(*fh.shape, fh.transform)
    extent = (left, right, bottom, top)

    site_label, als_dates = als_mod.acquisition_label(als_dir)
    subtitle = (
        f"Site: {site_label}   |   ALS: {args.chm} ({args.primary_stat} per cell), "
        f"acquired {als_dates}\n"
        f"BIOMASS scene: {item_id}   |   {config.BIOMASS_PRODUCT_TYPE}, {sensed}"
    )
    ref_label = f"ALS canopy height, {args.primary_stat} per cell (m)"

    # One height range for every scatter and residual plot of this scene, so
    # the BIOMASS and ETH panels can be read against each other directly.
    scene_series = [paired[ref_col].to_numpy(), paired["fh"].to_numpy()]
    if common is not None:
        scene_series += [common[eth_ref_col].to_numpy(), common["fh"].to_numpy()]
    scene_lims = viz.height_limits(*scene_series)

    viz.plot_maps(als_grid, fh_grid, extent, out("fig01_maps.png"),
                  als_label=args.primary_stat, subtitle=subtitle)
    viz.plot_scatter(paired[ref_col].to_numpy(), paired["fh"].to_numpy(),
                     out("fig02_scatter.png"), ref_label=ref_label,
                     subtitle=subtitle, lims=scene_lims)
    viz.plot_residuals(paired[ref_col].to_numpy(), paired["fh"].to_numpy(),
                       out("fig03_residuals.png"), ref_label=ref_label,
                       bins=config.HEIGHT_BINS, subtitle=subtitle,
                       lims=scene_lims)
    viz.plot_distributions(paired[ref_col].to_numpy(), paired["fh"].to_numpy(),
                           out("fig04_distributions.png"),
                           ref_label=args.primary_stat, subtitle=subtitle)
    viz.plot_bland_altman(paired[ref_col].to_numpy(), paired["fh"].to_numpy(),
                          out("fig05_bland_altman.png"), subtitle=subtitle)
    if shift is not None:
        viz.plot_shift_surface(shift, out("fig06_shift_search.png"),
                               subtitle=subtitle)
    viz.plot_stat_comparison(stat_table, out("fig07_als_statistic.png"),
                             subtitle=subtitle)
    if not qual_table.empty:
        viz.plot_quality_strata(qual_table, out("fig08_quality_classes.png"),
                                subtitle=subtitle)

    if common is not None:
        eth_label = f"ETH canopy height, {args.primary_stat} per cell (m)"
        viz.plot_scatter(common[eth_ref_col].to_numpy(), common["fh"].to_numpy(),
                         out("fig09_biomass_vs_eth.png"),
                         ref_label=eth_label,
                         prod_label="BIOMASS L2A forest height (m)",
                         subtitle=subtitle, prod_short="BIOMASS",
                         lims=scene_lims)
        viz.plot_scatter(common[ref_col].to_numpy(),
                         common[eth_ref_col].to_numpy(),
                         out("fig10_eth_vs_als.png"), ref_label=ref_label,
                         prod_label="ETH canopy height (m)", subtitle=subtitle,
                         prod_short="ETH", lims=scene_lims)
        eth_map = np.full(fh.shape, np.nan)
        eth_map[common["row"], common["col"]] = common[eth_ref_col]
        viz.plot_maps(als_grid, eth_map, extent, out("fig11_eth_maps.png"),
                      als_label=args.primary_stat, subtitle=subtitle,
                      prod_name="ETH canopy height (aggregated)",
                      prod_short="ETH", grid_name="BIOMASS")

    viz.plot_exclusions(drop_counts, coreg.EXCLUSION_LABELS_BY_KEY,
                        out("fig12_exclusions.png"), subtitle=subtitle)
    viz.plot_mask_map(reason_grid, drop_counts, coreg.EXCLUSION_LABELS_BY_KEY,
                      coreg.EXCLUSION_KEYS, extent,
                      out("fig13_mask_map.png"), subtitle=subtitle)
    pd.DataFrame([{"site": site, "scene": item_id, "biomass_date": date,
                   "cells_in_window": in_window, **drop_counts}]).to_csv(
        out("exclusion_counts.csv"), index=False)

    return {
        "site": site_label,
        "product": item_id,
        "biomass_date": date,
        "als_chm": args.chm,
        "als_stat": args.primary_stat,
        "als_dates": als_dates,
        "biomass_cell_m_x": px_x,
        "biomass_cell_m_y": px_y,
        "n_paired_cells": len(paired),
        "min_coverage": args.min_coverage,
        "shift_applied": bool(shift and shift.accepted),
        "shift_east_m": shift.shift_east_m if (shift and shift.accepted) else 0.0,
        "shift_north_m": shift.shift_north_m if (shift and shift.accepted) else 0.0,
        "shift_best_east_m": shift.shift_east_m if shift else np.nan,
        "shift_best_north_m": shift.shift_north_m if shift else np.nan,
        "shift_verdict": shift.reason if shift else "shift search not run",
        "pearson_before_shift": shift.pearson_at_zero if shift else np.nan,
        "rmse_before_shift": shift.rmse_at_zero if shift else np.nan,
        **headline,
        **eth_summary,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _float_or_none(value: str) -> float | None:
    """argparse type accepting a number or the literal 'none'."""
    if value.strip().lower() in ("none", "off", ""):
        return None
    return float(value)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--biomass-dir", type=Path, default=config.BIOMASS_DIR)
    p.add_argument("--out-dir", type=Path, default=config.OUTPUT_DIR)
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable); default: all")
    p.add_argument("--chm", default=config.ALS_CHM, help="ALS CHM to use as reference")
    p.add_argument("--chm-secondary", default=config.ALS_CHM_SECONDARY,
                   help="second CHM for a sensitivity check ('none' to skip)")
    p.add_argument("--primary-stat", default=config.PRIMARY_ALS_STAT,
                   choices=config.ALS_STATS,
                   help="ALS aggregate used for the headline results")
    p.add_argument("--min-coverage", type=float, default=config.MIN_ALS_COVERAGE)
    p.add_argument("--fine-cell", type=float, default=config.FINE_CELL_M)
    p.add_argument("--shift-search", type=float, default=config.SHIFT_SEARCH_M)
    p.add_argument("--shift-step", type=float, default=config.SHIFT_STEP_M)
    p.add_argument("--shift-criterion", default=config.SHIFT_CRITERION,
                   choices=["pearson", "rmse"])
    p.add_argument("--no-shift-search", action="store_true",
                   help="compare at the nominal geolocation, without co-registration")
    p.add_argument("--shift-min-gain", type=float, default=config.SHIFT_MIN_GAIN,
                   help="minimum correlation gain before a shift is applied "
                        "(default: %(default)s)")
    p.add_argument("--shift-min-cells", type=int, default=config.SHIFT_MIN_CELLS,
                   help="minimum cells in the co-registration set before a shift "
                        "may be applied (default: %(default)s)")
    p.add_argument("--force-shift", action="store_true",
                   help="apply the best offset even when it fails the peak checks")
    p.add_argument("--quality-max", type=_float_or_none,
                   default=config.FH_QUALITY_MAX,
                   help="keep only BIOMASS pixels with quality <= this value; "
                        "'none' disables the filter (default: %(default)s)")
    p.add_argument("--eth-dir", type=Path, default=config.ETH_DIR,
                   help="where the ETH clips are (default: %(default)s)")
    p.add_argument("--no-eth", action="store_true",
                   help="skip the ETH global canopy height comparison")
    p.add_argument("--suffix", default=None,
                   help="append this to every output filename; default: named "
                        "after the quality filter (_q2, _q20, _allquality), so "
                        "a run never overwrites another run's outputs")
    p.add_argument("--strict-mask", action="store_true",
                   help="additionally apply mask_pd04 and mask_steep to the ALS CHM")
    args = p.parse_args(argv)
    if args.suffix is None:
        args.suffix = config.quality_suffix(args.quality_max)
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    viz.use_style()
    if args.chm_secondary and args.chm_secondary == args.chm:
        args.chm_secondary = "none"

    site_dirs = config.site_dirs()
    if args.site:
        wanted = {s.lower() for s in args.site}
        site_dirs = [d for d in site_dirs if config.site_name(d).lower() in wanted]
    if not site_dirs:
        available = [config.site_name(d) for d in config.site_dirs()]
        print(f"No site folder matches {args.site}. Available: {available}")
        return 1

    extra = config.ALS_MASK_EXTRA + (["mask_pd04.tif", "mask_steep.tif"]
                                     if args.strict_mask else [])
    summaries = []

    for als_dir in site_dirs:
        site = config.site_name(als_dir)
        products = find_products(args.biomass_dir, site)
        print()
        print("=" * 72)
        print(f"{site}  --  {len(products)} BIOMASS forest height product(s)")
        print("=" * 72)
        if not products:
            print("  nothing downloaded for this site; run 01_download_biomass.py")
            continue

        chm = als_mod.load_chm(als_dir, args.chm, config.ALS_MASK, extra)
        chm_secondary = None
        if args.chm_secondary and args.chm_secondary.lower() != "none":
            chm_secondary = als_mod.load_chm(
                als_dir, args.chm_secondary, config.ALS_MASK, extra
            )

        for product in products:
            result = compare_product(product, chm, chm_secondary, als_dir, args)
            if result:
                summaries.append(result)

        # Release the site's rasters before loading the next one.
        del chm, chm_secondary

    if not summaries:
        print("\nNo product produced a usable comparison.")
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    combined = pd.DataFrame(summaries)
    combined.to_csv(_suffixed(args.out_dir, args.suffix)("summary_all_products.csv"),
                    index=False)

    print("\n" + "=" * 72)
    print("Summary across all sites and scenes")
    print(combined[["site", "biomass_date", "n_paired_cells", "shift_applied",
                    "bias", "mae", "rmse", "pearson_r"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nOutputs written to {args.out_dir}")
    print("Next: python 03_combine.py")
    print(f"Run finished {datetime.now():%Y-%m-%d %H:%M}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
