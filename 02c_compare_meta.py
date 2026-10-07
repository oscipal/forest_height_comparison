"""Step 2c -- compare the Meta canopy height map v2 against the ALS reference.

The Meta product (``01c_download_meta.py``) is ~1.2 m, as fine as the ALS
itself, so unlike ETH it is not judged on its own grid: at that resolution any
comparison measures crown-level geolocation and shadowing, not height. Both
CHMs are instead aggregated onto one coarser grid laid over the ALS footprint,
in the ALS CRS and aligned to its pixels -- ``--cell`` metres, 10 m by default,
so the result sits at the resolution of the ETH comparison.

Both go through the same machinery: one area-average warp onto a 1 m fine grid
that divides the target grid exactly, then a block reduction into every
statistic in ``config.ALS_STATS``. Each ALS statistic is paired with **the same
statistic of Meta** -- p90 with p90, mean with mean -- since both are canopy
height models sampled the same way. ``metrics_by_als_stat.csv`` and fig07
therefore compare matched aggregates, where for ETH they set one product value
against each ALS candidate.

Per site, in ``outputs_meta/<site>/``, it writes the figures and tables the ETH
step writes (fig01-05, fig07, ``paired_cells_meta.csv``,
``metrics_by_als_stat.csv``, ``metrics_by_height_bin.csv``) and, where there is
an AGBD map, the biomass scatter of step 4 for this product:
``scatter_agbd_meta.png``, ``paired_cells_agbd_meta.csv`` and
``summary_agbd_meta.csv``. That one alone runs on the 50 m AGBD grid, the only
grid the biomass map exists on; the AGBD map plays no part in the height
comparison.

Examples
--------
    python 02c_compare_meta.py --site Amacayacu
    python 02c_compare_meta.py --site Amacayacu --primary-stat p95
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
from bgt import agbd as agbd_mod
from bgt import als as als_mod
from bgt import coreg, meta, metrics, viz

#: Width, in t/ha, of the biomass bins the median profile is taken over.
AGBD_BIN_WIDTH = 25.0

#: Meta statistic the biomass scatter uses. The mean, because it is the height
#: statistic the AGBD map is itself predicted from (its ``zmean_chm`` band).
AGBD_META_STAT = "mean"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable)")
    p.add_argument("--meta-dir", type=Path, default=config.META_DIR,
                   help="where the Meta clips are (default: %(default)s)")
    p.add_argument("--out-dir", type=Path, default=config.META_OUTPUT_DIR,
                   help="destination for results (default: %(default)s)")
    p.add_argument("--chm", default=config.ALS_CHM)
    p.add_argument("--primary-stat", default=config.PRIMARY_ALS_STAT,
                   choices=config.ALS_STATS)
    p.add_argument("--fine-cell", type=float, default=config.FINE_CELL_M,
                   help="target fine cell size in metres (default: %(default)s)")
    p.add_argument("--cell", type=float, default=config.META_CELL_M,
                   help="cell size, in metres, of the grid ALS and Meta are "
                        "compared on (default: %(default)s)")
    p.add_argument("--min-coverage", type=float, default=config.MIN_ALS_COVERAGE)
    p.add_argument("--max-masked", type=float, default=config.AGBD_MAX_MASKED,
                   help="biomass scatter only: drop cells whose ALS mask "
                        "fraction exceeds this (default: %(default)s)")
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


def aggregate(chm: als_mod.AlsChm, target: meta.TargetGrid, fine_cell: float,
              stats: list[str] = config.ALS_STATS) -> dict[str, np.ndarray]:
    """Per-cell ``stats`` of ``chm`` on ``target``, at nominal position."""
    grid = coreg.build_fine_grid(chm, target, fine_cell, shift_search_m=0.0)
    return coreg.BlockAggregator(grid, target.shape).full_stats(0, 0, stats)


def imagery_note(site: str, meta_dir: Path) -> str:
    """The imagery dates under the ALS outline, for the figure captions."""
    path = meta.dates_path(site, meta_dir)
    if not path.is_file():
        return "imagery date unknown"
    dates = pd.read_csv(path).sort_values("share_of_outline", ascending=False)
    return "Maxar imagery " + ", ".join(
        f"{row.acq_date} ({row.share_of_outline:.1%})" for row in dates.itertuples()
    ) + " of the ALS outline"


def matched_stat_table(paired: pd.DataFrame) -> pd.DataFrame:
    """One metrics row per statistic, ALS and Meta aggregated the same way."""
    rows = []
    for stat in config.ALS_STATS:
        row = {"als_stat": stat}
        row.update(metrics.compute_metrics(paired[f"als_{stat}"].to_numpy(),
                                           paired[f"meta_{stat}"].to_numpy()))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("rmse").reset_index(drop=True)


def biomass_figure(site: str, chm: als_mod.AlsChm, meta_chm: als_mod.AlsChm,
                   args: argparse.Namespace, out_dir: Path, dates: str) -> None:
    """Step 4's biomass scatter, with the Meta height as the observable.

    The AGBD map exists only on its own 50 m grid, so both CHMs are aggregated
    onto that grid here, separately from the height comparison. Paired by the
    same rule as the tomographic and L1b scatters -- an AGBD estimate and at
    most ``--max-masked`` of the ALS area masked -- so the three figures stand
    on the same cells.
    """
    agbd = agbd_mod.load_agbd(site)
    target = meta.TargetGrid(agbd.transform, agbd.crs, agbd.shape)
    print(f"  biomass: the 50 m AGBD grid, {target.shape[0]} x {target.shape[1]} cells")
    stats = [AGBD_META_STAT]
    als_mean = aggregate(chm, target, args.fine_cell, stats)[AGBD_META_STAT]
    meta_mean = aggregate(meta_chm, target, args.fine_cell, stats)[AGBD_META_STAT]

    keep = agbd.valid(args.max_masked) & np.isfinite(meta_mean)
    rows, cols = np.where(keep)
    meta_col = f"meta_{AGBD_META_STAT}"
    sample = pd.DataFrame({
        "row": rows,
        "col": cols,
        "agbd_t_ha": agbd.agbd[keep],
        "masked_fraction": agbd.masked_fraction[keep],
        f"als_{AGBD_META_STAT}": als_mean[keep],
        meta_col: meta_mean[keep],
    })
    sample.to_csv(out_dir / "paired_cells_agbd_meta.csv", index=False)
    print(f"  {len(sample):,} cells with both AGBD and a Meta height")
    viz.plot_relation(
        sample["agbd_t_ha"].to_numpy(),
        sample[meta_col].to_numpy(),
        out_dir / "scatter_agbd_meta.png",
        x_label="ALS above-ground biomass density (t/ha)",
        y_label=f"Meta canopy height v2, {AGBD_META_STAT} per cell (m)",
        title="Meta canopy height vs ALS biomass",
        subtitle=(
            f"{agbd.path.name} band {config.AGBD_BAND_NAME} against "
            f"{meta.clip_path(site).name}, aggregated to the "
            f"{AGBD_META_STAT} per cell on the 50 m AGBD grid.\n{dates}. Cells "
            f"with more than {args.max_masked:.0%} of their ALS area masked are "
            "excluded."
        ),
        x_bin_width=AGBD_BIN_WIDTH,
    )
    # The ALS row is the ceiling: the AGBD map is predicted from that very
    # height, so no height product can relate to it more closely.
    pd.DataFrame([
        {"observable": f"Meta CHM v2 {AGBD_META_STAT} height",
         **metrics.association(sample["agbd_t_ha"], sample[meta_col])},
        {"observable": f"ALS {AGBD_META_STAT} height (AGBD predictor)",
         **metrics.association(sample["agbd_t_ha"], sample[f"als_{AGBD_META_STAT}"])},
    ]).to_csv(out_dir / "summary_agbd_meta.csv", index=False)


def compare_site(als_dir: Path, args: argparse.Namespace) -> dict | None:
    """Aggregate ALS and Meta onto one grid for one site and score the agreement."""
    site = config.site_name(als_dir)
    print(f"\n=== {site} ===")

    meta_chm = meta.load_clip(site, args.meta_dir)
    if meta_chm is None:
        return None
    chm = als_mod.load_chm(als_dir, chm_name=args.chm)

    target = meta.footprint_grid(chm, args.cell)
    grid_label = f"a {args.cell:g} m grid over the ALS footprint"
    print(f"  target: {grid_label}, {target.shape[0]} x {target.shape[1]} cells")

    print("  ALS:")
    als_stats = aggregate(chm, target, args.fine_cell)
    print("  Meta:")
    meta_stats = aggregate(meta_chm, target, args.fine_cell)

    rows, cols = np.indices(target.shape)
    xs, ys = rasterio.transform.xy(target.transform, rows.ravel(), cols.ravel())
    df = pd.DataFrame({
        "site": site,
        "product": "Meta_GlobalCanopyHeight_v2_DINOv3",
        "row": rows.ravel(),
        "col": cols.ravel(),
        "x": np.asarray(xs),
        "y": np.asarray(ys),
        "coverage": als_stats["coverage"].ravel(),
        "meta_coverage": meta_stats["coverage"].ravel(),
        "n_fine": als_stats["n_fine"].ravel(),
    })
    for stat in config.ALS_STATS:
        df[f"als_{stat}"] = als_stats[stat].ravel()
        df[f"meta_{stat}"] = meta_stats[stat].ravel()
    # The column the figures and the stratified table read, as "eth" is for ETH.
    df["meta"] = df[f"meta_{args.primary_stat}"]

    ref_col = f"als_{args.primary_stat}"
    keep = (
        np.isfinite(df["meta"])
        & (df["coverage"] >= args.min_coverage)
        & (df["meta_coverage"] >= args.min_coverage)
        & np.isfinite(df[ref_col])
    )
    paired = df[keep].copy()
    print(f"  {len(paired):,} paired cells at coverage >= {args.min_coverage:.0%} "
          f"(of {target.shape[0] * target.shape[1]:,} cells in the grid)")
    if len(paired) < 20:
        print("  too few paired cells for meaningful statistics -- skipped")
        return None

    headline = metrics.compute_metrics(paired[ref_col].to_numpy(),
                                       paired["meta"].to_numpy())
    print()
    print(metrics.format_summary(
        headline,
        ref_label=f"ALS {args.primary_stat} ({args.chm})",
        prod_label=f"Meta CHM v2 {args.primary_stat}",
    ))

    out_dir = args.out_dir / site
    out_dir.mkdir(parents=True, exist_ok=True)
    stat_table = matched_stat_table(paired)
    bin_table = metrics.binned_metrics(paired, ref_col, product_col="meta")
    paired.to_csv(out_dir / "paired_cells_meta.csv", index=False)
    stat_table.to_csv(out_dir / "metrics_by_als_stat.csv", index=False)
    bin_table.to_csv(out_dir / "metrics_by_height_bin.csv", index=False)

    # ---------------------------------------------------------------- #
    # Figures
    # ---------------------------------------------------------------- #
    print("\n  figures:")
    als_grid = np.full(target.shape, np.nan)
    meta_grid = np.full(target.shape, np.nan)
    als_grid[paired["row"], paired["col"]] = paired[ref_col]
    meta_grid[paired["row"], paired["col"]] = paired["meta"]
    left, bottom, right, top = rasterio.transform.array_bounds(*target.shape,
                                                               target.transform)
    extent = (left, right, bottom, top)

    site_label, als_dates = als_mod.acquisition_label(als_dir)
    dates = imagery_note(site, args.meta_dir)
    subtitle = (
        f"Site: {site_label}   |   ALS: {args.chm}, acquired {als_dates}\n"
        f"Meta canopy height v2 (DINOv3, ~1.2 m), {dates}\n"
        f"Both aggregated to {args.primary_stat} per cell on {grid_label}"
    )
    stat_label = args.primary_stat
    ref_label = f"ALS canopy height, {stat_label} per cell (m)"
    prod_label = f"Meta canopy height, {stat_label} per cell (m)"
    lims = viz.height_limits(paired[ref_col].to_numpy(), paired["meta"].to_numpy())

    viz.plot_maps(als_grid, meta_grid, extent, out_dir / "fig01_maps.png",
                  als_label=stat_label, subtitle=subtitle,
                  prod_name=f"Meta canopy height v2 ({stat_label})",
                  prod_short="Meta", grid_name=f"{target.cell_size_m()[0]:g} m")
    viz.plot_scatter(paired[ref_col].to_numpy(), paired["meta"].to_numpy(),
                     out_dir / "fig02_scatter.png", ref_label=ref_label,
                     prod_label=prod_label, subtitle=subtitle,
                     prod_short="Meta", lims=lims)
    viz.plot_residuals(paired[ref_col].to_numpy(), paired["meta"].to_numpy(),
                       out_dir / "fig03_residuals.png", ref_label=ref_label,
                       bins=config.HEIGHT_BINS, subtitle=subtitle,
                       prod_short="Meta", lims=lims)
    viz.plot_distributions(paired[ref_col].to_numpy(), paired["meta"].to_numpy(),
                           out_dir / "fig04_distributions.png",
                           ref_label=stat_label, subtitle=subtitle,
                           prod_short="Meta", prod_name="Meta canopy height")
    viz.plot_bland_altman(paired[ref_col].to_numpy(), paired["meta"].to_numpy(),
                          out_dir / "fig05_bland_altman.png", subtitle=subtitle,
                          prod_short="Meta")
    viz.plot_stat_comparison(
        stat_table, out_dir / "fig07_als_statistic.png",
        subtitle=subtitle + "\nEach bar pairs the ALS statistic with the same "
                            "statistic of Meta", prod_short="Meta")

    if agbd_mod.agbd_path(site).is_file():
        biomass_figure(site, chm, meta_chm, args, out_dir, dates)

    px_x, px_y = target.cell_size_m()
    return {
        "site": site_label,
        "product": "Meta_GlobalCanopyHeight_v2_DINOv3",
        "meta_imagery": dates,
        "als_chm": args.chm,
        "als_stat": args.primary_stat,
        "als_dates": als_dates,
        "cell_m_x": px_x,
        "cell_m_y": px_y,
        "n_paired_cells": len(paired),
        "min_coverage": args.min_coverage,
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
        print("\nNothing compared. Run 01c_download_meta.py first.")
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out_dir / "summary_meta_by_site.csv", index=False)

    print(f"\n{'=' * 70}\nMeta CHM v2 vs ALS, by site\n{'=' * 70}")
    print(summary[["site", "n_paired_cells", "bias", "mae", "rmse",
                   "pearson_r"]].to_string(index=False))
    print(f"\nWritten to {args.out_dir}  ({datetime.now():%Y-%m-%d %H:%M})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
