"""Step 3 -- pool the per-scene results into combined statistics and figures.

Reads every ``paired_cells.csv`` written by ``02_compare.py``, concatenates them
into one sample, and reports the comparison over all sites and scenes at once.
Results go to ``outputs/_combined/``.

Maps are deliberately not produced here: the scenes sit on different grids in
different parts of Africa, so there is no shared space to draw them in. Every
other figure has a pooled counterpart.

Examples
--------
    python 03_combine.py
    python 03_combine.py --group-by scene
    python 03_combine.py --site Luki2025 --site Mbalmayo
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

import config
from bgt import metrics, viz
from bgt import viz_combined as vizc


def load_paired(out_dir: Path, sites: list[str] | None = None) -> pd.DataFrame:
    """Concatenate every per-scene ``paired_cells.csv`` under ``out_dir``."""
    frames = []
    for path in sorted(out_dir.rglob("paired_cells.csv")):
        if config.COMBINED_DIR_NAME in path.parts:
            continue
        frame = pd.read_csv(path)
        if sites and frame["site"].iloc[0] not in sites:
            continue
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def scene_label(row: pd.Series) -> str:
    """Compact axis label: site plus the scene's sensing start date."""
    return f"{row['site']}  {row['biomass_date']}"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=config.OUTPUT_DIR)
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable); default: all")
    p.add_argument("--primary-stat", default=config.PRIMARY_ALS_STAT,
                   choices=config.ALS_STATS,
                   help="ALS aggregate used for the headline results")
    p.add_argument("--group-by", default="site", choices=["site", "scene"],
                   help="what the identity colour encodes (default: %(default)s)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    viz.use_style()

    paired = load_paired(args.out_dir, args.site)
    if paired.empty:
        print(f"No paired_cells.csv found under {args.out_dir}.")
        print("Run 02_compare.py first.")
        return 1

    ref_col = f"als_{args.primary_stat}"
    if ref_col not in paired.columns:
        print(f"Column {ref_col} is missing; re-run 02_compare.py with "
              f"--primary-stat {args.primary_stat}.")
        return 1

    paired["scene_label"] = paired.apply(scene_label, axis=1)
    sites = sorted(paired["site"].unique())
    scenes = sorted(paired["scene"].unique())
    print(f"Pooled {len(paired):,} paired cells from {len(scenes)} scene(s) "
          f"across {len(sites)} site(s): {', '.join(sites)}")

    combined_dir = args.out_dir / config.COMBINED_DIR_NAME
    combined_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Tables
    # ------------------------------------------------------------------ #
    pooled = metrics.compute_metrics(paired[ref_col].to_numpy(),
                                     paired["fh"].to_numpy())
    print()
    print(metrics.format_summary(
        pooled,
        ref_label=f"ALS {args.primary_stat}, all sites pooled",
        prod_label="BIOMASS L2A FH",
    ))

    per_scene = []
    for (site, scene, date), group in paired.groupby(
        ["site", "scene", "biomass_date"], observed=True
    ):
        row = {"site": site, "scene": scene, "biomass_date": date,
               "scene_label": f"{site}  {date}"}
        row.update(metrics.compute_metrics(group[ref_col].to_numpy(),
                                           group["fh"].to_numpy()))
        per_scene.append(row)
    per_scene = pd.DataFrame(per_scene)

    per_site = []
    for site, group in paired.groupby("site", observed=True):
        row = {"site": site, "n_scenes": group["scene"].nunique()}
        row.update(metrics.compute_metrics(group[ref_col].to_numpy(),
                                           group["fh"].to_numpy()))
        per_site.append(row)
    per_site = pd.DataFrame(per_site)

    stat_table = metrics.metrics_table(paired, config.ALS_STATS)
    bin_table = metrics.binned_metrics(paired, ref_col)
    qual_table = metrics.quality_metrics(paired, ref_col)

    pd.DataFrame([{"scope": "all sites pooled", "n_sites": len(sites),
                   "n_scenes": len(scenes), **pooled}]).to_csv(
        combined_dir / "metrics_pooled.csv", index=False)
    per_scene.to_csv(combined_dir / "metrics_by_scene.csv", index=False)
    per_site.to_csv(combined_dir / "metrics_by_site.csv", index=False)
    stat_table.to_csv(combined_dir / "metrics_by_als_stat.csv", index=False)
    bin_table.to_csv(combined_dir / "metrics_by_height_bin.csv", index=False)
    if not qual_table.empty:
        qual_table.to_csv(combined_dir / "metrics_by_quality_class.csv", index=False)
    paired.to_csv(combined_dir / "paired_cells_all.csv", index=False)

    # ------------------------------------------------------------------ #
    # Figures
    # ------------------------------------------------------------------ #
    als_dates = ", ".join(
        f"{s}" for s in sorted(paired["site"].unique())
    )
    date_span = f"{paired['biomass_date'].min()} to {paired['biomass_date'].max()}"
    subtitle = (
        f"Sites: {als_dates}   |   ALS: {config.ALS_CHM} "
        f"({args.primary_stat} per cell)\n"
        f"{len(scenes)} BIOMASS {config.BIOMASS_PRODUCT_TYPE} scenes, "
        f"sensing starts {date_span}   |   {len(paired):,} paired cells pooled"
    )
    ref_label = f"ALS canopy height, {args.primary_stat} per cell (m)"
    group_col = args.group_by

    print("\n  figures:")
    vizc.plot_scatter_grouped(paired, ref_col,
                              combined_dir / "fig02_scatter_combined.png",
                              ref_label=ref_label, group_col=group_col,
                              subtitle=subtitle)
    vizc.plot_residuals_grouped(paired, ref_col,
                                combined_dir / "fig03_residuals_combined.png",
                                ref_label=ref_label, group_col=group_col,
                                bins=config.HEIGHT_BINS, subtitle=subtitle)
    vizc.plot_distributions_grouped(paired, ref_col,
                                    combined_dir / "fig04_distributions_combined.png",
                                    ref_label=ref_label, group_col=group_col,
                                    subtitle=subtitle)
    vizc.plot_bland_altman_grouped(paired, ref_col,
                                   combined_dir / "fig05_bland_altman_combined.png",
                                   group_col=group_col, subtitle=subtitle)
    viz.plot_stat_comparison(stat_table,
                             combined_dir / "fig07_als_statistic_combined.png",
                             subtitle=subtitle)
    if not qual_table.empty:
        viz.plot_quality_strata(qual_table,
                                combined_dir / "fig08_quality_classes_combined.png",
                                subtitle=subtitle)
    vizc.plot_scene_metrics(per_scene, combined_dir / "fig09_per_scene.png",
                            subtitle=subtitle)
    vizc.plot_site_metrics(per_site, combined_dir / "fig10_per_site.png",
                           subtitle=subtitle)

    print("\n" + "=" * 72)
    print("Per site (scenes pooled within each site)")
    print(per_site[["site", "n_scenes", "n", "ref_mean", "prod_mean", "bias",
                    "mae", "rmse", "pearson_r"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("\nPer scene")
    print(per_scene[["scene_label", "n", "bias", "mae", "rmse", "pearson_r"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\nCombined outputs written to {combined_dir}")
    print(f"Run finished {datetime.now():%Y-%m-%d %H:%M}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
