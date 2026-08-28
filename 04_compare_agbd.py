"""Step 4 -- the ALS biomass map against two P-band observables.

Two scatterplots on the 50 m AGBD grid of one site:

1. the **normalised tomographic HV power** at a chosen height above the ground
   (30 m by default), read from ``data/Tomo/`` and paired cell by cell, since
   those layers are written on the AGBD grid itself;
2. the **HV backscatter of a BIOMASS L1b scene**, geocoded onto the same grid
   from the scene's own geolocation LUT and expressed as gamma-nought in dB.

The L1b scene is picked by acquisition date: the covering scene closest to an
anchor date. Neither input carries an acquisition date of its own -- the AGBD
map is derived from ALS and the tomographic layers arrive as bare rasters -- so
the default anchor is the earliest L2A forest-height acquisition the catalogue
returns near the site, which is only a proxy for the BIOMASS epoch of interest.
**If the acquisition date of the tomographic stack is known, pass it as
``--l1b-date``**; ``--l1b-scene`` overrides the choice entirely, and
``--list-scenes`` shows what there is to choose from.

Only the window covering the AGBD footprint is pulled out of the 128 MB
measurement file, so a run transfers roughly 50 MB per scene, not 130.

Examples
--------
    python 04_compare_agbd.py --site Amacayacu
    python 04_compare_agbd.py --site Amacayacu --tomo-height 20
    python 04_compare_agbd.py --site Amacayacu --list-scenes
    python 04_compare_agbd.py --site Amacayacu --skip-l1b
"""

from __future__ import annotations

import argparse
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform_bounds
from shapely.geometry import box

import config
from bgt import agbd, l1b, maap, metrics, viz

#: Width, in t/ha, of the biomass bins the median profile is taken over.
AGBD_BIN_WIDTH = 25.0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site", default="Amacayacu",
                   help="site whose AGBD map is compared (default: %(default)s)")
    p.add_argument("--tomo-height", type=int, default=config.TOMO_HEIGHT_M,
                   help="height above ground, in metres, of the tomographic "
                        "layer (default: %(default)s)")
    p.add_argument("--max-masked", type=float, default=config.AGBD_MAX_MASKED,
                   help="drop cells whose ALS mask fraction exceeds this "
                        "(default: %(default)s)")
    p.add_argument("--polarisation", default=config.L1B_POLARISATION,
                   help="L1b channel to compare (default: %(default)s)")
    p.add_argument("--radiometry", default=config.L1B_RADIOMETRY,
                   choices=["gammaNought", "sigmaNought"],
                   help="radiometric convention (default: %(default)s)")
    p.add_argument("--l1b-scene", default=None,
                   help="use this L1b product id instead of the nearest in time")
    p.add_argument("--l1b-date", default=None,
                   help="anchor date for the scene choice, YYYY-MM-DD; default "
                        "is the site's earliest L2A forest-height acquisition")
    p.add_argument("--list-scenes", action="store_true",
                   help="list the L1b scenes covering the site and stop")
    p.add_argument("--skip-l1b", action="store_true",
                   help="draw only the tomographic figure; needs no token")
    p.add_argument("--token", help="MAAP offline token (else MAAP_OFFLINE_TOKEN "
                                   "or a .maap_token file is used)")
    p.add_argument("--out-dir", type=Path, default=config.AGBD_OUTPUT_DIR)
    return p.parse_args(argv)


def footprint_wgs84(grid: agbd.AgbdGrid,
                    buffer_deg: float = config.L1B_BUFFER_DEG):
    """The AGBD grid outline in lon/lat, and the same box grown by a buffer."""
    with rasterio.open(grid.path) as src:
        bounds = transform_bounds(src.crs, "EPSG:4326", *src.bounds)
    search = (bounds[0] - buffer_deg, bounds[1] - buffer_deg,
              bounds[2] + buffer_deg, bounds[3] + buffer_deg)
    return box(*bounds), search


def write_grid(values: np.ndarray, grid: agbd.AgbdGrid, path: Path,
               description: str) -> Path:
    """Save a layer computed on the AGBD grid, so it can be mapped or re-checked.

    The geocoded backscatter is the one product of this step that is not simply
    read back out of an input, so it is worth keeping: it is what a shift test
    against the AGBD map is run on. It lives beside the annotation and LUT it
    was derived from, under ``data/``, rather than in the output tree -- like
    every other raster in this repository it is an intermediate, and only
    figures and tables are outputs.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff", "height": grid.shape[0], "width": grid.shape[1],
        "count": 1, "dtype": "float32", "crs": grid.crs,
        "transform": grid.transform, "nodata": np.nan, "compress": "deflate",
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(values.astype("float32"), 1)
        dst.set_band_description(1, description)
    print(f"    wrote {path.name}")
    return path


def anchor_date(site: str, explicit: str | None) -> datetime:
    """When the L1b scene should be closest to.

    Neither the AGBD map nor the tomographic layers carry an acquisition date,
    so without ``--l1b-date`` the anchor falls back to the earliest L2A
    forest-height acquisition in the search table step 1 wrote for this site.
    That is a proxy, not a match: the search box is the ALS outline plus a
    buffer, so a scene can appear in the table while its own footprint stops
    short of the site. It fixes the epoch to the right campaign, and the exact
    date should be given explicitly whenever it is known.
    """
    if explicit:
        return datetime.fromisoformat(explicit).replace(tzinfo=timezone.utc)

    table = config.BIOMASS_DIR / site / "search_results.csv"
    if not table.is_file():
        raise SystemExit(
            f"No L2A search table at {table} to date the comparison from. Run "
            "01_download_biomass.py --list-only first, or pass --l1b-date."
        )
    dates = pd.read_csv(table)["datetime"]
    return min(l1b.parse_datetime(d) for d in dates)


def tomographic_figure(grid: agbd.AgbdGrid, args: argparse.Namespace,
                       out_dir: Path) -> pd.DataFrame:
    """Scatter of the tomographic HV layer against AGBD, and the paired cells."""
    with warnings.catch_warnings():
        # The tomographic layers carry no geotransform by design; they are
        # paired with the AGBD map by pixel index.
        warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)
        tomo = agbd.load_tomo(args.tomo_height, grid.shape)

    keep = grid.valid(args.max_masked) & np.isfinite(tomo)
    rows, cols = np.where(keep)
    paired = pd.DataFrame({
        "row": rows,
        "col": cols,
        "agbd_t_ha": grid.agbd[keep],
        "masked_fraction": grid.masked_fraction[keep],
        "tomo_hv": tomo[keep],
    })

    print(f"  {len(paired):,} cells with both AGBD and a tomographic value")
    viz.plot_relation(
        paired["agbd_t_ha"].to_numpy(),
        paired["tomo_hv"].to_numpy(),
        out_dir / f"scatter_tomo_hv_{args.tomo_height}m.png",
        x_label="ALS above-ground biomass density (t/ha)",
        y_label=f"normalised tomographic HV power at {args.tomo_height} m"
                f"\n (fraction of full scale)",
        title=f"Tomographic HV power at {args.tomo_height} m vs ALS biomass",
        subtitle=(
            f"{grid.path.name} band {config.AGBD_BAND_NAME} against "
            f"{agbd.tomo_path(args.tomo_height).name}, paired on the 50 m grid "
            f"they share.\nCells with more than {args.max_masked:.0%} of their "
            "ALS area masked are excluded."
        ),
        x_bin_width=AGBD_BIN_WIDTH,
    )
    return paired


def l1b_figure(grid: agbd.AgbdGrid, args: argparse.Namespace, out_dir: Path,
               tokens: maap.TokenManager) -> pd.DataFrame | None:
    """Scatter of L1b backscatter against AGBD, and the paired cells."""
    outline, search_box = footprint_wgs84(grid)

    print("  searching the L1b catalogue")
    items = l1b.covering_scenes(l1b.search_scenes(search_box), outline)
    if not items:
        raise SystemExit("No L1b scene fully covers the AGBD footprint.")

    table = pd.DataFrame(l1b.scene_summary(it) for it in items)
    if args.list_scenes:
        print(table.to_string(index=False))
        return None

    if args.l1b_scene:
        chosen = next((it for it in items if it["id"] == args.l1b_scene), None)
        if chosen is None:
            raise SystemExit(
                f"{args.l1b_scene} is not among the covering scenes; "
                "run with --list-scenes to see them."
            )
    else:
        when = anchor_date(args.site, args.l1b_date)
        chosen = l1b.nearest_in_time(items, when)
        gap = abs(l1b.parse_datetime(chosen["properties"]["datetime"]) - when)
        source = "given" if args.l1b_date else "from the L2A search table"
        print(f"  anchor date {when:%Y-%m-%d} ({source}); nearest covering "
              f"scene is {gap.days} days away")

    print(f"  scene {chosen['id']}")
    annotation_path, lut_path = l1b.fetch_annotation(chosen, tokens)
    geometry = l1b.read_geometry(annotation_path)
    luts = l1b.read_luts(lut_path, args.radiometry)

    print("  reading the measurement window")
    patch = l1b.read_window(
        l1b.asset_href(chosen, "enclosure_tiff"), geometry, luts, search_box,
        tokens, polarisation=args.polarisation, radiometry=args.radiometry,
    )
    print(f"    window {patch.window.height} x {patch.window.width} pixels")

    power, count = l1b.to_grid(patch, grid.transform, grid.crs, grid.shape)
    backscatter_db = l1b.to_db(power)

    write_grid(
        backscatter_db, grid,
        config.L1B_DIR / chosen["id"] / f"{args.site}_{args.polarisation.lower()}_db.tif",
        description=f"{args.polarisation} {args.radiometry} (dB)",
    )

    keep = grid.valid(args.max_masked) & np.isfinite(backscatter_db)
    rows, cols = np.where(keep)
    label = f"{args.polarisation} {args.radiometry.replace('Nought', '0')}"
    paired = pd.DataFrame({
        "row": rows,
        "col": cols,
        "agbd_t_ha": grid.agbd[keep],
        "masked_fraction": grid.masked_fraction[keep],
        "backscatter_db": backscatter_db[keep],
        "l1b_pixels": count[keep],
    })

    print(f"  {len(paired):,} cells with both AGBD and backscatter "
          f"({paired['l1b_pixels'].mean():.1f} L1b pixels per cell)")
    acquired = l1b.parse_datetime(chosen["properties"]["datetime"])
    viz.plot_relation(
        paired["agbd_t_ha"].to_numpy(),
        paired["backscatter_db"].to_numpy(),
        out_dir / f"scatter_l1b_{args.polarisation.lower()}.png",
        x_label="ALS above-ground biomass density (t/ha)",
        y_label=f"BIOMASS L1b {label} (dB)",
        title=f"BIOMASS L1b {label} vs ALS biomass",
        subtitle=(
            f"{chosen['id']}\nacquired {acquired:%Y-%m-%d}, "
            f"{chosen['properties'].get('sat:orbit_state')}, "
            f"{chosen['properties'].get('grid:code')}. Beta-nought amplitude "
            "squared and converted with the product LUT, averaged in linear "
            f"power over each 50 m cell.\nCells with more than "
            f"{args.max_masked:.0%} of their ALS area masked are excluded."
        ),
        x_bin_width=AGBD_BIN_WIDTH,
    )

    table.to_csv(out_dir / "l1b_scenes.csv", index=False)
    return paired


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    viz.use_style()

    grid = agbd.load_agbd(args.site)
    out_dir = args.out_dir / args.site
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"{args.site}: {grid.shape[0]} x {grid.shape[1]} cells from {grid.path.name}")

    tomo_pairs = tomographic_figure(grid, args, out_dir)
    l1b_pairs = None

    if not args.skip_l1b:
        tokens = maap.TokenManager(maap.find_offline_token(args.token))
        tokens.refresh()  # fail fast on a bad token, before any transfer
        l1b_pairs = l1b_figure(grid, args, out_dir, tokens)
        if l1b_pairs is None:  # --list-scenes
            return 0

    paired = tomo_pairs
    if l1b_pairs is not None:
        paired = tomo_pairs.merge(
            l1b_pairs.drop(columns=["agbd_t_ha", "masked_fraction"]),
            on=["row", "col"], how="outer",
        )
    suffix = f"_{args.tomo_height}m"
    paired.to_csv(out_dir / f"paired_cells_agbd{suffix}.csv", index=False)
    print(f"    wrote paired_cells_agbd{suffix}.csv ({len(paired):,} rows)")

    rows = [{"observable": f"tomographic HV at {args.tomo_height} m",
             **metrics.association(tomo_pairs["agbd_t_ha"],
                                   tomo_pairs["tomo_hv"])}]
    if l1b_pairs is not None:
        rows.append({"observable": f"L1b {args.polarisation} {args.radiometry}",
                     **metrics.association(l1b_pairs["agbd_t_ha"],
                                           l1b_pairs["backscatter_db"])})
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / f"summary_agbd{suffix}.csv", index=False)
    print()
    print(summary[["observable", "n", "pearson_r", "spearman_rho",
                   "ols_slope", "r2_fit"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
