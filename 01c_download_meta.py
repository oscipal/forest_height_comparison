"""Step 1c -- clip the Meta global canopy height map v2 to each ALS site.

A third spaceborne height product: ~1.2 m canopy height predicted from Maxar
imagery with a DINOv3 backbone, published openly by Meta on AWS.

The tiles are 32768 x 32768 cloud-optimised GeoTIFFs named by zoom-10
quadkey. As with step 1b, only the window covering
each ALS outline is read over HTTP, and no token or account is needed.

Each tile is a mosaic of imagery from different dates. Alongside the height
clip the script therefore writes ``meta_imagery_dates_<site>.csv``: the share
of the ALS outline lying under each acquisition date, read from the tile's
metadata GeoJSON. Clips land in ``data/meta/<site>/``.

Examples
--------
    python 01c_download_meta.py --site Amacayacu --list-only
    python 01c_download_meta.py --site Amacayacu
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio
import requests
from rasterio.merge import merge
from rasterio.warp import transform_bounds

import config
from bgt import als, meta


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=config.META_DIR,
                   help="destination for the clips (default: %(default)s)")
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable); default: every "
                        "03_processed_* folder found")
    p.add_argument("--buffer", type=float, default=config.META_BUFFER_DEG,
                   help="degrees added around the ALS outline (default: %(default)s)")
    p.add_argument("--list-only", action="store_true",
                   help="report the tiles each site needs, download nothing")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(argv)


def clip_tiles(tiles: list[str], bbox: tuple[float, float, float, float],
               out_path: Path) -> Path | None:
    """Read the lon/lat ``bbox`` out of the remote tiles into ``out_path``.

    The clip stays on the product's own Web Mercator grid: resampling happens
    once, onto the comparison grid, and not here as well.
    """
    sources = []
    try:
        for tile in tiles:
            try:
                sources.append(rasterio.open(f"/vsicurl/{meta.tile_url(tile)}"))
            except rasterio.RasterioIOError as exc:
                print(f"    tile {tile} unreadable, skipped: {exc}")
        if not sources:
            return None
        bounds = transform_bounds("EPSG:4326", sources[0].crs, *bbox)
        data, transform = merge(sources, bounds=bounds)
        profile = sources[0].profile
    finally:
        for src in sources:
            src.close()

    profile.update(driver="GTiff", height=data.shape[1], width=data.shape[2],
                   count=1, transform=transform, compress="deflate",
                   tiled=True, blockxsize=512, blockysize=512)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data[0], 1)
        dst.update_tags(
            source="Meta global canopy height v2, DINOv3 (dinov3_global_chm_v2_ml3)",
            url=config.META_BASE_URL,
            tiles=",".join(tiles),
        )
    size_mb = out_path.stat().st_size / 1024 ** 2
    print(f"    wrote {out_path.name}  {data.shape[2]} x {data.shape[1]} px, "
          f"{size_mb:,.1f} MB")
    return out_path


def imagery_dates(tiles: list[str], als_dir: Path) -> pd.DataFrame:
    """Share of the ALS outline under each imagery acquisition date."""
    outline = als.load_outline(als_dir).to_crs("EPSG:4326")
    # Areas in the site's UTM zone, so the shares are shares of ground rather
    # than of square degrees.
    outline = outline.to_crs(outline.estimate_utm_crs())
    total = float(outline.area.sum())

    parts = []
    for tile in tiles:
        response = requests.get(meta.metadata_url(tile), timeout=120)
        response.raise_for_status()
        polygons = gpd.GeoDataFrame.from_features(response.json(), crs="EPSG:4326")
        parts.append(gpd.overlay(polygons.to_crs(outline.crs), outline[["geometry"]],
                                 how="intersection"))
    covered = pd.concat(parts, ignore_index=True)
    table = (covered.assign(area_km2=covered.area / 1e6)
             .groupby("acq_date", as_index=False)["area_km2"].sum())
    table["share_of_outline"] = table["area_km2"] * 1e6 / total
    return table.sort_values("acq_date")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    dirs = config.site_dirs()
    by_name = {config.site_name(d): d for d in dirs}
    if args.site:
        missing = [s for s in args.site if s not in by_name]
        if missing:
            raise SystemExit(f"unknown site(s): {', '.join(missing)}; "
                             f"available: {', '.join(by_name)}")
        dirs = [by_name[s] for s in args.site]

    for als_dir in dirs:
        site = config.site_name(als_dir)
        print(f"\n=== {site} ===")
        bbox = als.search_bbox(als_dir, buffer_deg=args.buffer)
        tiles = meta.tiles_for_bbox(bbox)
        print(f"  outline bbox {tuple(round(v, 4) for v in bbox)}")
        print(f"  tile(s): {', '.join(tiles)}")
        if args.list_only:
            for tile in tiles:
                print(f"    {meta.tile_url(tile)}")
            continue

        path = meta.clip_path(site, args.out_dir)
        if path.is_file() and not args.overwrite:
            print(f"    {path.name} exists, skipped (--overwrite to refetch)")
        else:
            clip_tiles(tiles, bbox, path)

        dates = imagery_dates(tiles, als_dir)
        dates.to_csv(meta.dates_path(site, args.out_dir), index=False)
        for row in dates.itertuples():
            print(f"    imagery {row.acq_date}: {row.share_of_outline:6.1%} "
                  f"of the outline")

    if not args.list_only:
        print(f"\nClips are in {args.out_dir}")
        print("Next: python 02c_compare_meta.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
