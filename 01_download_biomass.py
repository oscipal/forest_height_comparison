"""Step 1 -- find and download the BIOMASS L2A forest height products covering
each ALS site.

Every ``03_processed_<site>`` folder in the repository is treated as a site. The
search box is derived from that site's scan outline, so no coordinates need to
be typed in. Discovery is anonymous; downloading requires a MAAP offline token.

Products land in ``data/biomass/<site>/<product-id>/``.

Examples
--------
    python 01_download_biomass.py --list-only
    python 01_download_biomass.py --token <offline-token>
    python 01_download_biomass.py --site Luki2025 --site Mbalmayo
    python 01_download_biomass.py --datetime 2026-01-01/2026-12-31
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import config
from bgt import als, maap


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--token", help="MAAP offline token (else MAAP_OFFLINE_TOKEN "
                                   "or a .maap_token file is used)")
    p.add_argument("--out-dir", type=Path, default=config.BIOMASS_DIR,
                   help="download destination (default: %(default)s)")
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable); default: every "
                        "03_processed_* folder found")
    p.add_argument("--collection", default=config.BIOMASS_COLLECTION)
    p.add_argument("--product-type", default=config.BIOMASS_PRODUCT_TYPE,
                   help="set to 'any' to keep every product type in the collection")
    p.add_argument("--datetime", default=None,
                   help="STAC datetime filter, e.g. 2026-01-01/2026-12-31")
    p.add_argument("--buffer", type=float, default=config.SEARCH_BUFFER_DEG,
                   help="degrees added around the ALS outline (default: %(default)s)")
    p.add_argument("--list-only", action="store_true",
                   help="search and report, download nothing (needs no token)")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--max-items", type=int, default=None,
                   help="download at most this many items per site (newest last)")
    return p.parse_args(argv)


def select_sites(requested: list[str] | None) -> list[Path]:
    """Resolve ``--site`` values against the discovered site folders."""
    dirs = config.site_dirs()
    if not requested:
        return dirs
    wanted = {s.lower() for s in requested}
    return [d for d in dirs if config.site_name(d).lower() in wanted]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    site_dirs = select_sites(args.site)
    if not site_dirs:
        available = [config.site_name(d) for d in config.site_dirs()]
        print(f"No site folder matches {args.site}. Available: {available}")
        return 1
    print(f"Sites: {', '.join(config.site_name(d) for d in site_dirs)}")

    product_type = None if args.product_type == "any" else args.product_type

    tokens = None
    if not args.list_only:
        tokens = maap.TokenManager(maap.find_offline_token(args.token))
        tokens.refresh()  # fail fast on a bad token, before any large transfer
        print("MAAP access token obtained.")

    all_rows: list[dict] = []
    for als_dir in site_dirs:
        site = config.site_name(als_dir)
        print()
        print("=" * 72)
        print(site)
        print("=" * 72)

        bbox = als.search_bbox(als_dir, args.buffer)
        print(f"  search box (WGS84): "
              f"{bbox[0]:.4f}, {bbox[1]:.4f}, {bbox[2]:.4f}, {bbox[3]:.4f}")

        items = maap.search_items(
            bbox=bbox,
            collection=args.collection,
            product_type=product_type,
            datetime_range=args.datetime,
        )
        if not items:
            print("  no matching products; widen --buffer or relax --datetime")
            continue

        summaries = [{"site": site, **maap.item_summary(it)} for it in items]
        all_rows.extend(summaries)
        table = pd.DataFrame(summaries)
        print(f"  {len(items)} matching product(s):")
        print(table[["title", "datetime", "mission_phase", "grid_code"]]
              .to_string(index=False))

        site_dir = args.out_dir / site
        site_dir.mkdir(parents=True, exist_ok=True)
        table.to_csv(site_dir / "search_results.csv", index=False)
        (site_dir / "search_results.json").write_text(
            json.dumps(summaries, indent=2), encoding="utf-8"
        )

        if args.list_only:
            continue

        chosen = items[-args.max_items:] if args.max_items else items
        for i, item in enumerate(chosen, 1):
            title = item["properties"].get("title", item["id"])
            print(f"\n  [{i}/{len(chosen)}] {title}")
            maap.download_item(item, site_dir, tokens, overwrite=args.overwrite)

    if not all_rows:
        print("\nNothing found for any site.")
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    combined = args.out_dir / "search_results_all_sites.csv"
    pd.DataFrame(all_rows).to_csv(combined, index=False)
    print(f"\nSearch results for all sites: {combined}")
    if not args.list_only:
        print(f"Products are in {args.out_dir}")
        print("Next: python 02_compare.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
