"""Access to the ESA MAAP catalogue and to token-protected BIOMASS downloads.

Discovery (searching the STAC catalogue) is open. Downloading the actual product
files requires a bearer token, which is derived from a long-lived *offline
token* issued by the ESA MAAP portal:

    https://portal.maap.eo.esa.int  ->  user profile  ->  offline token

Provide it via the ``MAAP_OFFLINE_TOKEN`` environment variable, a ``.maap_token``
file next to this repository, or the ``--token`` command line flag.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import requests

import config

#: Access tokens are short lived; refresh a little before they actually expire.
_TOKEN_MARGIN_S = 30.0


class MaapError(RuntimeError):
    """Raised when the MAAP catalogue or token service misbehaves."""


# --------------------------------------------------------------------------- #
# Authentication
# --------------------------------------------------------------------------- #


def find_offline_token(explicit: str | None = None) -> str:
    """Locate the MAAP offline token.

    Resolution order: explicit argument, ``MAAP_OFFLINE_TOKEN`` environment
    variable, ``.maap_token`` file in the repository root or the home directory.
    """
    if explicit:
        return explicit.strip()

    env = os.environ.get("MAAP_OFFLINE_TOKEN")
    if env:
        return env.strip()

    for candidate in (config.ROOT / ".maap_token", Path.home() / ".maap_token"):
        if candidate.is_file():
            text = candidate.read_text(encoding="utf-8").strip()
            if text:
                return text

    raise MaapError(
        "No MAAP offline token found. Pass --token, set MAAP_OFFLINE_TOKEN, or "
        f"write the token to {config.ROOT / '.maap_token'}."
    )


@dataclass
class TokenManager:
    """Trades an offline token for short-lived access tokens, refreshing as needed.

    BIOMASS scenes take minutes to download and the access token typically lives
    for only five, so every request asks this object for a currently valid token
    rather than caching one at the start.
    """

    offline_token: str
    token_url: str = config.MAAP_TOKEN_URL
    client_id: str = config.MAAP_CLIENT_ID
    client_secret: str = config.MAAP_CLIENT_SECRET

    _access_token: str | None = field(default=None, init=False, repr=False)
    _expires_at: float = field(default=0.0, init=False, repr=False)

    @property
    def access_token(self) -> str:
        if self._access_token is None or time.time() >= self._expires_at:
            self.refresh()
        assert self._access_token is not None
        return self._access_token

    def refresh(self) -> str:
        payload = {
            "client_id": self.client_id,
            "client_secret": self.client_secret,
            "grant_type": "refresh_token",
            "refresh_token": self.offline_token,
            "scope": "offline_access openid",
        }
        response = requests.post(self.token_url, data=payload, timeout=60)
        if not response.ok:
            raise MaapError(
                "Token exchange failed "
                f"({response.status_code}): {response.text[:400]}\n"
                "Check that the offline token is current -- they are revoked when "
                "unused for 30 days."
            )
        body = response.json()
        token = body.get("access_token")
        if not token:
            raise MaapError(f"Token response contained no access_token: {body}")

        self._access_token = token
        self._expires_at = time.time() + float(body.get("expires_in", 300)) - _TOKEN_MARGIN_S
        return token

    def headers(self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.access_token}"}


# --------------------------------------------------------------------------- #
# Catalogue search
# --------------------------------------------------------------------------- #


def search_items(
    bbox: tuple[float, float, float, float],
    collection: str = config.BIOMASS_COLLECTION,
    product_type: str | None = config.BIOMASS_PRODUCT_TYPE,
    datetime_range: str | None = None,
    limit: int = 500,
    stac_url: str = config.MAAP_STAC_URL,
) -> list[dict[str, Any]]:
    """Return STAC items of ``collection`` intersecting ``bbox``.

    ``bbox`` is ``(lon_min, lat_min, lon_max, lat_max)`` in WGS84.
    ``datetime_range`` follows the STAC convention, e.g.
    ``"2026-01-01T00:00:00Z/2026-12-31T23:59:59Z"``.
    """
    params: dict[str, Any] = {
        "bbox": ",".join(f"{v:.6f}" for v in bbox),
        "limit": limit,
    }
    if datetime_range:
        params["datetime"] = datetime_range

    url = f"{stac_url.rstrip('/')}/collections/{collection}/items"
    response = requests.get(url, params=params, timeout=120)
    if not response.ok:
        raise MaapError(
            f"Catalogue search failed ({response.status_code}): {response.text[:400]}"
        )

    items = response.json().get("features", [])
    if product_type:
        items = [
            it for it in items if it["properties"].get("product:type") == product_type
        ]
    items.sort(key=lambda it: it["properties"].get("datetime", ""))
    return items


def item_summary(item: dict[str, Any]) -> dict[str, Any]:
    """Condense a STAC item to the fields worth printing or logging."""
    props = item["properties"]
    return {
        "id": item.get("id", props.get("title", "")),
        "title": props.get("title", ""),
        "product_type": props.get("product:type"),
        "datetime": props.get("datetime"),
        "start_datetime": props.get("start_datetime"),
        "end_datetime": props.get("end_datetime"),
        "mission_phase": props.get("eofeos:mission_phase"),
        "grid_code": props.get("grid:code"),
        "orbit_state": props.get("sat:orbit_state"),
        "version": props.get("version"),
        "processor": str(props.get("processing:software", "")),
        "bbox": item.get("bbox"),
    }


# --------------------------------------------------------------------------- #
# Download
# --------------------------------------------------------------------------- #


def download_asset(
    href: str,
    dest: Path,
    tokens: TokenManager,
    chunk_size: int = 1 << 20,
    overwrite: bool = False,
) -> Path:
    """Stream a token-protected asset to ``dest``.

    Downloads to a ``.part`` file first so an interrupted run never leaves a
    truncated raster that later steps would silently read.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if os.path.exists(config.long_path(dest)) and not overwrite:
        print(f"    exists, skipping: {dest.name}")
        return dest

    part = dest.with_suffix(dest.suffix + ".part")
    with requests.get(href, headers=tokens.headers(), stream=True, timeout=300) as r:
        if r.status_code in (401, 403):
            raise MaapError(
                f"Access denied for {dest.name} ({r.status_code}). The offline "
                "token is invalid, expired, or the account lacks BIOMASS access."
            )
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        written = 0
        with open(config.long_path(part), "wb") as fh:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                fh.write(chunk)
                written += len(chunk)
                if total:
                    pct = 100.0 * written / total
                    print(
                        f"\r    {dest.name}: {written/1e6:7.1f} / {total/1e6:.1f} MB"
                        f" ({pct:5.1f}%)",
                        end="",
                        flush=True,
                    )
        print()

    os.replace(config.long_path(part), config.long_path(dest))
    return dest


def download_item(
    item: dict[str, Any],
    out_dir: Path,
    tokens: TokenManager,
    assets: Iterable[str] = config.BIOMASS_ASSETS,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Download the requested assets of one STAC item into ``out_dir/<item id>``."""
    item_id = item["properties"].get("title") or item["id"]
    target = out_dir / item_id
    downloaded: dict[str, Path] = {}

    for key in assets:
        asset = item["assets"].get(key)
        if asset is None:
            print(f"    asset '{key}' not present on {item_id}")
            continue
        name = asset.get("file:local_path") or asset["href"].rsplit("/", 1)[-1]
        downloaded[key] = download_asset(
            asset["href"], target / name, tokens, overwrite=overwrite
        )

    return downloaded
