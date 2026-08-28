"""BIOMASS Level-1b HV backscatter, put on the AGBD grid.

The L1b DGM product is a detected, ground-range image: a four-band float32
GeoTIFF (HH, HV, VH, VV) holding beta-nought *amplitude*, with no geotransform,
because its rows are azimuth times and its columns ground range. Everything
needed to place it on the ground travels beside it:

* the annotation XML gives the azimuth time step and the polynomial that turns
  ground range into slant-range time;
* the LUT netCDF gives latitude, longitude and the beta-to-gamma (or -sigma)
  conversion factor on a regular grid in (azimuth time, slant-range time).

So a pixel is geolocated by turning its row and column into those two times and
interpolating the LUT there. Only the window covering the AGBD map is read --
a few hundred pixels a side, pulled out of the 128 MB scene by range request --
and the pixels landing in each AGBD cell are averaged in linear power.
"""

from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
import requests
from pyproj import Transformer
from rasterio.crs import CRS
from rasterio.windows import Window
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import shape

import config
from bgt import maap

#: How often a dropped auxiliary download is retried, and the base pause between
#: attempts, in seconds. The pause grows with the attempt number.
_DOWNLOAD_ATTEMPTS = 4
_RETRY_PAUSE_S = 5.0

#: Extra pixels kept around the window the AGBD footprint projects to, so that
#: the interpolation and the rounding to cell indices stay inside the read.
_WINDOW_MARGIN_PX = 12


class L1bError(RuntimeError):
    """Raised when an L1b product cannot be located or read."""


# --------------------------------------------------------------------------- #
# Catalogue
# --------------------------------------------------------------------------- #


def parse_datetime(value: str) -> datetime:
    """Parse a STAC timestamp into an aware UTC datetime."""
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc)


def search_scenes(bbox: tuple[float, float, float, float],
                  collection: str = config.L1B_COLLECTION,
                  datetime_range: str | None = None) -> list[dict[str, Any]]:
    """Every L1b item intersecting ``bbox``, oldest first."""
    items = maap.search_items(bbox=bbox, collection=collection,
                              product_type=None, datetime_range=datetime_range)
    return sorted(items, key=lambda it: it["properties"]["datetime"])


def covering_scenes(items: list[dict[str, Any]], footprint,
                    min_coverage: float = 0.999) -> list[dict[str, Any]]:
    """Items whose own footprint contains (almost all of) ``footprint``.

    A scene that clips the AGBD map would compare only part of it, so partial
    coverage is dropped rather than quietly reducing the sample.
    """
    kept = []
    for item in items:
        geom = shape(item["geometry"])
        if geom.intersection(footprint).area >= min_coverage * footprint.area:
            kept.append(item)
    return kept


def nearest_in_time(items: list[dict[str, Any]], when: datetime) -> dict[str, Any]:
    """The item acquired closest in time to ``when``."""
    if not items:
        raise L1bError("No L1b scene covers the AGBD footprint")
    return min(items,
               key=lambda it: abs(parse_datetime(it["properties"]["datetime"]) - when))


def scene_summary(item: dict[str, Any]) -> dict[str, Any]:
    """The handful of item fields worth printing or writing to a table."""
    props = item["properties"]
    return {
        "id": item["id"],
        "product_type": props.get("product:type"),
        "datetime": props.get("datetime"),
        "orbit_state": props.get("sat:orbit_state"),
        "grid_code": props.get("grid:code"),
        "mission_phase": props.get("eofeos:mission_phase"),
        "polarisations": ",".join(props.get("sar:polarizations", [])),
    }


def asset_href(item: dict[str, Any], key: str) -> str:
    """Download URL of one of an item's assets."""
    asset = item.get("assets", {}).get(key)
    if not asset or not asset.get("href"):
        raise L1bError(f"Item {item['id']} has no {key!r} asset")
    return asset["href"]


def fetch_annotation(item: dict[str, Any], tokens: maap.TokenManager,
                     cache_dir: Path = config.L1B_DIR) -> tuple[Path, Path]:
    """Download this scene's annotation XML and LUT netCDF, caching both.

    Together they are about 45 MB; the measurement GeoTIFF beside them is 128 MB
    and is never downloaded whole -- see :func:`read_window`.
    """
    scene_dir = cache_dir / item["id"]
    scene_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for key in ("enclosure_annot_xml", "enclosure_nc"):
        href = asset_href(item, key)
        path = scene_dir / href.rsplit("/", 1)[-1]
        if not path.is_file():
            print(f"  fetching {path.name}")
            _download_with_retry(href, path, tokens)
        paths.append(path)
    return paths[0], paths[1]


def _download_with_retry(href: str, path: Path, tokens: maap.TokenManager,
                         attempts: int = _DOWNLOAD_ATTEMPTS) -> Path:
    """Download, retrying a dropped connection.

    The catalogue truncates a transfer often enough that a single failure is not
    worth losing a run over. ``download_asset`` writes through a ``.part`` file,
    so a truncated attempt leaves nothing a later step could mistake for a
    complete download; each retry simply starts over.
    """
    for attempt in range(1, attempts + 1):
        try:
            return maap.download_asset(href, path, tokens)
        except (requests.RequestException, OSError) as error:
            if attempt == attempts:
                raise L1bError(
                    f"{path.name} could not be downloaded in {attempts} "
                    f"attempts: {error}"
                ) from error
            print(f"    transfer failed ({error}); retrying "
                  f"{attempt + 1}/{attempts}")
            time.sleep(_RETRY_PAUSE_S * attempt)
    raise AssertionError("unreachable")


# --------------------------------------------------------------------------- #
# Image geometry
# --------------------------------------------------------------------------- #


def _text(root: ET.Element, path: str) -> str:
    node = root.find(path)
    if node is None or node.text is None:
        raise L1bError(f"Annotation has no {path}")
    return node.text.strip()


@dataclass(frozen=True)
class ImageGeometry:
    """How an L1b image row and column map to azimuth and slant-range time."""

    n_lines: int
    n_samples: int
    azimuth_time_interval_s: float
    range_pixel_spacing_m: float
    ground_to_slant: np.ndarray   #: polynomial coefficients, lowest order first
    polarisations: tuple[str, ...]
    start_time: str

    def slant_range_time(self, sample: np.ndarray) -> np.ndarray:
        """Slant-range time, in seconds, of ground-range column ``sample``."""
        ground_range = np.asarray(sample, dtype=float) * self.range_pixel_spacing_m
        return np.polyval(self.ground_to_slant[::-1], ground_range)

    def azimuth_time(self, line: np.ndarray) -> np.ndarray:
        """Time, in seconds since the first line, of image row ``line``."""
        return np.asarray(line, dtype=float) * self.azimuth_time_interval_s

    def sample_of_slant_range_time(self, times: np.ndarray) -> np.ndarray:
        """Inverse of :meth:`slant_range_time`, by interpolation.

        The forward polynomial is monotonic across the swath, so inverting it on
        the image's own sample axis is exact at the samples and smooth between
        them -- and it avoids the separate slant-to-ground polynomial, whose
        coefficients are stated in a different argument.
        """
        samples = np.arange(self.n_samples, dtype=float)
        return np.interp(times, self.slant_range_time(samples), samples,
                         left=np.nan, right=np.nan)


def read_geometry(annotation_path: Path) -> ImageGeometry:
    """Parse the image geometry out of the L1b annotation XML."""
    root = ET.parse(annotation_path).getroot()
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]  # drop the namespace

    image = root.find("sarImage")
    if image is None:
        raise L1bError(f"{annotation_path.name} has no sarImage block")
    conversion = image.find("rangeCoordinateConversion/coordinateConversion")
    if conversion is None:
        raise L1bError(f"{annotation_path.name} has no range coordinate conversion")

    pols = tuple(
        node.text.strip() for node
        in root.findall("acquisitionInformation/polarisationList/polarisation")
        if node.text
    )
    return ImageGeometry(
        n_lines=int(_text(image, "numberOfLines")),
        n_samples=int(_text(image, "numberOfSamples")),
        azimuth_time_interval_s=float(_text(image, "azimuthTimeInterval")),
        range_pixel_spacing_m=float(_text(image, "rangePixelSpacing")),
        ground_to_slant=np.array(
            [float(v) for v in _text(conversion, "groundToSlantCoefficients").split()]
        ),
        polarisations=pols or config.L1B_POLARISATIONS,
        start_time=_text(image, "firstLineAzimuthTime"),
    )


# --------------------------------------------------------------------------- #
# Look-up tables
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Luts:
    """Geolocation and radiometry on the LUT (azimuth, slant-range) grid."""

    azimuth_time_s: np.ndarray
    slant_range_time_s: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    calibration: np.ndarray   #: beta-nought to gamma- or sigma-nought factor

    def interpolator(self, values: np.ndarray) -> RegularGridInterpolator:
        return RegularGridInterpolator(
            (self.azimuth_time_s, self.slant_range_time_s), values,
            method="linear", bounds_error=False, fill_value=np.nan,
        )


def read_luts(lut_path: Path, radiometry: str = config.L1B_RADIOMETRY) -> Luts:
    """Read the geolocation and radiometric-conversion LUTs.

    The netCDF is HDF5 underneath, and h5py reads its two coordinate axes as
    well as the grids; the GDAL netCDF driver exposes only the grids.
    """
    import h5py

    with h5py.File(lut_path, "r") as handle:
        def grid(name: str) -> np.ndarray:
            data = np.asarray(handle[name][...], dtype=float)
            fill = handle[name].attrs.get("_FillValue")
            if fill is not None:
                data[data == float(np.ravel(fill)[0])] = np.nan
            return data

        return Luts(
            azimuth_time_s=np.asarray(handle["relativeAzimuthTimeRGC"][...],
                                      dtype=float),
            slant_range_time_s=np.asarray(handle["slantRangeTimeRGC"][...],
                                          dtype=float),
            latitude=grid("geometry/latitude"),
            longitude=grid("geometry/longitude"),
            calibration=grid(f"radiometry/{radiometry}"),
        )


# --------------------------------------------------------------------------- #
# Reading the window that covers a target grid
# --------------------------------------------------------------------------- #


def _window_for_bbox(geometry: ImageGeometry, luts: Luts,
                     bbox: tuple[float, float, float, float]) -> Window:
    """The image window whose pixels can fall inside the lon/lat ``bbox``.

    The LUT is coarse enough to scan whole, so the window is found from the LUT
    nodes inside the box rather than by inverting the geolocation.
    """
    lon_min, lat_min, lon_max, lat_max = bbox
    inside = (
        (luts.latitude >= lat_min) & (luts.latitude <= lat_max)
        & (luts.longitude >= lon_min) & (luts.longitude <= lon_max)
    )
    rows, cols = np.where(inside)
    if rows.size == 0:
        raise L1bError("The scene geolocation grid does not reach the AGBD footprint")

    lines = luts.azimuth_time_s[rows] / geometry.azimuth_time_interval_s
    samples = geometry.sample_of_slant_range_time(luts.slant_range_time_s[cols])
    samples = samples[np.isfinite(samples)]
    if samples.size == 0:
        raise L1bError("The AGBD footprint falls outside the scene range swath")

    row_off = max(int(np.floor(lines.min())) - _WINDOW_MARGIN_PX, 0)
    col_off = max(int(np.floor(samples.min())) - _WINDOW_MARGIN_PX, 0)
    row_end = min(int(np.ceil(lines.max())) + _WINDOW_MARGIN_PX, geometry.n_lines)
    col_end = min(int(np.ceil(samples.max())) + _WINDOW_MARGIN_PX, geometry.n_samples)
    if row_end <= row_off or col_end <= col_off:
        raise L1bError("The AGBD footprint and the scene do not overlap")
    return Window(col_off, row_off, col_end - col_off, row_end - row_off)


def _band_index(geometry: ImageGeometry, tags: dict[str, str],
                polarisation: str) -> int:
    """1-based band holding ``polarisation``, from the file's own channel order."""
    sequence = tags.get("PolarisationsSequence", "").split()
    order = sequence or list(geometry.polarisations)
    if polarisation not in order:
        raise L1bError(f"{polarisation} is not among the scene channels {order}")
    return order.index(polarisation) + 1


@dataclass(frozen=True)
class BackscatterWindow:
    """A geolocated patch of one polarimetric channel, in linear power."""

    power: np.ndarray        #: gamma- or sigma-nought, linear
    latitude: np.ndarray
    longitude: np.ndarray
    window: Window
    radiometry: str
    polarisation: str


def read_window(href: str, geometry: ImageGeometry, luts: Luts,
                bbox: tuple[float, float, float, float],
                tokens: maap.TokenManager,
                polarisation: str = config.L1B_POLARISATION,
                radiometry: str = config.L1B_RADIOMETRY) -> BackscatterWindow:
    """Read and calibrate the part of ``href`` that covers ``bbox``.

    A remote measurement file is read with an HTTP range request, so only the
    tiles intersecting the window travel. Pixels are beta-nought amplitude:
    squaring gives beta-nought power, and the LUT factor takes that to
    gamma-nought (or sigma-nought).
    """
    window = _window_for_bbox(geometry, luts, bbox)

    with rasterio.Env(
        GDAL_HTTP_HEADERS=f"Authorization: Bearer {tokens.access_token}",
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tiff,.tif",
    ):
        with rasterio.open(_vsi_path(href)) as src:
            band = _band_index(geometry, src.tags(), polarisation)
            amplitude = src.read(band, window=window, masked=True)
    amplitude = amplitude.filled(np.nan).astype(float)

    rows = window.row_off + np.arange(window.height, dtype=float)
    cols = window.col_off + np.arange(window.width, dtype=float)
    points = np.stack(
        np.meshgrid(geometry.azimuth_time(rows), geometry.slant_range_time(cols),
                    indexing="ij"),
        axis=-1,
    )

    return BackscatterWindow(
        power=amplitude**2 * luts.interpolator(luts.calibration)(points),
        latitude=luts.interpolator(luts.latitude)(points),
        longitude=luts.interpolator(luts.longitude)(points),
        window=window,
        radiometry=radiometry,
        polarisation=polarisation,
    )


def _vsi_path(href: str) -> str:
    """A local path stays a path; a URL is read through the GDAL HTTP driver."""
    return href if "://" not in href else "/vsicurl/" + href


def to_grid(patch: BackscatterWindow, transform, crs: CRS,
            shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Average the patch onto a target raster grid.

    Returns the mean linear power per cell (nan where no pixel landed) and the
    number of contributing pixels. The average is taken in linear power and only
    the result turned into dB: averaging dB would bias the mean low.
    """
    finite = np.isfinite(patch.power) & np.isfinite(patch.latitude)
    to_target = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = to_target.transform(patch.longitude[finite], patch.latitude[finite])

    col, row = ~transform * (x, y)
    col = np.floor(col).astype(int)
    row = np.floor(row).astype(int)

    n_rows, n_cols = shape
    keep = (row >= 0) & (row < n_rows) & (col >= 0) & (col < n_cols)
    flat = row[keep] * n_cols + col[keep]

    total = np.zeros(n_rows * n_cols)
    count = np.zeros(n_rows * n_cols)
    np.add.at(total, flat, patch.power[finite][keep])
    np.add.at(count, flat, 1.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / count, np.nan)
    return mean.reshape(shape), count.reshape(shape)


def to_db(power: np.ndarray) -> np.ndarray:
    """Linear power to dB, dropping non-positive values instead of warning."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(power > 0, 10.0 * np.log10(power), np.nan)
