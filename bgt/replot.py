"""Rebuilding figures from a finished run's tables, without recomputing it.

The expensive half of a comparison -- reading the 1 m ALS rasters, masking them,
warping onto the fine grid, block-reducing and searching the shift -- produces
exactly one thing the figures need: the per-cell pairs. Those are written to
``paired_cells*.csv``. So a change that only affects *drawing* (an axis range, a
colour scale, a label) does not need the arithmetic repeated; it needs the
tables read back and the plotting calls made again.

That is what ``--replot`` does, and this module holds the pieces it needs to
recover the map geometry, which the tables carry only implicitly.

Two figures cannot be rebuilt this way and are left as they are:

* ``fig06_shift_search`` draws the objective surface over thousands of candidate
  offsets, which is never written to disk -- only its optimum and verdict are;
* ``fig13_mask_map`` needs the exclusion code of *every* cell in the window, and
  ``paired_cells`` by definition holds only the cells that were kept.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def grid_from_pairs(paired: pd.DataFrame, column: str,
                    shape: tuple[int, int]) -> np.ndarray:
    """One value per cell, laid back out on the analysis window."""
    grid = np.full(shape, np.nan)
    grid[paired["row"].to_numpy(), paired["col"].to_numpy()] = \
        paired[column].to_numpy(dtype=float)
    return grid


def _axis_fit(index: np.ndarray, coord: np.ndarray) -> tuple[float, float]:
    """Cell size and origin of one axis, from cell indices and centres."""
    if np.ptp(index) == 0:
        raise ValueError("a single row or column cannot fix the cell size")
    step, centre0 = np.polyfit(index.astype(float), coord.astype(float), 1)
    return float(step), float(centre0 - 0.5 * step)


def window_from_pairs(paired: pd.DataFrame) -> tuple[tuple[int, int],
                                                     tuple[float, float, float, float]]:
    """``(shape, extent)`` of the map window, recovered from the paired cells.

    Each row carries its row and column index and the coordinate of its centre,
    which together fix the cell size and the window origin exactly. What they
    cannot give back is how far the window ran *past* the last paired cell, so a
    scene whose right or bottom margin held no paired cell is redrawn slightly
    cropped. Only blank margin is lost: the maps draw paired cells and nothing
    else.
    """
    rows = paired["row"].to_numpy()
    cols = paired["col"].to_numpy()
    cell_x, left = _axis_fit(cols, paired["x"].to_numpy())
    cell_y, top = _axis_fit(rows, paired["y"].to_numpy())
    shape = (int(rows.max()) + 1, int(cols.max()) + 1)
    right = left + shape[1] * cell_x
    bottom = top + shape[0] * cell_y          # cell_y is negative: y runs down
    return shape, (left, right, bottom, top)


def read_table(path) -> pd.DataFrame:
    """A table written by an earlier run, or an empty frame if it is absent."""
    return pd.read_csv(path) if path.is_file() else pd.DataFrame()


def exclusion_counts(path) -> dict[str, int]:
    """The per-reason cell counts of one scene, keyed as the figure expects."""
    table = read_table(path)
    if table.empty:
        return {}
    row = table.iloc[0]
    skip = {"site", "scene", "biomass_date", "cells_in_window", "scene_label"}
    return {k: int(v) for k, v in row.items() if k not in skip}
