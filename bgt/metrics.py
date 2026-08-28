"""Validation statistics for a paired reference / product height sample.

The metric set is the one conventionally reported for canopy height product
validation: bias, MAE, RMSE and its centred (bias-free) counterpart, relative
errors, correlation, the agreement-based coefficient of determination against
the 1:1 line, and both an ordinary least-squares and a reduced major axis fit.
The last distinction matters here because *both* variables carry error, which
biases the OLS slope towards zero.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

import config


def compute_metrics(reference: np.ndarray, product: np.ndarray) -> dict[str, float]:
    """Compare ``product`` against ``reference`` (both 1-D, in metres).

    ``reference`` is the ALS aggregate, ``product`` the BIOMASS forest height.
    Residuals are defined as ``product - reference``, so a positive bias means
    BIOMASS overestimates relative to ALS.
    """
    ref = np.asarray(reference, dtype=float)
    prd = np.asarray(product, dtype=float)
    ok = np.isfinite(ref) & np.isfinite(prd)
    ref, prd = ref[ok], prd[ok]
    n = ref.size

    out: dict[str, float] = {"n": float(n)}
    if n < 3:
        return out

    resid = prd - ref
    bias = float(resid.mean())
    rmse = float(np.sqrt(np.mean(resid**2)))
    ref_mean = float(ref.mean())

    out.update(
        {
            "ref_mean": ref_mean,
            "ref_sd": float(ref.std(ddof=1)),
            "prod_mean": float(prd.mean()),
            "prod_sd": float(prd.std(ddof=1)),
            "bias": bias,
            "bias_pct": 100.0 * bias / ref_mean if ref_mean else np.nan,
            "median_bias": float(np.median(resid)),
            "mae": float(np.mean(np.abs(resid))),
            "rmse": rmse,
            "rmse_pct": 100.0 * rmse / ref_mean if ref_mean else np.nan,
            # Centred RMSE removes the constant offset: the part of the error a
            # simple recalibration could not fix.
            "rmse_centred": float(np.sqrt(max(rmse**2 - bias**2, 0.0))),
            "resid_sd": float(resid.std(ddof=1)),
            # Cells with a non-positive reference are excluded rather than
            # allowed to turn the whole statistic into nan.
            "mape_pct": float(
                np.nanmean(np.abs(resid) / np.where(ref > 0, ref, np.nan)) * 100.0
            ),
        }
    )

    out["pearson_r"] = float(np.corrcoef(ref, prd)[0, 1])
    out["r2_pearson"] = out["pearson_r"] ** 2
    out["spearman_rho"] = float(sps.spearmanr(ref, prd).statistic)

    # Nash-Sutcliffe style R^2 against the 1:1 line: how much of the reference
    # variance the product explains *as an estimator*, not after refitting.
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((ref - ref_mean) ** 2))
    out["r2_one_to_one"] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    fit = sps.linregress(ref, prd)
    out["ols_slope"] = float(fit.slope)
    out["ols_intercept"] = float(fit.intercept)
    out["ols_slope_stderr"] = float(fit.stderr)
    out["ols_p_value"] = float(fit.pvalue)

    # Reduced major axis: symmetric in the two variables, appropriate when both
    # sides carry measurement error.
    sd_ref, sd_prd = ref.std(ddof=1), prd.std(ddof=1)
    rma_slope = float(np.sign(out["pearson_r"]) * sd_prd / sd_ref) if sd_ref > 0 else np.nan
    out["rma_slope"] = rma_slope
    out["rma_intercept"] = float(prd.mean() - rma_slope * ref_mean)

    # Lin's concordance correlation: agreement, not just association.
    out["ccc"] = float(
        2 * out["pearson_r"] * sd_ref * sd_prd
        / (sd_ref**2 + sd_prd**2 + (prd.mean() - ref_mean) ** 2)
    )

    # Bland-Altman limits of agreement.
    out["loa_lower"] = bias - 1.96 * out["resid_sd"]
    out["loa_upper"] = bias + 1.96 * out["resid_sd"]

    return out


def metrics_table(
    df: pd.DataFrame,
    als_stats: list[str],
    product_col: str = "fh",
) -> pd.DataFrame:
    """One metrics row per candidate ALS statistic, best agreement first."""
    rows = []
    for stat in als_stats:
        col = f"als_{stat}"
        if col not in df.columns:
            continue
        row = {"als_stat": stat}
        row.update(compute_metrics(df[col].to_numpy(), df[product_col].to_numpy()))
        rows.append(row)
    table = pd.DataFrame(rows)
    if "rmse" in table:
        table = table.sort_values("rmse").reset_index(drop=True)
    return table


def binned_metrics(
    df: pd.DataFrame,
    ref_col: str,
    product_col: str = "fh",
    bins: list[float] | None = None,
) -> pd.DataFrame:
    """Error statistics stratified by reference height.

    Canopy height products are rarely unbiased across the whole range -- SAR and
    lidar products alike tend to overestimate short stands and saturate over
    tall ones -- so the stratified table is the one that shows the real
    behaviour.
    """
    bins = bins or config.HEIGHT_BINS
    work = df[[ref_col, product_col]].dropna()
    work = work.assign(bin=pd.cut(work[ref_col], bins=bins, right=False))

    rows = []
    for interval, group in work.groupby("bin", observed=True):
        if len(group) < 3:
            continue
        row = {
            "bin": str(interval),
            "bin_lower": interval.left,
            "bin_upper": interval.right,
            "bin_centre": (interval.left + interval.right) / 2.0,
        }
        row.update(compute_metrics(group[ref_col].to_numpy(), group[product_col].to_numpy()))
        rows.append(row)
    return pd.DataFrame(rows)


def quality_metrics(
    df: pd.DataFrame,
    ref_col: str,
    product_col: str = "fh",
    quality_col: str = "fh_quality",
    max_classes: int = 12,
) -> pd.DataFrame:
    """Error statistics stratified by the BIOMASS quality layer value.

    The semantics of the quality codes are defined in the product format
    specification; this table lets the values speak for themselves, so a
    threshold can be chosen from evidence rather than assumption.
    """
    if quality_col not in df.columns or df[quality_col].isna().all():
        return pd.DataFrame()

    work = df[[ref_col, product_col, quality_col]].dropna()
    classes = np.sort(work[quality_col].unique())
    if len(classes) > max_classes:
        # Too many distinct codes to be class labels -- bin them instead.
        work = work.assign(
            **{quality_col: pd.cut(work[quality_col], bins=max_classes).astype(str)}
        )
        classes = np.sort(work[quality_col].unique())

    rows = []
    for value in classes:
        group = work[work[quality_col] == value]
        if len(group) < 3:
            continue
        row = {"quality": value, "n_cells": len(group)}
        row.update(compute_metrics(group[ref_col].to_numpy(), group[product_col].to_numpy()))
        rows.append(row)
    return pd.DataFrame(rows)


def format_summary(m: dict[str, float], ref_label: str, prod_label: str) -> str:
    """Human-readable block for the console and the run log."""
    if m.get("n", 0) < 3:
        return f"{prod_label} vs {ref_label}: too few paired cells (n={m.get('n', 0):.0f})"
    return "\n".join(
        [
            f"{prod_label} vs {ref_label}   (residual = product - reference)",
            f"  n cells        {m['n']:.0f}",
            f"  reference      {m['ref_mean']:6.2f} +/- {m['ref_sd']:.2f} m",
            f"  product        {m['prod_mean']:6.2f} +/- {m['prod_sd']:.2f} m",
            f"  bias           {m['bias']:+6.2f} m  ({m['bias_pct']:+.1f} %)",
            f"  MAE            {m['mae']:6.2f} m",
            f"  RMSE           {m['rmse']:6.2f} m  ({m['rmse_pct']:.1f} %)",
            f"  RMSE centred   {m['rmse_centred']:6.2f} m",
            f"  Pearson r      {m['pearson_r']:6.3f}   (r2 {m['r2_pearson']:.3f})",
            f"  Spearman rho   {m['spearman_rho']:6.3f}",
            f"  R2 vs 1:1      {m['r2_one_to_one']:6.3f}",
            f"  CCC            {m['ccc']:6.3f}",
            f"  OLS fit        y = {m['ols_slope']:.3f} x {m['ols_intercept']:+.2f}",
            f"  RMA fit        y = {m['rma_slope']:.3f} x {m['rma_intercept']:+.2f}",
            f"  95 % LoA       [{m['loa_lower']:+.2f}, {m['loa_upper']:+.2f}] m",
        ]
    )


def association(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Strength and shape of the relation between two quantities.

    Unlike :func:`compute_metrics` this makes no assumption that the two sides
    share a unit, so it reports association and the fitted line only -- no bias,
    no RMSE, no 1:1 line.
    """
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]

    out: dict[str, float] = {"n": float(a.size)}
    if a.size < 3:
        return out

    fit = sps.linregress(a, b)
    out.update(
        {
            "pearson_r": float(np.corrcoef(a, b)[0, 1]),
            "spearman_rho": float(sps.spearmanr(a, b).statistic),
            "ols_slope": float(fit.slope),
            "ols_intercept": float(fit.intercept),
            "ols_slope_stderr": float(fit.stderr),
            "ols_p_value": float(fit.pvalue),
            "x_mean": float(a.mean()),
            "x_sd": float(a.std(ddof=1)),
            "y_mean": float(b.mean()),
            "y_sd": float(b.std(ddof=1)),
        }
    )
    out["r2_fit"] = out["pearson_r"] ** 2
    # Scatter about the fitted line: what the fit leaves unexplained, in the
    # unit of y.
    out["rmse_fit"] = float(np.sqrt(np.mean((b - (fit.slope * a + fit.intercept)) ** 2)))
    return out


def binned_profile(x: np.ndarray, y: np.ndarray, edges: np.ndarray) -> pd.DataFrame:
    """Median of ``y`` and its quartiles inside each bin of ``x``.

    The profile is what shows saturation: a relation can keep a respectable
    correlation while its median flattens out above some biomass.
    """
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]

    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = (a >= lo) & (a < hi)
        if not inside.any():
            continue
        values = b[inside]
        rows.append(
            {
                "x_lo": float(lo),
                "x_hi": float(hi),
                "x_centre": float((lo + hi) / 2.0),
                "n": int(inside.sum()),
                "y_q25": float(np.percentile(values, 25)),
                "y_median": float(np.median(values)),
                "y_q75": float(np.percentile(values, 75)),
            }
        )
    return pd.DataFrame(rows)
