"""
chess_scape_iha_analysis.py
============================

IHA analysis on CHESS-SCAPE HBV-simulated discharge for all 671
CAMELS-GB v2 catchments, for a given RCP and ensemble member (passed
as command-line arguments -- see Run section below).

Calendar
--------
CHESS-SCAPE uses a 360-day calendar (12 × 30-day months).  Standard
pandas DatetimeIndex cannot represent this (Feb 30, Nov 30, etc. do
not exist), so dates are handled as plain integer indices within each
water year rather than as datetime objects.  A bespoke 360-day IHA
class is used instead of the generic iha_pkg.IHA.

Water year definition
---------------------
Oct 1 of year Y-1  ->  Sep 30 of year Y
Month labels in Group 1 follow water-year order:
  month_01 = October, month_02 = November, ..., month_12 = September

Row indices are located by direct date-string lookup rather than
arithmetic, so they are robust to any upstream padding in the file.

Periods
-------
Pre-alteration  : WY 1982-2010  (Oct 1 1981 - Sep 30 2010, 29 years)
Post-alteration : WY 2051-2080  (Oct 1 2050 - Sep 30 2080, 30 years)

Inputs
------
  ../hbv-model/calibrated_parameters.csv        (KGE filter, via catchment_filter.py)
  <CHESS_SCAPE_ROOT>/<rcp>_<ens>_hbv_discharge.csv

Output
------
Three wide-format CSVs (one row per catchment), written to IHA_Results/:

  rva_results_rcp<rcp>_<ens>.csv    – hydrologic alteration factor + category
                                      for each of the 32 IHA parameters
  iha_pre_rcp<rcp>_<ens>.csv        – per-parameter mean across the pre period
  iha_post_rcp<rcp>_<ens>.csv       – per-parameter mean across the post period
  rva_summary_by_parameter_rcp<rcp>_<ens>.csv

Run
---
    python chess_scape_iha_analysis.py --rcp 85 --ensemble 15
    python chess_scape_iha_analysis.py --rcp 26 --ensemble 01

Requires numpy, pandas, and catchment_filter.py (from this repository).
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from catchment_filter import get_valid_catchments, filter_array_columns
from iha_categories import classify_ha, CAT_LABELS_SIMPLE

CHESS_SCAPE_ROOT = "./chess_scape_output"
OUT_DIR = "IHA_Results"
KGE_PATH = "../hbv-model/calibrated_parameters.csv"

# Row indices (0-based, inclusive) verified by direct date lookup in the raw file
PRE_START,  PRE_END  = 300,   10739   # Oct 1 1981 – Sep 30 2010  (WY 1982–2010)
POST_START, POST_END = 25140, 35939   # Oct 1 2050 – Sep 30 2080  (WY 2051–2080)
N_PRE_YEARS  = 29    # baseline is 29 water years
N_POST_YEARS = 30    # future is 30 water years

DAYS_PER_YEAR = 360          # 360-day calendar
MIN_DAYS      = 340          # minimum days for a year to be included

RVA_LOW_PCTL  = 25
RVA_HIGH_PCTL = 75
EXPECTED_FREQ = 0.5   # RVA null expectation -- matches the 50% of the
                      # baseline distribution enclosed by the 25th/75th
                      # percentile target range (Equation 1, Methods 3.4).


# ---------------------------------------------------------------------------
# 360-day IHA engine
# ---------------------------------------------------------------------------
# Month number for each day-of-year (1-indexed) in the 360-day calendar.
# Day 1-30 = month 1, day 31-60 = month 2, …, day 331-360 = month 12.
_DOY_TO_MONTH = np.repeat(np.arange(1, 13), 30)   # length 360


def _run_lengths(mask: np.ndarray) -> list:
    """Return lengths of consecutive True runs in a boolean 1-D array."""
    durations = []
    count = 0
    for v in mask:
        if v:
            count += 1
        else:
            if count:
                durations.append(count)
                count = 0
    if count:
        durations.append(count)
    return durations


def _iha_one_year(values: np.ndarray,
                  low_thr: float,
                  high_thr: float) -> dict:
    """
    Compute all 32 IHA parameters for a single water year.

    Parameters
    ----------
    values   : 1-D array of length 360 (daily discharge, mm/day)
    low_thr  : low-pulse threshold (computed from full period)
    high_thr : high-pulse threshold (computed from full period)

    Returns
    -------
    dict of {parameter_name: value}
    """
    row = {}

    # --- Group 1: monthly means in water-year order ------------------------
    # Oct=month_01, Nov=month_02, ..., Sep=month_12
    # In the 360-day block each month is 30 contiguous rows.
    WY_MONTHS = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
                 "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
    for wy_pos, label in enumerate(WY_MONTHS):
        s = wy_pos * 30
        row[f"month_{wy_pos+1:02d}_{label}"] = values[s:s + 30].mean()

    # --- Group 2: N-day rolling minima and maxima --------------------------
    for w in (1, 3, 7, 30, 90):
        if w == 1:
            rolled = values
        else:
            # manual rolling mean to avoid pandas overhead in a tight loop
            cs = np.cumsum(np.concatenate(([0.0], values)))
            rolled = (cs[w:] - cs[:-w]) / w
        row[f"min_{w}day"] = rolled.min()
        row[f"max_{w}day"] = rolled.max()

    # base flow index: 7-day min / annual mean
    ann_mean = values.mean()
    row["base_flow_index"] = (
        row["min_7day"] / ann_mean if ann_mean > 0 else np.nan
    )

    # --- Group 3: timing (day-of-360-year, 1-indexed) ----------------------
    row["doy_min"] = int(values.argmin()) + 1
    row["doy_max"] = int(values.argmax()) + 1

    # --- Group 4: pulses ---------------------------------------------------
    low_runs  = _run_lengths(values < low_thr)
    high_runs = _run_lengths(values > high_thr)

    row["low_pulse_count"]           = len(low_runs)
    row["low_pulse_duration_mean"]   = float(np.mean(low_runs))  if low_runs  else 0.0
    row["high_pulse_count"]          = len(high_runs)
    row["high_pulse_duration_mean"]  = float(np.mean(high_runs)) if high_runs else 0.0

    # --- Group 5: rate and frequency of change -----------------------------
    diffs = np.diff(values)
    rises = diffs[diffs > 0]
    falls = diffs[diffs < 0]

    signs = np.sign(diffs)
    signs = signs[signs != 0]
    reversals = int(np.sum(np.diff(signs) != 0)) if len(signs) > 1 else 0

    row["mean_rise_rate"]      = float(rises.mean()) if len(rises) > 0 else 0.0
    row["mean_fall_rate"]      = float(-falls.mean()) if len(falls) > 0 else 0.0
    row["number_of_reversals"] = reversals

    return row


def compute_iha_period(period_array: np.ndarray, n_years: int,
                        low_thr: float, high_thr: float) -> pd.DataFrame:
    """
    Compute IHA parameters for one period (pre or post) for a single catchment.

    Parameters
    ----------
    period_array : 2-D array of shape (10800, ) — one catchment's daily
                   discharge for 30 complete water years (10800 = 30 × 360)
    low_thr, high_thr : low/high pulse thresholds (mm/day). These must be
                   computed once from the BASELINE period only (standard IHA
                   convention, Richter et al. 1996; manuscript Table 1) and
                   passed in unchanged for both the pre and post calls, so
                   that pulse counts in both periods are referenced to the
                   same fixed flow thresholds.

    Returns
    -------
    DataFrame of shape (30, n_params) — one row per water year.
    """
    assert len(period_array) == n_years * DAYS_PER_YEAR, \
        f"Expected {n_years * DAYS_PER_YEAR} rows, got {len(period_array)}"

    records = []
    for y in range(n_years):
        yr_vals = period_array[y * DAYS_PER_YEAR : (y + 1) * DAYS_PER_YEAR]
        if np.isnan(yr_vals).sum() > (DAYS_PER_YEAR - MIN_DAYS):
            records.append({})   # too many NaN — skip year
            continue
        records.append(_iha_one_year(yr_vals, low_thr, high_thr))

    return pd.DataFrame(records)


def compute_rva_catchment(pre_df: pd.DataFrame,
                           post_df: pd.DataFrame) -> dict:
    """
    Apply the RVA to compare pre and post IHA DataFrames for one catchment.

    Returns a flat dict: {param_ha: value, param_cat: label, …} for all params.
    """
    result = {}
    for col in pre_df.columns:
        pre_vals  = pre_df[col].dropna().values
        post_vals = post_df[col].dropna().values

        if len(pre_vals) < 2 or len(post_vals) < 1:
            result[f"{col}_ha"]  = np.nan
            result[f"{col}_cat"] = "insufficient data"
            continue

        rva_low  = np.percentile(pre_vals, RVA_LOW_PCTL)
        rva_high = np.percentile(pre_vals, RVA_HIGH_PCTL)

        n_post   = len(post_vals)
        n_within = int(((post_vals >= rva_low) & (post_vals <= rva_high)).sum())
        obs_freq = n_within / n_post

        ha = (obs_freq - EXPECTED_FREQ) / EXPECTED_FREQ

        # Five-category directional scheme (manuscript Section 3.4), which
        # distinguishes departure (negative HA, future flows shifting BELOW
        # the historical interquartile range) from concentration (positive
        # HA, future flows clustering more tightly WITHIN it) as physically
        # distinct signals, rather than collapsing both into a single
        # "alteration magnitude" as the standard 3-category RVA scheme does.
        # Boundary convention: every category boundary (0.33, 0.67, -0.33,
        # -0.67) belongs to the lower-magnitude-of-alteration category, so
        # bins are contiguous with no gaps or overlaps.
        cat = CAT_LABELS_SIMPLE[classify_ha(ha)]

        result[f"{col}_ha"]  = round(ha, 4)
        result[f"{col}_cat"] = cat

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="IHA/RVA analysis on CHESS-SCAPE HBV discharge for one RCP, "
                    "either for a single ensemble member or the ensemble mean."
    )
    parser.add_argument("--rcp", required=True, choices=["26", "45", "60", "85"],
                        help="RCP scenario, e.g. 85 for RCP8.5.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ensemble", choices=["01", "04", "06", "15"],
                        help="Single ensemble member, e.g. 01 or 15.")
    group.add_argument("--mean", action="store_true",
                        help="Use the ensemble-mean discharge instead of a single member "
                            "(requires calculate_ensemble_mean_q.py to have been run first).")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    rcp = args.rcp

    if args.mean:
        label = "ENSEMBLE"
        data_path = f"./ensemble_mean_discharge/rcp{rcp}_hbv_ensemble_mean_discharge.csv"
    else:
        label = args.ensemble
        data_path = f"{CHESS_SCAPE_ROOT}/rcp{rcp}_{args.ensemble}_hbv_discharge.csv"

    os.makedirs(OUT_DIR, exist_ok=True)

    # -----------------------------------------------------------------------
    # 1. Load data — read all rows, keep date for verification only
    # -----------------------------------------------------------------------
    print("Loading data …")
    raw = pd.read_csv(data_path)

    catchment_ids = [c for c in raw.columns if c != "date"]
    n_catchments  = len(catchment_ids)
    print(f"  {n_catchments} catchments, {len(raw)} total rows")

    valid_ids     = get_valid_catchments(kge_path=KGE_PATH)
    raw           = filter_array_columns(raw, valid_ids)
    catchment_ids = [c for c in raw.columns if c != "date"]
    n_catchments  = len(catchment_ids)
    print(f"  {n_catchments} catchments retained after KGE filter")

    # Verify period boundaries
    pre_dates  = raw["date"].iloc[[PRE_START,  PRE_END]]
    post_dates = raw["date"].iloc[[POST_START, POST_END]]
    print(f"  Pre  period rows {PRE_START}–{PRE_END}:   "
          f"{pre_dates.iloc[0][:10]} -> {pre_dates.iloc[1][:10]}  (WY 1982-2010, Oct-Sep)")
    print(f"  Post period rows {POST_START}–{POST_END}: "
          f"{post_dates.iloc[0][:10]} -> {post_dates.iloc[1][:10]}  (WY 2051-2080, Oct-Sep)")

    # Extract the two period blocks as numpy arrays: shape (10800, 671)
    discharge = raw[catchment_ids].values   # (36000, 671)
    pre_block  = discharge[PRE_START  : PRE_END  + 1]   # (10800, 671)
    post_block = discharge[POST_START : POST_END + 1]   # (10800, 671)

    del discharge, raw   # free memory
    print(f"  Pre block shape:  {pre_block.shape}")
    print(f"  Post block shape: {post_block.shape}")

    # -----------------------------------------------------------------------
    # 2. Loop over catchments
    # -----------------------------------------------------------------------
    rva_rows  = {}
    pre_means  = {}
    post_means = {}

    print(f"\nProcessing {n_catchments} catchments …")
    for i, cid in enumerate(catchment_ids):
        if (i + 1) % 100 == 0 or i == 0:
            print(f"  [{i+1:3d}/{n_catchments}] catchment {cid}")

        # Pulse thresholds fixed from the baseline period only (standard IHA
        # convention; manuscript Table 1), then reused unchanged for the
        # post period so pulse counts in both periods are referenced to the
        # same fixed flow thresholds.
        low_thr  = np.percentile(pre_block[:, i], RVA_LOW_PCTL)
        high_thr = np.percentile(pre_block[:, i], RVA_HIGH_PCTL)

        pre_iha  = compute_iha_period(pre_block[:, i],  N_PRE_YEARS,  low_thr, high_thr)
        post_iha = compute_iha_period(post_block[:, i], N_POST_YEARS, low_thr, high_thr)

        rva_rows[cid]   = compute_rva_catchment(pre_iha, post_iha)
        pre_means[cid]  = pre_iha.mean().round(4).to_dict()
        post_means[cid] = post_iha.mean().round(4).to_dict()

    # -----------------------------------------------------------------------
    # 3. Assemble and save output tables
    # -----------------------------------------------------------------------
    print("\nWriting outputs …")

    rva_df  = pd.DataFrame.from_dict(rva_rows,  orient="index")
    pre_df  = pd.DataFrame.from_dict(pre_means,  orient="index")
    post_df = pd.DataFrame.from_dict(post_means, orient="index")

    rva_df.index.name  = "catchment_id"
    pre_df.index.name  = "catchment_id"
    post_df.index.name = "catchment_id"

    rva_path  = os.path.join(OUT_DIR, f"rva_results_rcp{rcp}_{label}.csv")
    pre_path  = os.path.join(OUT_DIR, f"iha_pre_rcp{rcp}_{label}.csv")
    post_path = os.path.join(OUT_DIR, f"iha_post_rcp{rcp}_{label}.csv")

    rva_df.to_csv(rva_path)
    pre_df.to_csv(pre_path)
    post_df.to_csv(post_path)

    print(f"  {rva_path}")
    print(f"  {pre_path}")
    print(f"  {post_path}")

    # -----------------------------------------------------------------------
    # 4. Quick summary: departure/concentration fraction per parameter
    #    ("total departure fraction" = moderate + high departure combined,
    #    matching the terminology used to rank parameters in the manuscript)
    # -----------------------------------------------------------------------
    ha_cols = [c for c in rva_df.columns if c.endswith("_ha")]
    cat_cols = [c for c in rva_df.columns if c.endswith("_cat")]

    param_names = [c[:-3] for c in ha_cols]   # strip "_ha"

    print("\n--- Departure / concentration fraction per parameter ---")
    summary_rows = []
    for ha_col, param in zip(ha_cols, param_names):
        cat_col = f"{param}_cat"
        cats = rva_df[cat_col]
        n_high_dep  = (cats == CAT_LABELS_SIMPLE["high_departure"]).sum()
        n_mod_dep   = (cats == CAT_LABELS_SIMPLE["mod_departure"]).sum()
        n_high_con  = (cats == CAT_LABELS_SIMPLE["high_concentration"]).sum()
        n_mod_con   = (cats == CAT_LABELS_SIMPLE["mod_concentration"]).sum()
        mean_ha = rva_df[ha_col].mean()
        summary_rows.append({
            "parameter":            param,
            "mean_HA":              round(mean_ha, 3),
            "frac_departure":       round((n_high_dep + n_mod_dep) / n_catchments, 3),
            "frac_high_departure":  round(n_high_dep / n_catchments, 3),
            "frac_concentration":   round((n_high_con + n_mod_con) / n_catchments, 3),
            "frac_high_concentration": round(n_high_con / n_catchments, 3),
        })

    summary_df = (pd.DataFrame(summary_rows)
                    .sort_values("frac_departure", ascending=False))
    print(summary_df.to_string(index=False))

    summary_path = os.path.join(OUT_DIR, f"rva_summary_by_parameter_rcp{rcp}_{label}.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  {summary_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
