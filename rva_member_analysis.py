"""
rva_member_analysis.py
======================

Range of Variability Approach (RVA) statistics computed separately for each
CHESS-SCAPE ensemble member against its own baseline, from the yearly IHA
tables written by chess_scape_iha_analysis.py.

For every RCP, ensemble member, catchment and IHA parameter:

  * Target range: 25th-75th percentile of the 29 baseline water years
    (WY 1982-2010).
  * Expected frequency f_exp: leave-one-out coverage of the baseline, i.e. the
    fraction of baseline years that fall inside the target range computed
    from the other 28 years. This is the expected in-range frequency of a new
    year drawn from an unchanged regime, and accounts for tied values in
    discrete (count) parameters.
  * HA = (f_obs - f_exp) / f_exp, with f_obs the fraction of the 30 future
    water years (WY 2051-2080) inside the target range; classified with the
    five-category scheme of iha_categories.py.
  * Direction: fractions of future years below and above the target range.
  * Magnitude: relative change in the median (%), and Cliff's delta between
    future and baseline years (P(future > baseline) - P(future < baseline)).

Day of annual minimum and maximum flow are circular quantities. They are
re-centred on each catchment's circular-mean baseline day before analysis, so
that a target range never straddles the water-year boundary; their magnitude
is reported as the shift in median day (days) rather than as a percentage.

The same statistics are computed for the member-pair baseline comparisons of
iha_null_member_pairs.py, giving the no-change reference level of each
statistic.

Inputs  (IHA_Results/annual/)
------
  iha_annual_rcp<rcp>_<member>.csv.gz
  iha_annual_null_ref<A>_cmp<B>.csv.gz

Outputs (IHA_Results/member_rva/)
-------
  rva_catchment_member.csv.gz   – one row per RCP x member x catchment x parameter
  rva_catchment_ensemble.csv.gz – one row per RCP x catchment x parameter:
                                  median HA across members, member spread
                                  (standard deviation of HA), number of members
                                  in departure and in concentration
  rva_null_catchment.csv.gz     – as rva_catchment_member, for the 12 ordered
                                  member pairs of the no-change comparison
  rva_summary_by_parameter.csv  – per RCP and parameter: percentage of
                                  catchments in each category (mean and range
                                  across members), departure below/above,
                                  median magnitude, and the no-change rates

Run
---
    python rva_member_analysis.py
"""

import itertools
import os

import numpy as np
import pandas as pd

from iha_categories import classify_ha

ANNUAL_DIR = os.path.join("IHA_Results", "annual")
OUT_DIR = os.path.join("IHA_Results", "member_rva")
MEMBERS = ["01", "04", "06", "15"]
RCPS = ["26", "45", "60", "85"]
CIRCULAR = {"doy_min", "doy_max"}
DAYS_PER_YEAR = 360
LOW_PCTL, HIGH_PCTL = 25, 75


# ---------------------------------------------------------------------------
# Core statistics (arrays are catchments x years)
# ---------------------------------------------------------------------------
def circular_mean_day(days: np.ndarray) -> np.ndarray:
    angle = 2 * np.pi * (days - 0.5) / DAYS_PER_YEAR
    mean = np.arctan2(np.nanmean(np.sin(angle), axis=1), np.nanmean(np.cos(angle), axis=1))
    return mean * DAYS_PER_YEAR / (2 * np.pi) + 0.5


def recentre(days: np.ndarray, centre: np.ndarray) -> np.ndarray:
    """Days relative to centre, wrapped to [-180, 180)."""
    half = DAYS_PER_YEAR / 2
    return ((days - centre[:, None] + half) % DAYS_PER_YEAR) - half


def target_range(base: np.ndarray):
    return (np.nanpercentile(base, LOW_PCTL, axis=1),
            np.nanpercentile(base, HIGH_PCTL, axis=1))


def leave_one_out_coverage(base: np.ndarray) -> np.ndarray:
    """Fraction of baseline years inside the target range of the other years."""
    n_years = base.shape[1]
    inside = np.zeros(base.shape[0])
    counted = np.zeros(base.shape[0])
    for i in range(n_years):
        rest = np.delete(base, i, axis=1)
        lo, hi = target_range(rest)
        value = base[:, i]
        ok = ~np.isnan(value)
        inside += ok & (value >= lo) & (value <= hi)
        counted += ok
    return inside / counted


def cliffs_delta(base: np.ndarray, fut: np.ndarray) -> np.ndarray:
    """P(future > baseline) - P(future < baseline), per catchment."""
    diff = fut[:, :, None] - base[:, None, :]
    valid = ~np.isnan(diff)
    gt = np.nansum(diff > 0, axis=(1, 2))
    lt = np.nansum(diff < 0, axis=(1, 2))
    return (gt - lt) / valid.sum(axis=(1, 2))


def rva_statistics(base: np.ndarray, fut: np.ndarray, circular: bool) -> pd.DataFrame:
    """RVA and magnitude statistics for one parameter (catchments x years)."""
    if circular:
        centre = circular_mean_day(base)
        base, fut = recentre(base, centre), recentre(fut, centre)

    lo, hi = target_range(base)
    n_fut = np.sum(~np.isnan(fut), axis=1)
    f_obs = np.sum((fut >= lo[:, None]) & (fut <= hi[:, None]), axis=1) / n_fut
    below = np.sum(fut < lo[:, None], axis=1) / n_fut
    above = np.sum(fut > hi[:, None], axis=1) / n_fut
    f_exp = leave_one_out_coverage(base)
    ha = (f_obs - f_exp) / f_exp

    med_base = np.nanmedian(base, axis=1)
    med_fut = np.nanmedian(fut, axis=1)
    if circular:
        change = med_fut - med_base                       # days
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            change = np.where(med_base != 0, 100 * (med_fut - med_base) / np.abs(med_base), np.nan)

    return pd.DataFrame({
        "f_exp": f_exp, "f_obs": f_obs, "ha": ha,
        "category": [classify_ha(h) for h in ha],
        "frac_below": below, "frac_above": above,
        "median_change": change, "cliffs_delta": cliffs_delta(base, fut),
        "target_low": lo, "target_high": hi,
    })


# ---------------------------------------------------------------------------
# Table handling
# ---------------------------------------------------------------------------
def read_annual(path: str) -> pd.DataFrame:
    return pd.read_csv(path).set_index(["catchment_id", "water_year"])


def parameters(table: pd.DataFrame) -> list:
    return [c for c in table.columns if c != "period"]


def wide(table: pd.DataFrame, param: str) -> pd.DataFrame:
    """catchments x water years for one parameter."""
    return table[param].unstack("water_year").astype(float)


def compare_tables(base_tbl: pd.DataFrame, fut_tbl: pd.DataFrame, params: list) -> pd.DataFrame:
    frames = []
    for p in params:
        b, f = wide(base_tbl, p), wide(fut_tbl, p)
        f = f.loc[b.index]
        stats = rva_statistics(b.values, f.values, circular=p in CIRCULAR)
        stats.insert(0, "parameter", p)
        stats.insert(0, "catchment_id", b.index.values)
        frames.append(stats)
    return pd.concat(frames, ignore_index=True)


def member_results(params=None) -> pd.DataFrame:
    frames = []
    for rcp, member in itertools.product(RCPS, MEMBERS):
        tbl = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_rcp{rcp}_{member}.csv.gz"))
        params = params or parameters(tbl)
        res = compare_tables(tbl[tbl.period == "baseline"], tbl[tbl.period == "future"], params)
        res.insert(0, "member", member)
        res.insert(0, "rcp", rcp)
        frames.append(res)
        print(f"  RCP{rcp} member {member}")
    return pd.concat(frames, ignore_index=True)


def null_results() -> pd.DataFrame:
    frames = []
    baselines = {}
    for m in MEMBERS:
        tbl = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_rcp85_{m}.csv.gz"))
        baselines[m] = tbl[tbl.period == "baseline"]
    params = parameters(baselines[MEMBERS[0]])
    for ref, cmp_ in itertools.permutations(MEMBERS, 2):
        comp = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_null_ref{ref}_cmp{cmp_}.csv.gz"))
        res = compare_tables(baselines[ref], comp, params)
        res.insert(0, "comparison_member", cmp_)
        res.insert(0, "reference_member", ref)
        frames.append(res)
        print(f"  no-change pair {ref} vs {cmp_}")
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
def is_dep(cat):
    return cat.isin(["mod_departure", "high_departure"])


def is_con(cat):
    return cat.isin(["mod_concentration", "high_concentration"])


def ensemble_summary(res: pd.DataFrame) -> pd.DataFrame:
    res = res.assign(dep=is_dep(res.category), con=is_con(res.category))
    g = res.groupby(["rcp", "catchment_id", "parameter"])
    out = g.agg(ha_median=("ha", "median"), ha_sd=("ha", "std"),
                n_departure=("dep", "sum"), n_concentration=("con", "sum"),
                frac_below_median=("frac_below", "median"),
                frac_above_median=("frac_above", "median"),
                median_change_median=("median_change", "median"),
                cliffs_delta_median=("cliffs_delta", "median")).reset_index()
    out["category_median_ha"] = [classify_ha(h) for h in out.ha_median]
    return out


def parameter_summary(res: pd.DataFrame, null: pd.DataFrame) -> pd.DataFrame:
    res = res.assign(dep=is_dep(res.category), con=is_con(res.category),
                     high_dep=res.category.eq("high_departure"))
    res["dep_below"] = res.dep & (res.frac_below > res.frac_above)
    res["dep_above"] = res.dep & (res.frac_above > res.frac_below)
    per_member = (res.groupby(["rcp", "parameter", "member"])
                  [["dep", "high_dep", "con", "dep_below", "dep_above"]].mean() * 100)
    mags = res.groupby(["rcp", "parameter", "member"])[["median_change", "cliffs_delta"]].median()
    per_member = per_member.join(mags).reset_index()

    agg = per_member.groupby(["rcp", "parameter"]).agg(
        pct_departure=("dep", "mean"), pct_departure_min=("dep", "min"),
        pct_departure_max=("dep", "max"), pct_high_departure=("high_dep", "mean"),
        pct_departure_below=("dep_below", "mean"), pct_departure_above=("dep_above", "mean"),
        pct_concentration=("con", "mean"), pct_concentration_min=("con", "min"),
        pct_concentration_max=("con", "max"),
        median_change=("median_change", "mean"), cliffs_delta=("cliffs_delta", "mean"),
    ).reset_index()

    null = null.assign(dep=is_dep(null.category), con=is_con(null.category))
    nul = (null.groupby("parameter")[["dep", "con"]].mean() * 100).rename(
        columns={"dep": "pct_departure_nochange", "con": "pct_concentration_nochange"})
    nul["ha_mean_nochange"] = null.groupby("parameter").ha.mean()
    return agg.merge(nul.reset_index(), on="parameter").round(3)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Member results ...")
    res = member_results()
    print("No-change results ...")
    null = null_results()

    res.to_csv(os.path.join(OUT_DIR, "rva_catchment_member.csv.gz"), index=False, float_format="%.6g", compression={"method": "gzip", "mtime": 0})
    null.to_csv(os.path.join(OUT_DIR, "rva_null_catchment.csv.gz"), index=False, float_format="%.6g", compression={"method": "gzip", "mtime": 0})
    ensemble_summary(res).to_csv(os.path.join(OUT_DIR, "rva_catchment_ensemble.csv.gz"),
                                 index=False, float_format="%.6g", compression={"method": "gzip", "mtime": 0})
    summary = parameter_summary(res, null)
    summary.to_csv(os.path.join(OUT_DIR, "rva_summary_by_parameter.csv"), index=False)
    print(f"Outputs written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
