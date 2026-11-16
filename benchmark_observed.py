"""
benchmark_observed.py
=====================

Benchmarks the modelling chain against observations over the baseline period
(WY 1982-2010), and reports low-flow performance of the calibrated HBV model.

For each retained catchment, HBV is run with the calibrated parameters and the
observed CAMELS-GB v2 forcing (the forcing used for calibration). Three outputs
are written to IHA_Results/benchmark/:

  kge_reproduction.csv
      Calibration (WY1991-2010) and validation (WY2011-2019) KGE recomputed
      from this run, next to the values stored in calibrated_parameters.csv.
      Close agreement confirms that the forcing and run set-up match the
      original calibration.

  lowflow_metrics.csv
      Per catchment and period (calibration, validation, baseline): KGE; KGE of
      log-transformed flows; percentage bias in Q95 and Q70 (flows exceeded
      95% and 70% of the time); percentage bias in low-flow volume of the flow
      duration curve (FLV, 70-100% exceedance; Yilmaz et al., 2008); and the
      median daily recession ratio Q(t)/Q(t-1) on recession days, observed and
      simulated.

  iha_annual_benchmark.csv.gz
      The 32 IHA parameters for each baseline water year, for the observed
      flows, the HBV simulation with observed forcing, and the HBV simulations
      with CHESS-SCAPE forcing (each ensemble member). Pulse parameters are
      computed with the pulse thresholds of each reference series (column
      thresholds_from: observed or hbv_observed_forcing), so that any pair of
      series can be compared with the RVA against a common threshold.

Observed and observed-forcing series use the Gregorian calendar: monthly means
use calendar months, and the day of annual minimum and maximum flow is
expressed on the 360-day scale of the CHESS-SCAPE series (day of water year
x 360 / year length). Water years with fewer than 94.4% valid daily values
(the 340-of-360 rule used for the projections) are excluded; shorter gaps in
the remaining years are filled by linear interpolation before computing IHA.

Run
---
    python benchmark_observed.py --camels-pattern "/path/to/timeseries/*_{id}_*.csv" \\
        --date-col date --precip-col precipitation --temp-col temperature \\
        --pet-col pet --q-col discharge_spec --run-start 1970-10-01

Requires the hbv_model package (hbv-model repository) on the Python path.
"""

import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from catchment_filter import get_valid_catchments
from chess_scape_iha_analysis import (
    CHESS_SCAPE_ROOT, DAYS_PER_YEAR, KGE_PATH, MIN_DAYS, N_PRE_YEARS, PRE_END,
    PRE_FIRST_WY, PRE_START, RVA_HIGH_PCTL, RVA_LOW_PCTL, _iha_one_year,
    compute_iha_period,
)

OUT_DIR = os.path.join("IHA_Results", "benchmark")
MEMBERS = ["01", "04", "06", "15"]
PARAM_COLS = ["TT", "CFMAX", "CFR", "CWH", "FC", "LP", "BETA",
              "K0", "K1", "K2", "UZL", "PERC", "MAXBAS"]
PERIODS = {
    "calibration": ("1990-10-01", "2010-09-30"),
    "validation": ("2010-10-01", "2019-09-30"),
    "baseline": ("1981-10-01", "2010-09-30"),
}
BASELINE_WY = range(PRE_FIRST_WY, PRE_FIRST_WY + N_PRE_YEARS)
MIN_FRACTION = MIN_DAYS / DAYS_PER_YEAR
WY_MONTH_ORDER = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]


# ---------------------------------------------------------------------------
# Low-flow metrics
# ---------------------------------------------------------------------------
def _kge(obs, sim):
    if np.std(obs) == 0 or np.std(sim) == 0 or np.mean(obs) == 0:
        return np.nan
    r = np.corrcoef(obs, sim)[0, 1]
    return 1 - np.sqrt((r - 1) ** 2 + (np.std(sim) / np.std(obs) - 1) ** 2
                       + (np.mean(sim) / np.mean(obs) - 1) ** 2)


def _flv(obs, sim):
    """Yilmaz et al. (2008) low-flow volume bias, 70-100% exceedance (%)."""
    eps = 1e-6
    lo = lambda q: np.log(np.maximum(np.sort(q)[: int(np.ceil(0.3 * len(q)))], eps))
    lo_o, lo_s = lo(obs), lo(sim)
    ref_o = np.sum(lo_o - lo_o.min())
    if ref_o == 0:
        return np.nan
    return -100 * (np.sum(lo_s - lo_s.min()) - ref_o) / ref_o


def _recession_ratio(q, min_length=5):
    """Median Q(t)/Q(t-1) over recession runs of >= min_length declining days,
    excluding the first two days of each run."""
    q = np.asarray(q)
    falling = np.diff(q) < 0
    ratios, run = [], []
    for i, f in enumerate(falling):
        if f and q[i] > 0:
            run.append(q[i + 1] / q[i])
        else:
            if len(run) >= min_length:
                ratios.extend(run[2:])
            run = []
    if len(run) >= min_length:
        ratios.extend(run[2:])
    return float(np.median(ratios)) if ratios else np.nan


def lowflow_metrics(obs, sim):
    ok = ~np.isnan(obs) & ~np.isnan(sim)
    obs, sim = obs[ok], sim[ok]
    if len(obs) < 365:
        return {}
    eps = 0.01 * np.mean(obs)
    pct_bias = lambda a, b: 100 * (b - a) / a if a > 0 else np.nan
    q95_o, q95_s = np.percentile(obs, 5), np.percentile(sim, 5)
    q70_o, q70_s = np.percentile(obs, 30), np.percentile(sim, 30)
    return {
        "n_days": int(len(obs)),
        "kge": _kge(obs, sim),
        "kge_log": _kge(np.log(obs + eps), np.log(sim + eps)),
        "q95_obs": q95_o, "q95_sim": q95_s, "q95_bias_pct": pct_bias(q95_o, q95_s),
        "q70_obs": q70_o, "q70_sim": q70_s, "q70_bias_pct": pct_bias(q70_o, q70_s),
        "flv_pct": _flv(obs, sim),
        "recession_ratio_obs": _recession_ratio(obs),
        "recession_ratio_sim": _recession_ratio(sim),
    }


# ---------------------------------------------------------------------------
# Gregorian-calendar IHA
# ---------------------------------------------------------------------------
def water_year(dates: pd.DatetimeIndex) -> np.ndarray:
    return dates.year + (dates.month >= 10)


def annual_iha_gregorian(q: pd.Series, low_thr: float, high_thr: float) -> pd.DataFrame:
    """IHA parameters for each baseline water year of a daily series indexed by date."""
    wy = water_year(q.index)
    records = []
    for y in BASELINE_WY:
        yr = q[wy == y]
        if len(yr) == 0 or yr.notna().mean() < MIN_FRACTION:
            records.append({})
            continue
        vals = yr.interpolate(limit_direction="both").to_numpy(dtype=float)
        row = _iha_one_year(vals, low_thr, high_thr)
        months = yr.index.month
        for pos, m in enumerate(WY_MONTH_ORDER):
            key = [k for k in row if k.startswith(f"month_{pos + 1:02d}_")][0]
            row[key] = vals[months == m].mean()
        scale = DAYS_PER_YEAR / len(vals)
        row["doy_min"] = (int(vals.argmin()) + 1) * scale
        row["doy_max"] = (int(vals.argmax()) + 1) * scale
        records.append(row)
    out = pd.DataFrame(records)
    out.insert(0, "water_year", list(BASELINE_WY))
    return out


def baseline_thresholds(q: pd.Series):
    sel = q[(q.index >= PERIODS["baseline"][0]) & (q.index <= PERIODS["baseline"][1])].dropna()
    return np.percentile(sel, RVA_LOW_PCTL), np.percentile(sel, RVA_HIGH_PCTL)


# ---------------------------------------------------------------------------
# Per-catchment work
# ---------------------------------------------------------------------------
def run_catchment(job):
    cid, params, path, cols, run_start, chess = job
    from hbv_model.hbv import HBVModel

    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df[cols["date"]])
    df = df[df["date"] >= pd.Timestamp(run_start)].set_index("date")
    forcing = df[[cols["precip"], cols["temp"], cols["pet"]]]
    if forcing.isna().any().any():
        forcing = forcing.interpolate(limit_direction="both")
    q_sim, _ = HBVModel(params).run(forcing[cols["precip"]].values,
                                   forcing[cols["temp"]].values,
                                   forcing[cols["pet"]].values)
    obs = df[cols["q"]].astype(float)
    sim = pd.Series(q_sim, index=df.index)

    kge_row = {"catchment_id": cid}
    metric_rows = []
    for name, (start, end) in PERIODS.items():
        sel = (df.index >= start) & (df.index <= end)
        m = lowflow_metrics(obs[sel].values, sim[sel].values)
        if name in ("calibration", "validation"):
            kge_row[f"{name}_kge_rerun"] = m.get("kge", np.nan)
        metric_rows.append({"catchment_id": cid, "period": name, **m})

    series = {"observed": obs, "hbv_observed_forcing": sim}
    thresholds = {k: baseline_thresholds(v) for k, v in series.items()}
    iha_frames = []
    for ref, (lo, hi) in thresholds.items():
        for name, q in series.items():
            t = annual_iha_gregorian(q, lo, hi)
            t.insert(0, "catchment_id", cid)
            t.insert(0, "thresholds_from", ref)
            t.insert(0, "series", name)
            iha_frames.append(t)
        for member, q360 in chess.items():
            t = compute_iha_period(q360, N_PRE_YEARS, lo, hi)
            t.insert(0, "water_year", list(BASELINE_WY))
            t.insert(0, "catchment_id", cid)
            t.insert(0, "thresholds_from", ref)
            t.insert(0, "series", f"chess_scape_{member}")
            iha_frames.append(t)
    return kge_row, metric_rows, pd.concat(iha_frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description="Benchmark the modelling chain against observations.")
    ap.add_argument("--camels-pattern", required=True,
                    help="Glob for one catchment's CAMELS-GB v2 time-series file, with {id} "
                         "in place of the gauge id.")
    ap.add_argument("--date-col", required=True)
    ap.add_argument("--precip-col", required=True)
    ap.add_argument("--temp-col", required=True)
    ap.add_argument("--pet-col", required=True)
    ap.add_argument("--q-col", required=True, help="Observed discharge in mm/day.")
    ap.add_argument("--run-start", required=True,
                    help="First date of the HBV run (as in the original calibration).")
    ap.add_argument("--params-csv", default=KGE_PATH)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None, help="Only the first N catchments (testing).")
    args = ap.parse_args()

    cols = {"date": args.date_col, "precip": args.precip_col, "temp": args.temp_col,
            "pet": args.pet_col, "q": args.q_col}
    valid = get_valid_catchments(kge_path=args.params_csv)
    pars = pd.read_csv(args.params_csv).set_index("gauge_id")
    ids = sorted(valid)[: args.limit] if args.limit else sorted(valid)

    print("Loading CHESS-SCAPE baseline discharge ...")
    chess = {}
    for m in MEMBERS:
        raw = pd.read_csv(f"{CHESS_SCAPE_ROOT}/rcp85_{m}_hbv_discharge.csv")
        chess[m] = raw[[str(c) for c in ids]].values[PRE_START:PRE_END + 1]

    jobs = []
    for j, cid in enumerate(ids):
        matches = glob.glob(args.camels_pattern.replace("{id}", str(cid)))
        if len(matches) != 1:
            raise FileNotFoundError(f"{len(matches)} files match catchment {cid}: "
                                    f"{args.camels_pattern.replace('{id}', str(cid))}")
        params = {p: float(pars.loc[cid, p]) for p in PARAM_COLS}
        jobs.append((cid, params, matches[0], cols, args.run_start,
                     {m: chess[m][:, j] for m in MEMBERS}))

    print(f"Running {len(jobs)} catchments ...")
    kge_rows, metric_rows, iha_frames = [], [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, (k, m, t) in enumerate(pool.map(run_catchment, jobs, chunksize=4)):
            kge_rows.append(k)
            metric_rows.extend(m)
            iha_frames.append(t)
            if (i + 1) % 100 == 0:
                print(f"  {i + 1}/{len(jobs)}")

    os.makedirs(OUT_DIR, exist_ok=True)
    kge = pd.DataFrame(kge_rows).set_index("catchment_id")
    kge["calibration_kge_file"] = pars.loc[kge.index, "calibration_kge"]
    kge["validation_kge_file"] = pars.loc[kge.index, "validation_kge"]
    kge.to_csv(os.path.join(OUT_DIR, "kge_reproduction.csv"))
    pd.DataFrame(metric_rows).to_csv(os.path.join(OUT_DIR, "lowflow_metrics.csv"),
                                     index=False, float_format="%.6g")
    pd.concat(iha_frames, ignore_index=True).to_csv(
        os.path.join(OUT_DIR, "iha_annual_benchmark.csv.gz"), index=False, float_format="%.10g",
        compression={"method": "gzip", "mtime": 0})

    d_cal = (kge.calibration_kge_rerun - kge.calibration_kge_file).abs()
    d_val = (kge.validation_kge_rerun - kge.validation_kge_file).abs()
    print("\nKGE reproduction (|rerun - stored|):")
    print(f"  calibration: median {d_cal.median():.4f}, max {d_cal.max():.4f}")
    print(f"  validation:  median {d_val.median():.4f}, max {d_val.max():.4f}")
    print(f"Outputs written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
