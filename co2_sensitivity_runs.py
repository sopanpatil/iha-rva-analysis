"""
co2_sensitivity_runs.py
=======================

Sensitivity of the projected flow-regime change to the stomatal response of
vegetation to rising CO2, for RCP8.5 and the four ensemble members.

The CHESS-SCAPE PET used in the projections is FAO-56 Penman-Monteith
reference evapotranspiration (Allen et al., 1998) with a fixed surface
resistance rs = 70 s/m, for which rs/ra = 0.34 u2 (u2: wind speed at 2 m).
Following Yang et al. (2019), surface resistance increases linearly with the
atmospheric CO2 concentration C (ppm) by 0.09% per ppm:

    rs(C) = 70 [1 + 0.0009 (C - 300)] / [1 + 0.0009 (C_ref - 300)]

normalised so that rs = 70 s/m at C_ref, the mean concentration of the
calibration period (WY 1991-2010), to which the calibrated HBV parameters are
adapted. Because rs enters the Penman-Monteith equation only through the
denominator, the CO2-adjusted PET follows from the original PET exactly as

    PET_CO2 = PET [D + g (1 + 0.34 u2)] / [D + g (1 + 0.34 u2 rs(C)/70)]

with D the slope of the saturation vapour pressure curve at the daily mean
temperature and g the psychrometric constant at the catchment's mean
elevation. u2 is taken, in order of preference, from the CHESS-SCAPE
catchment-mean daily wind speed, from the catchment's mean observed wind speed
over 1991-2010 in CAMELS-GB v2 (--camels-pattern), or as the FAO-56 default of
2 m/s; wind speeds at 10 m are converted to 2 m with the FAO-56 log profile.
The log records which source was used.

CO2 concentrations: RCP8.5 (Meinshausen et al., 2011), read from the RCP
database file if available, otherwise linearly interpolated from the decadal
values embedded below.

For each member, HBV is run twice from the start of the CHESS-SCAPE record
with the calibrated parameters: with the original PET (control, which must
reproduce the stored discharge) and with the CO2-adjusted PET. The yearly IHA
parameters of the CO2 run are written in the same format as
IHA_Results/annual/, with pulse thresholds from the run's own baseline.

Outputs (IHA_Results/co2_sensitivity/)
  iha_annual_rcp85_<member>_co2pet.csv.gz   yearly IHA parameters, CO2-adjusted run
  control_reproduction.csv                   control vs stored discharge, per member
  period_means_co2pet.csv.gz                 baseline/future means of pr, pet, q, aet
                                             for the original and CO2-adjusted runs
  co2_concentrations_used.csv                year, CO2 (ppm), rs/70

Run
---
    python co2_sensitivity_runs.py --topo-csv /path/to/camels_gb_v2_topographic_attributes.csv
        [--limit N] [--workers 8] [--co2-file RCP85_MIDYEAR_CONCENTRATIONS.DAT]
        [--camels-pattern PATTERN] [--u2 U2] [--out-dir DIR] [--write-control]

No wind data were used for the published results: they are for the FAO-56 default
u2 = 2 m/s (--u2 2; IHA_Results/co2_sensitivity/) and for u2 = 4 m/s, the upper end
for Great Britain (--u2 4 --out-dir IHA_Results/co2_sensitivity_u2_4). The two bracket
the effect, since the CO2 effect on PET increases with wind speed.

With --write-control, the yearly IHA parameters of the control run are also
written (iha_annual_rcp85_<member>_control.csv.gz).

Requires the hbv_model package (hbv-model repository) on the Python path.
"""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from benchmark_observed import PARAM_COLS
from catchment_filter import get_valid_catchments
from chess_scape_iha_analysis import (
    CHESS_SCAPE_ROOT, KGE_PATH, N_POST_YEARS, N_PRE_YEARS, POST_END, POST_START, PRE_END,
    PRE_START, RVA_HIGH_PCTL, RVA_LOW_PCTL, annual_table, compute_iha_period,
)
from extract_period_means import SEASON_ROWS, period_means

OUT_DIR = os.path.join("IHA_Results", "co2_sensitivity")
MEMBERS = ["01", "04", "06", "15"]
S_R = 0.0009                 # fractional increase in rs per ppm CO2 (Yang et al., 2019)
CAL_YEARS = range(1991, 2011)  # calendar years spanning the calibration period
# RCP8.5 CO2 (ppm), Meinshausen et al. (2011); used if the database file is not available
CO2_DECADAL = {1980: 338.4, 1990: 353.9, 2000: 368.9, 2010: 389.3, 2020: 415.8, 2030: 448.8,
               2040: 489.4, 2050: 540.5, 2060: 603.5, 2070: 677.1, 2080: 758.2}


def co2_series(co2_file):
    """Annual RCP8.5 CO2 (ppm) for 1980-2080, and its source."""
    years = np.arange(1980, 2081)
    if co2_file and os.path.exists(co2_file):
        rows = []
        for line in open(co2_file):
            parts = line.split()
            if len(parts) > 4 and parts[0].isdigit() and 1700 < int(parts[0]) < 2600:
                rows.append((int(parts[0]), float(parts[3])))
        s = pd.Series(dict(rows))
        ref = pd.Series(CO2_DECADAL)
        dev = (s.reindex(ref.index) - ref).abs().max()
        if not dev < 3:
            raise ValueError(f"CO2 column in {co2_file} differs from the expected values by {dev:.1f} ppm; "
                             "check the column order.")
        return s.reindex(years), f"file {co2_file}"
    s = pd.Series(CO2_DECADAL).reindex(years).interpolate()
    return s, "decadal values (Meinshausen et al., 2011), linearly interpolated"


def svp_slope(t):
    """Slope of the saturation vapour pressure curve (kPa/degC), FAO-56 Eq. 13."""
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    return 4098 * es / (t + 237.3) ** 2


def psychrometric_constant(z):
    """kPa/degC at elevation z (m), FAO-56 Eqs. 7-8."""
    p = 101.3 * ((293 - 0.0065 * z) / 293) ** 5.26
    return 0.000665 * p


def adjust_pet(pet, tas, u2, gamma, rs_ratio):
    d = svp_slope(tas)
    return pet * (d + gamma * (1 + 0.34 * u2)) / (d + gamma * (1 + 0.34 * u2 * rs_ratio))


def camels_wind(pattern, ids):
    """Mean observed wind speed at 10 m (m/s) over 1991-2010 per catchment, or None."""
    import glob
    out = []
    for cid in ids:
        files = glob.glob(pattern.replace("{id}", str(cid)))
        if len(files) != 1:
            return None, f"{len(files)} CAMELS files for catchment {cid}"
        df = pd.read_csv(files[0])
        cols = [c for c in df.columns if "wind" in c.lower()]
        if len(cols) != 1:
            return None, f"wind columns in CAMELS file: {cols}"
        sel = df[df.columns[0]].astype(str).str[:4].astype(int).between(1991, 2010)
        out.append(df.loc[sel, cols[0]].mean())
    return np.array(out, float), f"CAMELS-GB v2 column '{cols[0]}', mean 1991-2010"


def run_one(job):
    cid, params, pr, tas, pet, pet_co2 = job
    from hbv_model.hbv import HBVModel
    q_ctrl, _ = HBVModel(params).run(pr, tas, pet)
    q_co2, _ = HBVModel(params).run(pr, tas, pet_co2)
    return cid, np.asarray(q_ctrl, float), np.asarray(q_co2, float)


def load(var, member, ids):
    name = {"pr": "pr_catchment_means_combined", "tas": "tas_catchment_means_combined",
            "pet": "pet_catchment_means_combined", "q": "hbv_discharge",
            "sfcWind": "sfcWind_catchment_means_combined"}[var]
    path = os.path.join(CHESS_SCAPE_ROOT, f"rcp85_{member}_{name}.csv")
    if not os.path.exists(path):
        return None, None
    raw = pd.read_csv(path)
    return raw[[str(c) for c in ids]].to_numpy(dtype=float), raw["date"].astype(str).to_numpy()


def main():
    ap = argparse.ArgumentParser(description="CO2 stomatal-response sensitivity of the RCP8.5 projections.")
    ap.add_argument("--limit", type=int, default=None, help="Only the first N catchments (trial).")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--co2-file", default=None, help="RCP85_MIDYEAR_CONCENTRATIONS.DAT, if downloaded.")
    ap.add_argument("--topo-csv", required=True,
                    help="CAMELS-GB v2 topographic attributes file (catchment mean elevation).")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--u2", type=float, default=None,
                    help="Use this constant wind speed at 2 m (m/s) instead of looking up wind data.")
    ap.add_argument("--write-control", action="store_true",
                    help="Also write the yearly IHA parameters of the control run.")
    ap.add_argument("--camels-pattern", default=None,
                    help="Glob for one catchment's CAMELS-GB v2 time series, {id} = gauge id (for wind).")
    args = ap.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    valid = get_valid_catchments(kge_path=KGE_PATH)
    ids = sorted(valid)[: args.limit] if args.limit else sorted(valid)
    pars = pd.read_csv(KGE_PATH).set_index("gauge_id")
    elev = pd.read_csv(args.topo_csv).set_index("gauge_id").loc[ids, "elev_mean"].to_numpy(float)
    gamma = psychrometric_constant(elev)

    avail = sorted({f.split("_", 2)[2].rsplit(".", 1)[0] for f in os.listdir(CHESS_SCAPE_ROOT)
                    if f.startswith("rcp85_01_")})
    print("CHESS-SCAPE files available for rcp85_01:", ", ".join(avail))

    co2, source = co2_series(args.co2_file)
    c_ref = co2.loc[list(CAL_YEARS)].mean()
    rs_ratio = (1 + S_R * (co2 - 300)) / (1 + S_R * (c_ref - 300))
    pd.DataFrame({"year": co2.index, "co2_ppm": co2.values, "rs_over_70": rs_ratio.values}).to_csv(
        os.path.join(out_dir, "co2_concentrations_used.csv"), index=False, float_format="%.5g")
    print(f"CO2 source: {source}; reference {c_ref:.1f} ppm (mean of {CAL_YEARS[0]}-{CAL_YEARS[-1]})")
    print(f"rs/70 in 2051-2080: {rs_ratio.loc[2051:2080].mean():.3f}")

    obs_wind, wind_note = camels_wind(args.camels_pattern, ids) if args.camels_pattern else (None, "not requested")
    print(f"Observed wind: {wind_note}")

    control_rows, mean_rows = [], []
    for member in MEMBERS:
        print(f"\nMember {member}")
        pr, dates = load("pr", member, ids)
        tas, _ = load("tas", member, ids)
        pet, _ = load("pet", member, ids)
        q_stored, _ = load("q", member, ids)
        wind, _ = load("sfcWind", member, ids)
        if not -10 < np.nanmean(tas) < 40:
            raise ValueError(f"Mean tas {np.nanmean(tas):.1f}: expected degC.")
        if args.u2 is not None:
            u2 = np.full_like(pet, args.u2)
            print(f"  u2 = {args.u2} m/s (set on the command line)")
        elif wind is not None:
            u2 = wind * 4.87 / np.log(67.8 * 10 - 5.42)
            print("  u2 from CHESS-SCAPE sfcWind (10 m, converted to 2 m)")
        elif obs_wind is not None:
            u2 = np.broadcast_to(obs_wind * 4.87 / np.log(67.8 * 10 - 5.42), pet.shape)
            print(f"  u2 from observed catchment-mean wind (median {np.median(u2[0]):.2f} m/s at 2 m)")
        else:
            u2 = np.full_like(pet, 2.0)
            print("  u2 = 2 m/s (FAO-56 default; no wind data found)")
        year_ratio = rs_ratio.reindex([int(d[:4]) for d in dates]).to_numpy()[:, None]
        pet_co2 = adjust_pet(pet, tas, u2, gamma[None, :], year_ratio)

        jobs = [(cid, {p: float(pars.loc[cid, p]) for p in PARAM_COLS},
                 pr[:, j], tas[:, j], pet[:, j], pet_co2[:, j]) for j, cid in enumerate(ids)]
        q_ctrl = np.full_like(pet, np.nan)
        q_co2 = np.full_like(pet, np.nan)
        col = {cid: j for j, cid in enumerate(ids)}
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for k, (cid, qc, qx) in enumerate(pool.map(run_one, jobs, chunksize=4)):
                q_ctrl[:, col[cid]], q_co2[:, col[cid]] = qc, qx
                if (k + 1) % 100 == 0:
                    print(f"  {k + 1}/{len(jobs)}")

        diff = np.abs(q_ctrl - q_stored)
        control_rows.append({"member": member, "n_catchments": len(ids),
                             "max_abs_diff_mm_day": float(np.nanmax(diff)),
                             "median_abs_diff_mm_day": float(np.nanmedian(diff))})
        print(f"  control vs stored discharge: max |diff| {np.nanmax(diff):.2e} mm/day")

        runs = {"co2pet": q_co2, **({"control": q_ctrl} if args.write_control else {})}
        for run, q_run in runs.items():
            frames = []
            for j, cid in enumerate(ids):
                pre, post = q_run[PRE_START:PRE_END + 1, j], q_run[POST_START:POST_END + 1, j]
                lo, hi = np.percentile(pre, RVA_LOW_PCTL), np.percentile(pre, RVA_HIGH_PCTL)
                frames.append(annual_table(cid, compute_iha_period(pre, N_PRE_YEARS, lo, hi),
                                           compute_iha_period(post, N_POST_YEARS, lo, hi)))
            pd.concat(frames, ignore_index=True).to_csv(
                os.path.join(out_dir, f"iha_annual_rcp85_{member}_{run}.csv.gz"), index=False,
                float_format="%.10g", compression={"method": "gzip", "mtime": 0})

        periods = {"baseline": (PRE_START, PRE_END, N_PRE_YEARS),
                   "future": (POST_START, POST_END, N_POST_YEARS)}
        for run, (pet_run, q_run) in {"original": (pet, q_ctrl), "co2pet": (pet_co2, q_co2)}.items():
            for p, bounds in periods.items():
                m_pr, m_pet, m_q = (period_means(v, *bounds) for v in (pr, pet_run, q_run))
                for season in SEASON_ROWS:
                    for var, vals in (("pr", m_pr), ("pet", m_pet), ("q", m_q)):
                        mean_rows.append(pd.DataFrame({"member": member, "run": run, "period": p,
                                                       "season": season, "variable": var,
                                                       "catchment_id": ids, "value": vals[season]}))
                mean_rows.append(pd.DataFrame({"member": member, "run": run, "period": p,
                                               "season": "annual", "variable": "aet",
                                               "catchment_id": ids,
                                               "value": m_pr["annual"] - m_q["annual"]}))

    pd.DataFrame(control_rows).to_csv(os.path.join(out_dir, "control_reproduction.csv"), index=False)
    pd.concat(mean_rows, ignore_index=True).to_csv(
        os.path.join(out_dir, "period_means_co2pet.csv.gz"), index=False, float_format="%.6g",
        compression={"method": "gzip", "mtime": 0})
    print(f"\nOutputs written to {out_dir}/")


if __name__ == "__main__":
    main()
