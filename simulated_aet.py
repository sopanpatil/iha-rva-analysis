"""
simulated_aet.py
================

Baseline (WY 1982-2010) and future (WY 2051-2080) means of the actual
evapotranspiration simulated by HBV, for every RCP and ensemble member.

HBV is re-run for each catchment with the stored catchment-mean CHESS-SCAPE
forcing and the calibrated parameters, and the simulated discharge is checked
against the stored discharge, so that the evapotranspiration comes from
exactly the projection runs analysed in this study.

Inputs (from <CHESS_SCAPE_ROOT>/, as for extract_period_means.py)
  rcp<rcp>_<ens>_{pr,tas,pet}_catchment_means_combined.csv
  rcp<rcp>_<ens>_hbv_discharge.csv
  ../hbv-model/calibrated_parameters.csv

Outputs (IHA_Results/period_means/)
  period_means_aet.csv.gz       rcp, member, catchment_id, period, season (annual, DJF,
                                JJA), variable (aet), value (mm/day)
  aet_run_reproduction.csv      maximum and median absolute difference between the
                                re-run and the stored discharge, per RCP and member

Requires the hbv_model package on the Python path.

    python simulated_aet.py [--workers 8] [--limit N]
"""

import argparse
import itertools
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from benchmark_observed import PARAM_COLS
from catchment_filter import get_valid_catchments
from chess_scape_iha_analysis import CHESS_SCAPE_ROOT, KGE_PATH
from extract_period_means import MEMBERS, OUT_DIR, PERIODS, RCPS, period_means

FILES = {"pr": "pr_catchment_means_combined", "tas": "tas_catchment_means_combined",
         "pet": "pet_catchment_means_combined", "q": "hbv_discharge"}


def load(rcp, member, var, ids):
    raw = pd.read_csv(os.path.join(CHESS_SCAPE_ROOT, f"rcp{rcp}_{member}_{FILES[var]}.csv"))
    return raw[[str(c) for c in ids]].to_numpy(dtype=float)


def run_one(job):
    cid, params, pr, tas, pet = job
    from hbv_model.hbv import HBVModel
    q, states = HBVModel(params).run(pr, tas, pet)
    return cid, np.asarray(q, float), np.asarray(states["ETact"], float)


def main():
    ap = argparse.ArgumentParser(description="Period means of HBV-simulated actual evapotranspiration.")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None, help="Only the first N catchments (trial).")
    args = ap.parse_args()

    valid = get_valid_catchments(kge_path=KGE_PATH)
    ids = sorted(valid)[: args.limit] if args.limit else sorted(valid)
    pars = pd.read_csv(KGE_PATH).set_index("gauge_id")
    col = {cid: j for j, cid in enumerate(ids)}

    rows, repro = [], []
    for rcp, member in itertools.product(RCPS, MEMBERS):
        print(f"RCP{rcp} member {member}")
        pr, tas, pet, q_stored = (load(rcp, member, v, ids) for v in ("pr", "tas", "pet", "q"))
        jobs = [(cid, {p: float(pars.loc[cid, p]) for p in PARAM_COLS}, pr[:, j], tas[:, j], pet[:, j])
                for j, cid in enumerate(ids)]
        q_sim = np.full_like(pet, np.nan)
        aet = np.full_like(pet, np.nan)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for cid, q, et in pool.map(run_one, jobs, chunksize=4):
                q_sim[:, col[cid]], aet[:, col[cid]] = q, et
        diff = np.abs(q_sim - q_stored)
        repro.append({"rcp": rcp, "member": member, "n_catchments": len(ids),
                      "max_abs_diff_mm_day": float(np.nanmax(diff)),
                      "median_abs_diff_mm_day": float(np.nanmedian(diff))})
        print(f"  re-run vs stored discharge: max |diff| {np.nanmax(diff):.2e} mm/day")
        for p, bounds in PERIODS.items():
            for season, v in period_means(aet, *bounds).items():
                rows.append(pd.DataFrame({"rcp": rcp, "member": member, "catchment_id": ids,
                                          "period": p, "season": season, "variable": "aet",
                                          "value": v}))

    os.makedirs(OUT_DIR, exist_ok=True)
    pd.DataFrame(repro).to_csv(os.path.join(OUT_DIR, "aet_run_reproduction.csv"), index=False)
    out = os.path.join(OUT_DIR, "period_means_aet.csv.gz")
    pd.concat(rows, ignore_index=True).to_csv(out, index=False, float_format="%.6g",
                                              compression={"method": "gzip", "mtime": 0})
    print(f"Written {out}")


if __name__ == "__main__":
    main()
