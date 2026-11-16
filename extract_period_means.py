"""
extract_period_means.py
=======================

Baseline (WY 1982-2010) and future (WY 2051-2080) means of the catchment
forcing and HBV soil moisture for every RCP and ensemble member, for mapping
projected changes in the climate inputs and internal model states.

Variables (catchment means, from <CHESS_SCAPE_ROOT>/):
  pr   precipitation (mm/day)       rcp<rcp>_<ens>_pr_catchment_means_combined.csv
  tas  air temperature (degC)       rcp<rcp>_<ens>_tas_catchment_means_combined.csv
  pet  potential ET (mm/day)        rcp<rcp>_<ens>_pet_catchment_means_combined.csv
  q    HBV discharge (mm/day)       rcp<rcp>_<ens>_hbv_discharge.csv
  sm   HBV soil moisture (mm)       rcp<rcp>_<ens>_hbv_sm.csv

Each is averaged over the whole water year (annual), December-February (DJF)
and June-August (JJA). Actual evapotranspiration is estimated from the
long-term water balance as aet = pr - q (annual means; storage change over a
30-year period is negligible).

Output
------
  IHA_Results/period_means/period_means.csv.gz
      rcp, member, catchment_id, period, season, variable, value

Run
---
    python extract_period_means.py
"""

import itertools
import os

import numpy as np
import pandas as pd

from catchment_filter import get_valid_catchments, filter_array_columns
from chess_scape_iha_analysis import (
    CHESS_SCAPE_ROOT, KGE_PATH, N_POST_YEARS, N_PRE_YEARS, POST_END, POST_START,
    PRE_END, PRE_START, DAYS_PER_YEAR,
)

OUT_DIR = os.path.join("IHA_Results", "period_means")
MEMBERS = ["01", "04", "06", "15"]
RCPS = ["26", "45", "60", "85"]
FILES = {
    "pr": "rcp{rcp}_{ens}_pr_catchment_means_combined.csv",
    "tas": "rcp{rcp}_{ens}_tas_catchment_means_combined.csv",
    "pet": "rcp{rcp}_{ens}_pet_catchment_means_combined.csv",
    "q": "rcp{rcp}_{ens}_hbv_discharge.csv",
    "sm": "rcp{rcp}_{ens}_hbv_sm.csv",
}
# Day-of-water-year rows (0-based) in the 360-day calendar: Oct = 0-29, ..., Sep = 330-359
SEASON_ROWS = {"annual": np.arange(DAYS_PER_YEAR),
               "DJF": np.arange(60, 150),
               "JJA": np.arange(240, 330)}
PERIODS = {"baseline": (PRE_START, PRE_END, N_PRE_YEARS),
           "future": (POST_START, POST_END, N_POST_YEARS)}


def period_means(values: np.ndarray, start: int, end: int, n_years: int) -> dict:
    block = values[start:end + 1].reshape(n_years, DAYS_PER_YEAR, -1)
    return {s: np.nanmean(block[:, rows, :], axis=(0, 1)) for s, rows in SEASON_ROWS.items()}


def main():
    valid = get_valid_catchments(kge_path=KGE_PATH)
    rows = []
    for rcp, ens in itertools.product(RCPS, MEMBERS):
        print(f"RCP{rcp} member {ens}")
        means, ids = {}, None
        for var, pattern in FILES.items():
            raw = filter_array_columns(
                pd.read_csv(os.path.join(CHESS_SCAPE_ROOT, pattern.format(rcp=rcp, ens=ens))), valid)
            cols = [c for c in raw.columns if c != "date"]
            if ids is None:
                ids = cols
            elif cols != ids:
                raise ValueError(f"Catchment columns differ in {var} for RCP{rcp}/{ens}")
            vals = raw[ids].to_numpy(dtype=float)
            means[var] = {p: period_means(vals, *bounds) for p, bounds in PERIODS.items()}
        for p in PERIODS:
            means.setdefault("aet", {})[p] = {"annual": means["pr"][p]["annual"] - means["q"][p]["annual"]}
        for var, by_period in means.items():
            for p, by_season in by_period.items():
                for season, v in by_season.items():
                    rows.append(pd.DataFrame({"rcp": rcp, "member": ens, "catchment_id": ids,
                                              "period": p, "season": season,
                                              "variable": var, "value": v}))
    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, "period_means.csv.gz")
    pd.concat(rows, ignore_index=True).to_csv(out, index=False, float_format="%.6g",
                                              compression={"method": "gzip", "mtime": 0})
    print(f"Written {out}")


if __name__ == "__main__":
    main()
