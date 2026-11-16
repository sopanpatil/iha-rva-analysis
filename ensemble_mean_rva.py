"""
ensemble_mean_rva.py
====================

RVA departure and concentration computed from ensemble-mean discharge, with
the same statistics as rva_member_analysis.py (leave-one-out expected
frequency, circular timing), for comparison with the member-wise analysis.
Because both use the same method, the difference between them is the effect
of averaging discharge across members before the analysis.

Inputs  (IHA_Results/annual/)
------
  iha_annual_rcp<rcp>_ENSEMBLE.csv.gz   (chess_scape_iha_analysis.py --mean)

Output (IHA_Results/member_rva/)
------
  rva_summary_ensemble_mean.csv – per RCP and parameter: percentage of
                                  catchments in departure and in concentration

Run
---
    python ensemble_mean_rva.py
"""

import os

import pandas as pd

from rva_member_analysis import (ANNUAL_DIR, OUT_DIR, RCPS, compare_tables, is_con, is_dep,
                                 parameters, read_annual)


def main():
    frames = []
    for rcp in RCPS:
        tbl = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_rcp{rcp}_ENSEMBLE.csv.gz"))
        res = compare_tables(tbl[tbl.period == "baseline"], tbl[tbl.period == "future"], parameters(tbl))
        res = res.assign(dep=is_dep(res.category) * 100.0, con=is_con(res.category) * 100.0)
        summ = (res.groupby("parameter", sort=False)
                   .agg(pct_departure=("dep", "mean"), pct_concentration=("con", "mean"))
                   .round(3).reset_index())
        summ.insert(0, "rcp", rcp)
        frames.append(summ)
        print(f"  RCP{rcp}")
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "rva_summary_ensemble_mean.csv")
    pd.concat(frames, ignore_index=True).to_csv(path, index=False)
    print(f"Output written to {path}")


if __name__ == "__main__":
    main()
