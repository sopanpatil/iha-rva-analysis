"""
nochange_expected_frequency.py
==============================

Bias of the hydrologic alteration (HA) factor under no change, for three
choices of the expected in-range frequency f_exp:

  half        the conventional f_exp = 0.5
  in_sample   the fraction of baseline years inside their own target range
  loo         the leave-one-out coverage used in rva_member_analysis.py

For each of the 12 ordered member pairs of iha_null_member_pairs.py, the
baseline water years of the comparison member are treated as the "future"
years and compared with the target range of the reference member's
baseline. Without forced change an unbiased f_exp gives a mean HA of zero.

Inputs  (IHA_Results/annual/)
------
  iha_annual_rcp85_<member>.csv.gz
  iha_annual_null_ref<A>_cmp<B>.csv.gz

Output (IHA_Results/member_rva/)
------
  nochange_expected_frequency.csv – per parameter and f_exp choice: mean HA
                                    and percentage of comparisons classed as
                                    departure and as concentration

Run
---
    python nochange_expected_frequency.py
"""

import itertools
import os

import numpy as np
import pandas as pd

from iha_categories import classify_ha
from rva_member_analysis import (ANNUAL_DIR, CIRCULAR, MEMBERS, OUT_DIR, circular_mean_day,
                                 is_con, is_dep, leave_one_out_coverage, parameters,
                                 read_annual, recentre, target_range, wide)


def coverage(values: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Fraction of non-missing years inside [lo, hi], per catchment."""
    inside = np.sum((values >= lo[:, None]) & (values <= hi[:, None]), axis=1)
    return inside / np.sum(~np.isnan(values), axis=1)


def pair_statistics(base_tbl: pd.DataFrame, cmp_tbl: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for p in parameters(base_tbl):
        b = wide(base_tbl, p)
        c = wide(cmp_tbl, p).loc[b.index]
        base, comp = b.values, c.values
        if p in CIRCULAR:
            centre = circular_mean_day(base)
            base, comp = recentre(base, centre), recentre(comp, centre)
        lo, hi = target_range(base)
        f_obs = coverage(comp, lo, hi)
        for choice, f_exp in [("half", 0.5), ("in_sample", coverage(base, lo, hi)),
                              ("loo", leave_one_out_coverage(base))]:
            ha = (f_obs - f_exp) / f_exp
            frames.append(pd.DataFrame({"parameter": p, "f_exp_choice": choice, "ha": ha,
                                        "category": [classify_ha(h) for h in ha]}))
    return pd.concat(frames, ignore_index=True)


def main():
    baselines = {}
    for m in MEMBERS:
        tbl = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_rcp85_{m}.csv.gz"))
        baselines[m] = tbl[tbl.period == "baseline"]
    frames = []
    for ref, cmp_ in itertools.permutations(MEMBERS, 2):
        comp = read_annual(os.path.join(ANNUAL_DIR, f"iha_annual_null_ref{ref}_cmp{cmp_}.csv.gz"))
        frames.append(pair_statistics(baselines[ref], comp))
        print(f"  no-change pair {ref} vs {cmp_}")
    res = pd.concat(frames, ignore_index=True)
    res = res.assign(dep=is_dep(res.category) * 100.0, con=is_con(res.category) * 100.0)
    out = (res.groupby(["parameter", "f_exp_choice"], sort=False)
              .agg(ha_mean=("ha", "mean"), pct_departure=("dep", "mean"),
                   pct_concentration=("con", "mean"))
              .round(4).reset_index())
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "nochange_expected_frequency.csv")
    out.to_csv(path, index=False)
    print(f"Output written to {path}")


if __name__ == "__main__":
    main()
