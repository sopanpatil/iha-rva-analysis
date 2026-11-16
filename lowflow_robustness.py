"""
lowflow_robustness.py
=====================

How well HBV represents interannual low-flow behaviour, and whether the
projected low-flow alteration depends on that.

1. Interannual agreement (WY 1982-2010, HBV with observed forcing against
   observed flows, from IHA_Results/benchmark/iha_annual_benchmark.csv.gz):
   per catchment and parameter, the Spearman correlation between observed and
   simulated annual values, and the ratio of simulated to observed relative
   spread ((P75 - P25) / median).

2. Projected change in catchment subsets: departure, Cliff's delta and relative
   median change (mean over members) for the low-flow parameters, for all
   catchments, for catchments with good low-flow performance (validation KGE of
   log flows >= 0.7 and |Q95 bias| <= 50%), and for catchments in which the
   relative spread of the 7-day minimum is reproduced within a factor of 1.5.

Outputs (IHA_Results/benchmark/)
  lowflow_interannual_agreement.csv
  lowflow_subset_summary.csv

    python lowflow_robustness.py
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BENCH = os.path.join("IHA_Results", "benchmark")
RVA = os.path.join("IHA_Results", "member_rva")
PARAMS = ["min_1day", "min_7day", "min_30day", "min_90day", "base_flow_index", "month_10_Jul",
          "month_11_Aug", "month_12_Sep", "max_1day", "max_30day", "month_03_Dec"]
LOW = ["min_1day", "min_3day", "min_7day", "min_30day", "min_90day", "base_flow_index",
       "month_10_Jul", "month_11_Aug", "month_12_Sep"]


def relative_spread(x):
    med = x.median()
    return (x.quantile(0.75) - x.quantile(0.25)) / med if med > 0 else np.nan


def interannual_agreement():
    b = pd.read_csv(os.path.join(BENCH, "iha_annual_benchmark.csv.gz"))
    b = b[b.thresholds_from == "observed"]
    obs = b[b.series == "observed"].set_index(["catchment_id", "water_year"])
    sim = b[b.series == "hbv_observed_forcing"].set_index(["catchment_id", "water_year"])
    rows = []
    for p in PARAMS:
        j = pd.concat([obs[p].rename("o"), sim[p].rename("s")], axis=1).dropna()
        for cid, g in j.groupby(level=0):
            if len(g) < 10:
                continue
            rho = spearmanr(g.o, g.s)[0] if g.o.nunique() > 1 and g.s.nunique() > 1 else np.nan
            ro, rs = relative_spread(g.o), relative_spread(g.s)
            rows.append({"catchment_id": cid, "parameter": p, "n_years": len(g), "spearman_rho": rho,
                         "relative_spread_ratio": rs / ro if ro and ro > 0 else np.nan})
    return pd.DataFrame(rows)


def subset_summary(agree):
    lf = pd.read_csv(os.path.join(BENCH, "lowflow_metrics.csv"))
    val = lf[lf.period == "validation"].set_index("catchment_id")
    spread = agree[agree.parameter == "min_7day"].set_index("catchment_id").relative_spread_ratio
    subsets = {
        "all": val.index,
        "good_lowflow_performance": val.index[(val.kge_log >= 0.7) & (val.q95_bias_pct.abs() <= 50)],
        "realistic_interannual_spread": spread.index[(spread >= 1 / 1.5) & (spread <= 1.5)],
    }
    m = pd.read_csv(os.path.join(RVA, "rva_catchment_member.csv.gz"), dtype={"member": str})
    m = m[m.parameter.isin(LOW)]
    m["dep"] = m.category.isin(["mod_departure", "high_departure"])
    rows = []
    for name, ids in subsets.items():
        x = m[m.catchment_id.isin(ids)]
        g = (x.groupby(["rcp", "parameter", "member"])
             .agg(dep=("dep", "mean"), cliffs_delta=("cliffs_delta", "median"),
                  median_change=("median_change", "median"))
             .groupby(["rcp", "parameter"]).mean().reset_index())
        g["pct_departure"] = 100 * g.pop("dep")
        g.insert(0, "n_catchments", x.catchment_id.nunique())
        g.insert(0, "subset", name)
        rows.append(g)
    return pd.concat(rows, ignore_index=True)


def main():
    agree = interannual_agreement()
    agree.round(4).to_csv(os.path.join(BENCH, "lowflow_interannual_agreement.csv"), index=False)
    subset_summary(agree).round(3).to_csv(os.path.join(BENCH, "lowflow_subset_summary.csv"), index=False)
    print(f"Outputs written to {BENCH}/")


if __name__ == "__main__":
    main()
