"""
benchmark_analysis.py
=====================

RVA statistics for the benchmark of the modelling chain over the baseline
period (WY 1982-2010), from IHA_Results/benchmark/iha_annual_benchmark.csv.gz
(benchmark_observed.py). Uses the same statistics as rva_member_analysis.py
(leave-one-out expected frequency, circular timing, Cliff's delta).

Comparisons (reference -> comparison; pulse thresholds of the reference):
  hydrological model : observed           -> hbv_observed_forcing
  climate forcing    : hbv_observed_forcing -> chess_scape_<member>  (each member)
  full chain         : observed           -> chess_scape_<member>  (each member)

Outputs (IHA_Results/benchmark/)
  benchmark_rva_catchment.csv.gz  one row per comparison x member x catchment x parameter
  benchmark_rva_summary.csv       per comparison and parameter: % of catchments with low
                                  alteration, departure (below/above) and concentration
                                  (mean and range over members), median Cliff's delta,
                                  and the no-change percentage with low alteration

    python benchmark_analysis.py
"""

import os

import pandas as pd

from rva_member_analysis import CIRCULAR, OUT_DIR as RVA_DIR, is_con, is_dep, rva_statistics

BENCH_DIR = os.path.join("IHA_Results", "benchmark")
MEMBERS = ["01", "04", "06", "15"]


def wide(df, param, ids):
    return df.pivot(index="catchment_id", columns="water_year", values=param).reindex(ids).astype(float)


def compare(tbl, ref, cmp_, thresholds_from, params, ids):
    base = tbl[(tbl.series == ref) & (tbl.thresholds_from == thresholds_from)]
    fut = tbl[(tbl.series == cmp_) & (tbl.thresholds_from == thresholds_from)]
    frames = []
    for p in params:
        b, f = wide(base, p, ids), wide(fut, p, ids)
        ok = b.notna().sum(axis=1) >= 10
        stats = rva_statistics(b.values[ok.values], f.values[ok.values], circular=p in CIRCULAR)
        stats.insert(0, "parameter", p)
        stats.insert(0, "catchment_id", b.index[ok.values])
        frames.append(stats)
    return pd.concat(frames, ignore_index=True)


def main():
    tbl = pd.read_csv(os.path.join(BENCH_DIR, "iha_annual_benchmark.csv.gz"))
    params = [c for c in tbl.columns if c not in ("series", "thresholds_from", "catchment_id", "water_year")]
    ids = sorted(tbl.catchment_id.unique())
    jobs = [("hydrological_model", "observed", "hbv_observed_forcing", "observed", "all")]
    for m in MEMBERS:
        jobs.append(("climate_forcing", "hbv_observed_forcing", f"chess_scape_{m}", "hbv_observed_forcing", m))
        jobs.append(("full_chain", "observed", f"chess_scape_{m}", "observed", m))
    frames = []
    for name, ref, cmp_, thr, member in jobs:
        res = compare(tbl, ref, cmp_, thr, params, ids)
        res.insert(0, "member", member)
        res.insert(0, "comparison", name)
        frames.append(res)
        print(f"  {name} ({member})")
    res = pd.concat(frames, ignore_index=True)
    res.to_csv(os.path.join(BENCH_DIR, "benchmark_rva_catchment.csv.gz"), index=False, float_format="%.6g",
                compression={"method": "gzip", "mtime": 0})

    res = res.assign(dep=is_dep(res.category), con=is_con(res.category), low=res.category.eq("low"))
    res["dep_below"] = res.dep & (res.frac_below > res.frac_above)
    res["dep_above"] = res.dep & (res.frac_above > res.frac_below)
    per = (res.groupby(["comparison", "member", "parameter"])
           [["low", "dep", "dep_below", "dep_above", "con"]].mean() * 100)
    per = per.join(res.groupby(["comparison", "member", "parameter"]).cliffs_delta.median())
    per = per.reset_index()
    summ = per.groupby(["comparison", "parameter"]).agg(
        pct_low=("low", "mean"), pct_low_min=("low", "min"), pct_low_max=("low", "max"),
        pct_departure=("dep", "mean"), pct_departure_below=("dep_below", "mean"),
        pct_departure_above=("dep_above", "mean"), pct_concentration=("con", "mean"),
        cliffs_delta=("cliffs_delta", "mean")).reset_index()
    null = pd.read_csv(os.path.join(RVA_DIR, "rva_null_catchment.csv.gz"))
    null_low = (null.category.eq("low").groupby(null.parameter).mean() * 100).rename("pct_low_nochange")
    summ = summ.merge(null_low.reset_index(), on="parameter")
    summ.round(3).to_csv(os.path.join(BENCH_DIR, "benchmark_rva_summary.csv"), index=False)
    print(f"Outputs written to {BENCH_DIR}/")


if __name__ == "__main__":
    main()
