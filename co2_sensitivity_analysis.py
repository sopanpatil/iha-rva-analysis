"""
co2_sensitivity_analysis.py
===========================

Compares the projected RCP8.5 alteration with the original PET and with the
CO2-adjusted PET (co2_sensitivity_runs.py), for each IHA parameter, using
the member-wise RVA statistics of rva_member_analysis.py.

Outputs (written to the directory of the CO2 runs, default IHA_Results/co2_sensitivity/)
  co2_summary_by_parameter.csv   per parameter and run (original, and the CO2-adjusted run):
                                 percentage of catchments in departure (mean of members, with
                                 range), high departure, departure below the range, median
                                 change and Cliff's delta (mean of the member medians)
  co2_water_balance.csv          median across catchments of the change (%) in annual and
                                 JJA PET, annual AET and annual and JJA discharge, per run

    python co2_sensitivity_analysis.py
    python co2_sensitivity_analysis.py --co2-dir IHA_Results/co2_sensitivity_u2_4 --run-name co2pet_u2_4
"""

import argparse
import os

import pandas as pd

from rva_member_analysis import (ANNUAL_DIR, MEMBERS, compare_tables, is_dep, parameters,
                                 read_annual)

CO2_DIR = os.path.join("IHA_Results", "co2_sensitivity")


def results(other_dir=CO2_DIR, other_file="co2pet", other_name="co2pet"):
    frames = []
    for m in MEMBERS:
        for run, path in {"original": os.path.join(ANNUAL_DIR, f"iha_annual_rcp85_{m}.csv.gz"),
                          other_name: os.path.join(other_dir, f"iha_annual_rcp85_{m}_{other_file}.csv.gz")}.items():
            tbl = read_annual(path)
            res = compare_tables(tbl[tbl.period == "baseline"], tbl[tbl.period == "future"], parameters(tbl))
            res.insert(0, "member", m)
            res.insert(0, "run", run)
            frames.append(res)
    return pd.concat(frames, ignore_index=True)


def summary(res):
    res = res.assign(dep=is_dep(res.category), high_dep=res.category.eq("high_departure"))
    res["dep_below"] = res.dep & (res.frac_below > res.frac_above)
    per = (res.groupby(["run", "parameter", "member"])[["dep", "high_dep", "dep_below"]].mean() * 100)
    per = per.join(res.groupby(["run", "parameter", "member"])[["median_change", "cliffs_delta"]].median())
    per = per.reset_index()
    return per.groupby(["run", "parameter"]).agg(
        pct_departure=("dep", "mean"), pct_departure_min=("dep", "min"), pct_departure_max=("dep", "max"),
        pct_high_departure=("high_dep", "mean"), pct_departure_below=("dep_below", "mean"),
        median_change=("median_change", "mean"), cliffs_delta=("cliffs_delta", "mean"),
    ).round(3).reset_index()


def water_balance(co2_dir=CO2_DIR):
    m = pd.read_csv(os.path.join(co2_dir, "period_means_co2pet.csv.gz"))
    w = m.pivot_table(index=["member", "run", "season", "variable", "catchment_id"],
                      columns="period", values="value").reset_index()
    w["change_pct"] = 100 * (w.future - w.baseline) / w.baseline
    med = w.groupby(["run", "season", "variable", "member"]).change_pct.median()
    return med.groupby(["run", "season", "variable"]).mean().round(2).unstack("run").reset_index()


def main():
    ap = argparse.ArgumentParser(description="Summarise the CO2 sensitivity runs.")
    ap.add_argument("--co2-dir", default=CO2_DIR, help="Output directory of co2_sensitivity_runs.py.")
    ap.add_argument("--run-name", default="co2pet", help="Label of the CO2-adjusted run in the summary.")
    args = ap.parse_args()
    summary(results(args.co2_dir, "co2pet", args.run_name)).to_csv(
        os.path.join(args.co2_dir, "co2_summary_by_parameter.csv"), index=False)
    water_balance(args.co2_dir).to_csv(os.path.join(args.co2_dir, "co2_water_balance.csv"), index=False)
    print(f"Outputs written to {args.co2_dir}/")


if __name__ == "__main__":
    main()
