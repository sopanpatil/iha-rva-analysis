# calculate_ensemble_mean_q.py
# ------------------------------
# Computes the ensemble-mean HBV discharge (averaged across the 4 CHESS-SCAPE
# ensemble members) for each RCP, used as an alternative input to the IHA/RVA
# analysis (chess_scape_iha_analysis.py --mean), as opposed to running IHA/RVA
# on each ensemble member individually.
#
# Inputs
# ------
#   <INPUT_DIR>/<rcp>_<ens>_hbv_discharge.csv   (one per RCP/ensemble member)
#
# Outputs (in <OUTPUT_DIR>/)
# ------
#   <rcp>_hbv_ensemble_mean_discharge.csv
#
# Usage
# -----
#   python calculate_ensemble_mean_q.py
#   python calculate_ensemble_mean_q.py --rcp rcp85   # single RCP

import argparse
from pathlib import Path

import pandas as pd

INPUT_DIR = Path("./chess_scape_output")
OUTPUT_DIR = Path("./ensemble_mean_discharge")

RCPS = ["rcp26", "rcp45", "rcp60", "rcp85"]
ENSEMBLES = ["01", "04", "06", "15"]
MODEL_NAME = "hbv"


def compute_one_rcp(rcp: str) -> None:
    out_path = OUTPUT_DIR / f"{rcp}_{MODEL_NAME}_ensemble_mean_discharge.csv"
    if out_path.exists():
        print(f"  [skip] {out_path.name} already exists")
        return

    dfs = []
    for ens in ENSEMBLES:
        fname = INPUT_DIR / f"{rcp}_{ens}_{MODEL_NAME}_discharge.csv"
        if not fname.exists():
            print(f"  ERROR: missing {fname}")
            return
        df = pd.read_csv(fname, index_col=0)   # col 0 = date
        dfs.append(df)

    # Stack along a new axis and take the mean -- all DataFrames must share
    # the same index and columns, which they should if the runs are aligned.
    ensemble_mean = pd.concat(dfs, axis=0, keys=ENSEMBLES) \
                      .groupby(level=1) \
                      .mean()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ensemble_mean.to_csv(out_path)
    print(f"  Written: {out_path}  shape={ensemble_mean.shape}")


def parse_args():
    p = argparse.ArgumentParser(description="Compute ensemble-mean HBV discharge per RCP.")
    p.add_argument("--rcp", choices=RCPS, default=None, help="Single RCP. Default: all four.")
    return p.parse_args()


def main():
    args = parse_args()
    rcps = [args.rcp] if args.rcp else RCPS
    for rcp in rcps:
        print(f"RCP: {rcp}")
        compute_one_rcp(rcp)
    print("\nDone.")


if __name__ == "__main__":
    main()
