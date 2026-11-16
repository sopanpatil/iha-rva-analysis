"""
iha_null_member_pairs.py
========================

Baseline-period IHA tables for building an empirical no-change distribution
of the hydrologic alteration (HA) factor.

For every ordered pair of CHESS-SCAPE ensemble members (reference A,
comparison B), the baseline-period (WY 1982-2010) discharge of member B is
summarised with the 32 IHA parameters, using the low and high pulse
thresholds of member A's baseline. Comparing these values against member A's
own baseline IHA values (written by chess_scape_iha_analysis.py) gives HA
values for two samples of the historical climate, i.e. HA in the absence of
forced change.

The four RCP scenarios share the same forcing before December 2010
(Robinson et al., 2023), so the baseline period is identical across RCPs.
This script verifies that for every member and records the result.

Inputs
------
  ../hbv-model/calibrated_parameters.csv        (KGE filter, via catchment_filter.py)
  <CHESS_SCAPE_ROOT>/rcp<rcp>_<ens>_hbv_discharge.csv   (all RCPs and members)

Outputs (IHA_Results/annual/)
-----------------------------
  iha_annual_null_ref<A>_cmp<B>.csv.gz  – catchment_id, water_year and the 32
                                          IHA parameters of member B's baseline,
                                          computed with member A's pulse thresholds
  baseline_identity_check.csv           – maximum absolute difference in
                                          baseline daily discharge between each
                                          RCP and the reference RCP, per member

Run
---
    python iha_null_member_pairs.py            # reference RCP8.5 files
    python iha_null_member_pairs.py --rcp 26

Requires numpy, pandas, and chess_scape_iha_analysis.py / catchment_filter.py
(from this repository).
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd

from catchment_filter import get_valid_catchments, filter_array_columns
from chess_scape_iha_analysis import (
    ANNUAL_DIR, CHESS_SCAPE_ROOT, KGE_PATH, N_PRE_YEARS, PRE_END, PRE_FIRST_WY,
    PRE_START, RVA_HIGH_PCTL, RVA_LOW_PCTL, compute_iha_period,
)

MEMBERS = ["01", "04", "06", "15"]
RCPS = ["26", "45", "60", "85"]


def load_baseline(rcp: str, member: str, valid_ids: set):
    """Return (catchment_ids, baseline daily discharge array) for one run."""
    path = f"{CHESS_SCAPE_ROOT}/rcp{rcp}_{member}_hbv_discharge.csv"
    raw = filter_array_columns(pd.read_csv(path), valid_ids)
    ids = [c for c in raw.columns if c != "date"]
    return ids, raw[ids].values[PRE_START:PRE_END + 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    parser.add_argument("--rcp", choices=RCPS, default="85",
                        help="RCP whose files supply the baseline discharge (default 85).")
    args = parser.parse_args()

    valid_ids = get_valid_catchments(kge_path=KGE_PATH)
    os.makedirs(ANNUAL_DIR, exist_ok=True)

    print("Loading baseline discharge ...")
    baseline, ids = {}, None
    for member in MEMBERS:
        member_ids, block = load_baseline(args.rcp, member, valid_ids)
        if ids is None:
            ids = member_ids
        elif member_ids != ids:
            raise ValueError(f"Catchment columns differ for member {member}")
        baseline[member] = block

    print("Checking that baseline discharge is identical across RCPs ...")
    checks = []
    for member in MEMBERS:
        for rcp in RCPS:
            if rcp == args.rcp:
                continue
            other_ids, block = load_baseline(rcp, member, valid_ids)
            same_ids = other_ids == ids
            max_diff = float(np.nanmax(np.abs(block - baseline[member]))) if same_ids else np.nan
            checks.append({"member": member, "rcp": rcp, "reference_rcp": args.rcp,
                           "same_catchments": same_ids, "max_abs_diff_mm_day": max_diff})
            print(f"  member {member}, RCP{rcp} vs RCP{args.rcp}: max |diff| = {max_diff:.3g}")
    pd.DataFrame(checks).to_csv(os.path.join(ANNUAL_DIR, "baseline_identity_check.csv"),
                                index=False)

    for ref, cmp_ in itertools.permutations(MEMBERS, 2):
        print(f"Reference member {ref}, comparison member {cmp_} ...")
        frames = []
        for i, cid in enumerate(ids):
            low_thr = np.percentile(baseline[ref][:, i], RVA_LOW_PCTL)
            high_thr = np.percentile(baseline[ref][:, i], RVA_HIGH_PCTL)
            iha = compute_iha_period(baseline[cmp_][:, i], N_PRE_YEARS, low_thr, high_thr)
            iha.insert(0, "water_year", PRE_FIRST_WY + np.arange(len(iha)))
            iha.insert(0, "catchment_id", cid)
            frames.append(iha)
        path = os.path.join(ANNUAL_DIR, f"iha_annual_null_ref{ref}_cmp{cmp_}.csv.gz")
        pd.concat(frames, ignore_index=True).to_csv(path, index=False, float_format="%.10g",
                                                       compression={"method": "gzip", "mtime": 0})
        print(f"  {path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
