"""
catchment_attribute_analysis.py
================================

Catchment attribute association analysis described in Section 3.5 of
Patil et al., "Climate change produces sharply uneven responses across
the full flow regime of British catchments" (submitted to Water
Resources Research).

For each of the 32 IHA parameters, each RCP scenario, and a set of eight
CAMELS-GB v2 catchment attributes, this script computes:
  (1) raw Pearson correlations between ensemble-mean hydrologic
      alteration (HA) and each attribute, and
  (2) partial correlations controlling for aridity, to isolate
      associations independent of the northwest-southeast climatic
      gradient described in Section 2.

The Benjamini-Hochberg false discovery rate procedure (alpha = 0.05) is
applied across the full set of tests in each case, matching the
1,024 raw tests (32 parameters x 8 attributes x 4 RCPs) and 896 partial
tests (32 parameters x 7 non-aridity attributes x 4 RCPs) reported in
Section 4.5.

Catchment filtering
--------------------
This script does not apply the validation-KGE filter itself. The
rva_results_rcp<rcp>_ENSEMBLE.csv files produced by
chess_scape_iha_analysis.py already contain only the 621 catchments
retained via catchment_filter.py (validation KGE >= 0.5; Methods 3.2),
so the inner join against the CAMELS-GB v2 attribute tables (which cover
all 671 catchments) naturally restricts every test to the same 621
catchments used throughout the rest of the analysis. main() prints the
merged catchment count per RCP so this can be checked directly rather
than assumed silently.

Inputs
------
  <rva-dir>/rva_results_rcp{26,45,60,85}_ENSEMBLE.csv   (from chess_scape_iha_analysis.py --mean)
  CAMELS-GB v2 attribute tables (hydrologic, hydrogeology, climatic,
  topographic), each keyed on `gauge_id`, matching the convention used
  by figures/figure1_catchment_map.py and figures/figures_3_4_5_main.py.

Output
------
  {out-dir}/raw_correlation_screen.csv       : all raw correlation tests
  {out-dir}/partial_correlation_screen.csv   : all partial correlation tests

Summary statistics and the frac_high_perc parameter ranking (Section 4.5)
are printed to stdout.

Run
---
    python catchment_attribute_analysis.py \\
        --hydrologic-csv camels_gb_v2_hydrologic_attributes.csv \\
        --hydrogeology-csv camels_gb_v2_hydrogeology_attributes.csv \\
        --climatic-csv camels_gb_v2_climatic_attributes.csv \\
        --topographic-csv camels_gb_v2_topographic_attributes.csv

    # Non-default locations for the RVA results / correlation outputs:
    python catchment_attribute_analysis.py \\
        --rva-dir /path/to/IHA_Results \\
        --hydrologic-csv camels_gb_v2_hydrologic_attributes.csv \\
        --hydrogeology-csv camels_gb_v2_hydrogeology_attributes.csv \\
        --climatic-csv camels_gb_v2_climatic_attributes.csv \\
        --topographic-csv camels_gb_v2_topographic_attributes.csv \\
        --out-dir /path/to/IHA_Results/attribute_correlations

Requires numpy, pandas, scipy, and statsmodels, plus
chess_scape_iha_analysis.py having already been run for all four RCPs
with --mean (see this repository's README).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

RCPS = ["26", "45", "60", "85"]

# Default location of the rva_results_rcp<rcp>_ENSEMBLE.csv files, matching
# chess_scape_iha_analysis.py's OUT_DIR, so this script can be run with no
# arguments beyond the four attribute CSVs once the main pipeline has been
# run. Correlation-screen outputs are written to a subfolder of the same
# directory (attribute_correlations/) rather than alongside rva_results_*/
# rva_summary_by_parameter_* directly: those files are one row per catchment,
# whereas raw_correlation_screen.csv/partial_correlation_screen.csv are one
# row per (RCP, parameter, attribute) correlation test -- a different unit
# of observation that is easy to confuse with the per-catchment files if the
# two are interleaved in a single directory listing.
DEFAULT_RVA_DIR = "IHA_Results"
DEFAULT_OUT_DIR = "IHA_Results/attribute_correlations"

# The eight CAMELS-GB v2 attributes used in Section 3.5. Where CAMELS-GB
# v2 provides multiple closely related measures of the same underlying
# characteristic (e.g., several baseflow index formulations, several
# elevation summaries), only one representative attribute is retained
# per characteristic to avoid redundant tests. Note that the model's
# own simulated base flow index (Section 3.4) is deliberately excluded
# here: it is a dynamic, per-water-year IHA parameter, not a static
# catchment attribute, and is tested separately (see Section 4.5).
ARIDITY_COL = "aridity"
NON_ARIDITY_ATTRS = [
    "frac_high_perc",   # fraction of catchment area on high-permeability geology
    "runoff_ratio",
    "area",
    "dpsbar",           # drainage path slope
    "elev_mean",
    "p_seasonality",
    "frac_snow",
]
ALL_ATTRS = NON_ARIDITY_ATTRS + [ARIDITY_COL]


def load_attributes(hydrologic_csv, hydrogeology_csv, climatic_csv, topographic_csv):
    """Merge the four CAMELS-GB v2 attribute files on gauge_id -> catchment_id."""
    hydro = pd.read_csv(hydrologic_csv)
    geol = pd.read_csv(hydrogeology_csv)
    clim = pd.read_csv(climatic_csv)
    topo = pd.read_csv(topographic_csv)

    attrs = hydro.merge(geol, on="gauge_id").merge(clim, on="gauge_id").merge(topo, on="gauge_id")
    attrs = attrs.rename(columns={"gauge_id": "catchment_id"})

    missing = [c for c in ALL_ATTRS if c not in attrs.columns]
    if missing:
        raise ValueError(f"Expected attribute columns not found: {missing}")

    return attrs[["catchment_id"] + ALL_ATTRS]


def load_ha_long(rva_dir):
    """
    Load ensemble-mean rva_results for all four RCPs and reshape to long
    format: one row per (rcp, catchment_id, parameter, ha_value).
    """
    rva_dir = Path(rva_dir)
    frames = []
    for rcp in RCPS:
        path = rva_dir / f"rva_results_rcp{rcp}_ENSEMBLE.csv"
        if not path.exists():
            raise FileNotFoundError(f"Expected file not found: {path}")
        df = pd.read_csv(path)
        ha_cols = [c for c in df.columns if c.endswith("_ha")]
        long = df.melt(
            id_vars="catchment_id",
            value_vars=ha_cols,
            var_name="parameter",
            value_name="ha",
        )
        long["parameter"] = long["parameter"].str.replace("_ha$", "", regex=True)
        long["rcp"] = rcp
        frames.append(long)
    return pd.concat(frames, ignore_index=True)


def partial_corr(x, y, z):
    """Pearson correlation between x and y after linearly residualising both on z."""
    bx = np.polyfit(z, x, 1)
    by = np.polyfit(z, y, 1)
    rx = x - np.polyval(bx, z)
    ry = y - np.polyval(by, z)
    return pearsonr(rx, ry)


def run_raw_screen(attrs, ha_long, min_n=50):
    """Raw Pearson correlation of HA against every attribute (including aridity)."""
    records = []
    merged = attrs.merge(ha_long, on="catchment_id", how="inner")
    for rcp in RCPS:
        for param, group in merged[merged["rcp"] == rcp].groupby("parameter"):
            for attr in ALL_ATTRS:
                sub = group[[attr, "ha"]].dropna()
                if len(sub) < min_n:
                    continue
                r, p = pearsonr(sub[attr], sub["ha"])
                records.append({"rcp": rcp, "parameter": param, "attribute": attr, "r": r, "p": p})
    df = pd.DataFrame(records)
    _, p_fdr, _, _ = multipletests(df["p"], alpha=0.05, method="fdr_bh")
    df["p_fdr"] = p_fdr
    df["sig_fdr"] = p_fdr < 0.05
    return df


def run_partial_screen(attrs, ha_long, min_n=50):
    """Partial Pearson correlation of HA against each non-aridity attribute,
    controlling for aridity."""
    records = []
    merged = attrs.merge(ha_long, on="catchment_id", how="inner")
    for rcp in RCPS:
        for param, group in merged[merged["rcp"] == rcp].groupby("parameter"):
            for attr in NON_ARIDITY_ATTRS:
                sub = group[[ARIDITY_COL, attr, "ha"]].dropna()
                if len(sub) < min_n:
                    continue
                r, p = partial_corr(sub[attr].values, sub["ha"].values, sub[ARIDITY_COL].values)
                records.append({"rcp": rcp, "parameter": param, "attribute": attr, "r_partial": r, "p": p})
    df = pd.DataFrame(records)
    _, p_fdr, _, _ = multipletests(df["p"], alpha=0.05, method="fdr_bh")
    df["p_fdr"] = p_fdr
    df["sig_fdr"] = p_fdr < 0.05
    return df


def print_summary(raw_df, partial_df):
    n_raw = len(raw_df)
    n_raw_sig = int(raw_df["sig_fdr"].sum())
    n_partial = len(partial_df)
    n_partial_sig = int(partial_df["sig_fdr"].sum())

    print(f"Raw correlations: {n_raw_sig}/{n_raw} significant after FDR correction "
          f"({100 * n_raw_sig / n_raw:.1f}%)")
    print(f"Partial correlations (controlling for aridity): {n_partial_sig}/{n_partial} "
          f"significant after FDR correction ({100 * n_partial_sig / n_partial:.1f}%)")

    print("\nRanking of all 32 parameters by |partial r| with frac_high_perc, by RCP:")
    for rcp in RCPS:
        sub = partial_df[(partial_df["attribute"] == "frac_high_perc") & (partial_df["rcp"] == rcp)].copy()
        sub["abs_r"] = sub["r_partial"].abs()
        sub = sub.sort_values("abs_r", ascending=False).reset_index(drop=True)
        sub["rank"] = sub.index + 1
        for param in ["low_pulse_count", "high_pulse_count"]:
            row = sub[sub["parameter"] == param]
            if not row.empty:
                row = row.iloc[0]
                print(f"  RCP{rcp}: {param:18s} rank={int(row['rank']):2d}  r={row['r_partial']:.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rva-dir", default=DEFAULT_RVA_DIR,
                        help=f"Directory containing rva_results_rcp{{26,45,60,85}}_ENSEMBLE.csv "
                             f"(default: {DEFAULT_RVA_DIR}, matching chess_scape_iha_analysis.py's output directory)")
    parser.add_argument("--hydrologic-csv", required=True, help="CAMELS-GB v2 hydrologic attributes CSV (keyed on gauge_id)")
    parser.add_argument("--hydrogeology-csv", required=True, help="CAMELS-GB v2 hydrogeology attributes CSV (keyed on gauge_id)")
    parser.add_argument("--climatic-csv", required=True, help="CAMELS-GB v2 climatic attributes CSV (keyed on gauge_id)")
    parser.add_argument("--topographic-csv", required=True, help="CAMELS-GB v2 topographic attributes CSV (keyed on gauge_id)")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR,
                        help=f"Directory to write output CSVs (default: {DEFAULT_OUT_DIR})")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading CAMELS-GB v2 attributes...")
    attrs = load_attributes(args.hydrologic_csv, args.hydrogeology_csv, args.climatic_csv, args.topographic_csv)

    print("Loading ensemble-mean HA values...")
    ha_long = load_ha_long(args.rva_dir)

    print("\nMerged catchment count by RCP (attributes inner-joined with HA values;")
    print("expect 621 for every RCP -- see 'Catchment filtering' in the module docstring):")
    merged_check = attrs.merge(ha_long, on="catchment_id", how="inner")
    for rcp in RCPS:
        n = merged_check.loc[merged_check["rcp"] == rcp, "catchment_id"].nunique()
        flag = "" if n == 621 else "  <-- unexpected count, check inputs"
        print(f"  RCP{rcp}: {n} catchments{flag}")
    del merged_check

    print("\nRunning raw correlation screen (32 parameters x 8 attributes x 4 RCPs)...")
    raw_df = run_raw_screen(attrs, ha_long)
    raw_df.to_csv(out_dir / "raw_correlation_screen.csv", index=False)

    print("Running aridity-controlled partial correlation screen "
          "(32 parameters x 7 non-aridity attributes x 4 RCPs)...")
    partial_df = run_partial_screen(attrs, ha_long)
    partial_df.to_csv(out_dir / "partial_correlation_screen.csv", index=False)

    print()
    print_summary(raw_df, partial_df)


if __name__ == "__main__":
    sys.exit(main())
