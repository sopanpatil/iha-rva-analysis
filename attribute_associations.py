"""
attribute_associations.py
=========================

Association between catchment attributes and the magnitude and direction of
projected change in each IHA parameter.

The response variable is Cliff's delta between future and baseline water
years (median of the four ensemble members; rva_member_analysis.py), which is
signed and bounded in [-1, 1]. The predictors are four weakly correlated
CAMELS-GB v2 attributes:

  aridity         PET / P; stands for the climatic-wetness gradient, with
                  which runoff ratio, drainage path slope, mean elevation and
                  precipitation seasonality are all strongly correlated
                  (|Spearman rho| >= 0.80 across the retained catchments)
  frac_high_perc  fraction of catchment area on highly productive (high
                  permeability) geology
  area            catchment area (log10 km2)
  frac_snow       fraction of precipitation falling as snow

For each RCP, IHA parameter and attribute the script reports the Spearman rank
correlation and, for the three non-aridity attributes, the partial Spearman
correlation controlling for aridity (Pearson correlation of rank residuals).
p-values are adjusted with the Benjamini-Hochberg procedure within each RCP.

Outputs (IHA_Results/attributes/)
-------
  attribute_rank_correlations.csv   rcp, parameter, attribute, rho, p, p_fdr,
                                    partial_rho, partial_p, partial_p_fdr
  attribute_intercorrelation.csv    Spearman correlation matrix of the eight
                                    attributes used in the original screen

Run
---
    python attribute_associations.py --attributes-dir /path/to/camels_gb_v2_attributes
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

ENS_CSV = os.path.join("IHA_Results", "member_rva", "rva_catchment_ensemble.csv.gz")
OUT_DIR = os.path.join("IHA_Results", "attributes")
PREDICTORS = ["aridity", "frac_high_perc", "log_area", "frac_snow"]
SCREEN_ATTRS = ["aridity", "runoff_ratio", "dpsbar", "elev_mean", "p_seasonality",
                "frac_high_perc", "log_area", "frac_snow"]
FILES = {
    "hydrologic": "camels_gb_v2_hydrologic_attributes.csv",
    "hydrogeology": "camels_gb_v2_hydrogeology_attributes.csv",
    "climatic": "camels_gb_v2_climatic_attributes.csv",
    "topographic": "camels_gb_v2_topographic_attributes.csv",
}


def load_attributes(directory: str) -> pd.DataFrame:
    tables = [pd.read_csv(os.path.join(directory, f)).set_index("gauge_id") for f in FILES.values()]
    attrs = pd.concat(tables, axis=1)
    attrs = attrs.loc[:, ~attrs.columns.duplicated()]
    attrs["log_area"] = np.log10(attrs["area"])
    return attrs


def partial_spearman(y, x, z):
    """Spearman correlation of y and x controlling for z (rank residuals)."""
    ry, rx, rz = rankdata(y), rankdata(x), rankdata(z)
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)
    return pearsonr(res(ry), res(rx))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--attributes-dir", required=True,
                    help="Directory holding the four CAMELS-GB v2 attribute CSV files.")
    args = ap.parse_args()

    attrs = load_attributes(args.attributes_dir)
    ens = pd.read_csv(ENS_CSV)
    ids = attrs.index.intersection(ens.catchment_id.unique())
    os.makedirs(OUT_DIR, exist_ok=True)
    attrs.loc[ids, SCREEN_ATTRS].corr("spearman").round(3).to_csv(
        os.path.join(OUT_DIR, "attribute_intercorrelation.csv"))

    rows = []
    for (rcp, param), g in ens.groupby(["rcp", "parameter"]):
        y = g.set_index("catchment_id").cliffs_delta_median
        d = pd.concat([y.rename("y"), attrs[PREDICTORS]], axis=1, join="inner").dropna()
        for a in PREDICTORS:
            rho, p = spearmanr(d.y, d[a])
            row = {"rcp": rcp, "parameter": param, "attribute": a, "n": len(d), "rho": rho, "p": p}
            if a != "aridity":
                row["partial_rho"], row["partial_p"] = partial_spearman(d.y.values, d[a].values,
                                                                        d.aridity.values)
            rows.append(row)
    out = pd.DataFrame(rows)
    for rcp, idx in out.groupby("rcp").groups.items():
        out.loc[idx, "p_fdr"] = multipletests(out.loc[idx, "p"], method="fdr_bh")[1]
        sub = out.loc[idx].dropna(subset=["partial_p"])
        out.loc[sub.index, "partial_p_fdr"] = multipletests(sub.partial_p, method="fdr_bh")[1]
    out.round(4).to_csv(os.path.join(OUT_DIR, "attribute_rank_correlations.csv"), index=False)
    print(f"Outputs written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
