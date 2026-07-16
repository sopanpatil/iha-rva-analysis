# iha-rva-analysis

Indicators of Hydrologic Alteration (IHA) and Range of Variability Approach (RVA) analysis pipeline for Great Britain, comparing baseline (WY1982–WY2010) and future (WY2051–WY2080) HBV-simulated discharge across four RCP scenarios and four CHESS-SCAPE ensemble members.

This repository depends on the separate [`hbv-model`](https://github.com/sopanpatil/hbv-model) repository for the HBV model itself and its calibrated parameters.

## External data dependencies (not bundled in this repository)

- **CAMELS-GB v2**: catchment attributes and observed discharge. Publicly available from the NERC Environmental Data Service; see Coxon et al. (2026). The topographic attributes file (gauge latitude/longitude, used by the `figures/` scripts) is part of this dataset, as are the hydrologic, hydrogeology, and climatic attributes files used by `catchment_attribute_analysis.py` below. All four are keyed on `gauge_id`.
- **CHESS-SCAPE**: climate projection forcing (precipitation, temperature, surface wind, shortwave radiation, relative humidity). Publicly available from the Centre for Environmental Data Analysis; see Robinson et al. (2023).

## Contents

- `catchment_filter.py` — filters catchments by validation KGE ≥ 0.5 (621 of 671 retained), reading directly from `hbv-model`'s `calibrated_parameters.csv`. Its default path (`../hbv-model/calibrated_parameters.csv`) assumes `hbv-model` is cloned as a sibling folder alongside this repository; if not, pass `kge_path` explicitly. This is the single source of truth for the KGE threshold and file location — the `figures/` scripts import from here directly rather than duplicating it (see `figures/README.md`).
- `calculate_ensemble_mean_q.py` — computes the ensemble-mean HBV discharge (averaged across the 4 CHESS-SCAPE ensemble members) for each RCP, used as the primary analytical unit for the IHA/RVA analysis. Writes to `./ensemble_mean_discharge/`.
- `iha_categories.py` — single shared module defining the five-category directional hydrologic alteration (HA) classification scheme (high/moderate departure, low alteration, moderate/high concentration; see Methods for exact thresholds), used by `chess_scape_iha_analysis.py` and, indirectly, by `figures/shared_utils.py`.
- `chess_scape_iha_analysis.py` — computes all 32 IHA parameters and applies the RVA comparison (baseline vs. future) for a given RCP, either for a single ensemble member (`--ensemble`) or the ensemble mean (`--mean`, the primary analytical unit used in the paper). Outputs per-catchment HA values and five-category classifications, plus a parameter-level departure/concentration fraction summary, to `IHA_Results/`.
- `catchment_attribute_analysis.py` — implements the catchment attribute association analysis of Methods 3.5 / Results 4.5. For each of the 32 IHA parameters, each RCP, and eight CAMELS-GB v2 attributes, computes raw Pearson correlations with ensemble-mean HA, plus partial correlations controlling for aridity, each with Benjamini-Hochberg FDR correction. Reads the `rva_results_rcp{26,45,60,85}_ENSEMBLE.csv` files written by `chess_scape_iha_analysis.py --mean` above (default: `IHA_Results/`, so it needs no `--rva-dir` argument if run from the repository root after the main pipeline), and writes `raw_correlation_screen.csv` / `partial_correlation_screen.csv` to `IHA_Results/attribute_correlations/` by default — kept in a subfolder rather than alongside `rva_results_*`/`rva_summary_by_parameter_*` since it's a different unit of observation (one row per correlation test, not per catchment).
- `figures/` — all figure-generation scripts for the manuscript (Figures 1–5) and supplementary information (Figures S1–S12). See `figures/README.md` for details. These read the `rva_results_*.csv` files produced by `chess_scape_iha_analysis.py` above, so should be run after the main analysis pipeline.

## Running the pipeline

```bash
# 1. Compute the ensemble-mean discharge for each RCP (requires hbv-model's
#    calibrated HBV run to already have produced per-ensemble-member discharge)
python calculate_ensemble_mean_q.py --rcp rcp85

# 2. Run the IHA/RVA analysis on the ensemble mean (primary analysis)
python chess_scape_iha_analysis.py --rcp 85 --mean

# Or, to characterise ensemble uncertainty, run each member individually:
python chess_scape_iha_analysis.py --rcp 85 --ensemble 01
python chess_scape_iha_analysis.py --rcp 85 --ensemble 04
python chess_scape_iha_analysis.py --rcp 85 --ensemble 06
python chess_scape_iha_analysis.py --rcp 85 --ensemble 15
```

Repeat across all four RCPs (`--rcp 26`, `45`, `60`, `85`) for the complete analysis. This produces the `rva_results_rcp<rcp>_<label>.csv` files (both `ENSEMBLE` and all four members) that `figures/` and `catchment_attribute_analysis.py` read.

```bash
# 3. Test catchment attribute associations with ensemble-mean HA (Methods 3.5),
#    once all four RCPs' ENSEMBLE rva_results files exist in IHA_Results/
python catchment_attribute_analysis.py \
    --hydrologic-csv camels_gb_v2_hydrologic_attributes.csv \
    --hydrogeology-csv camels_gb_v2_hydrogeology_attributes.csv \
    --climatic-csv camels_gb_v2_climatic_attributes.csv \
    --topographic-csv camels_gb_v2_topographic_attributes.csv
```

```bash
# 4. Generate all manuscript and supplementary figures (run from figures/)
cd figures
python figure1_catchment_map.py
python figure2_rcp_divergence.py
python figures_3_4_5_main.py
python supplementary_figures_S1_S12.py
```

## Output

`IHA_Results/` (including its `attribute_correlations/` subfolder) is committed to this repository. At under 20 MB across 80-odd CSVs, it is the small, final, per-catchment result set that directly underpins every figure, table, and reported percentage in the manuscript — committing it means a tagged release of this repository (and its Zenodo archive) captures the exact code-and-results pair used in the paper, rather than requiring a separately versioned data record. If you re-run the pipeline locally, these files will be regenerated in place; check `git diff` before committing over them, since the ones currently in the repo are the specific run archived at the release DOI cited in the manuscript's Open Research section.

By contrast, `chess_scape_output/` and `ensemble_mean_discharge/` (the per-ensemble-member and ensemble-mean daily discharge time series) are excluded via `.gitignore`. These are multi-GB across the 16 RCP-ensemble-member combinations and fully regenerable from the public CHESS-SCAPE forcing plus this repository's code and `hbv-model`'s `calibrated_parameters.csv` — see the "software over output" principle in AGU's data/software guidance for the reasoning.

Each run of `chess_scape_iha_analysis.py` writes three CSVs to `IHA_Results/`:

- `rva_results_rcp<rcp>_<label>.csv` — one row per catchment, with `<param>_ha` (the hydrologic alteration value) and `<param>_cat` (five-category classification) for each of the 32 IHA parameters. `<label>` is `ENSEMBLE` for the ensemble-mean run, or the member ID (`01`/`04`/`06`/`15`) for an individual member.
- `iha_pre_rcp<rcp>_<label>.csv` / `iha_post_rcp<rcp>_<label>.csv` — mean IHA parameter values for the baseline and future periods respectively.
- `rva_summary_by_parameter_rcp<rcp>_<label>.csv` — parameter-level summary: mean HA, and the fraction of catchments in departure / high departure / concentration / high concentration, ranked by total departure fraction (matching the parameter ranking reported in the paper's Results).

`catchment_attribute_analysis.py` writes two further CSVs to `IHA_Results/attribute_correlations/` by default:

- `raw_correlation_screen.csv` — one row per (RCP, IHA parameter, attribute): raw Pearson `r`, `p`, and FDR-corrected `p_fdr`/`sig_fdr` (1,024 rows: 32 parameters x 8 attributes x 4 RCPs).
- `partial_correlation_screen.csv` — as above, but partial correlation `r_partial` controlling for aridity, for the 7 non-aridity attributes (896 rows: 32 parameters x 7 attributes x 4 RCPs).

`figures/` scripts write 300 dpi PNGs (Figures 1–5, S1–S12) to the path set in `OUTPUTS`/`OUT_PATH` near the top of each script.

## Requirements

```bash
pip install numpy pandas scipy statsmodels --break-system-packages
```

(`figures/` has its own additional requirements — see `figures/README.md`.)

## Citation

If you use this code, please cite the archived release (see `CITATION.cff`).
