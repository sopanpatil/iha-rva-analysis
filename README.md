# iha-rva-analysis

Indicators of Hydrologic Alteration (IHA) and Range of Variability Approach (RVA) analysis of projected flow-regime change in Great Britain. The pipeline compares baseline (WY1982–WY2010) and future (WY2051–WY2080) HBV-simulated daily discharge for 621 CAMELS-GB v2 catchments, forced with four CHESS-SCAPE ensemble members under four RCP scenarios.

This repository contains the code and the per-catchment results behind Patil, Dallison and Jahanshahi, *Climate change drives sharply uneven responses across the full flow regime of British catchments*.

## Version 2.0.0

Version 2.0.0 accompanies the revised manuscript and changes the analysis as follows:

- **Per-member analysis.** IHA parameters and RVA statistics are computed separately for each ensemble member against its own baseline, then summarised across the four members. Version 1 averaged discharge across members before the analysis.
- **Leave-one-out expected frequency.** The expected in-range frequency of the RVA is the leave-one-out coverage of the baseline, instead of the conventional 0.5, which is biased for finite baselines and for count parameters with tied values.
- **No-change reference.** The 12 comparisons between the baseline periods of different ensemble members give the level of each statistic expected without forced change.
- **Direction and effect size.** Each departure is split into shifts below and above the target range, and the magnitude of change is given by Cliff's delta. The days of the annual minimum and maximum flow are treated as circular quantities.
- **Benchmark against observations.** The modelling chain is evaluated over the baseline period against observed flows, including its low-flow performance.
- **Catchment attributes.** Rank (Spearman) and partial rank correlations of Cliff's delta with four weakly correlated attributes replace the Pearson screen of HA against eight attributes.
- **CO2 sensitivity.** The RCP8.5 projections are repeated with potential evapotranspiration adjusted for the stomatal response of vegetation to rising CO2.

The version 1 scripts that these replace (`catchment_attribute_analysis.py`, `figures/figures_3_4_5_main.py`, `figures/supplementary_figures_S1_S12.py`) have been removed; they remain available in the archived v1.0.1 release (https://doi.org/10.5281/zenodo.21862874).

## External data and code (not included)

- **hbv-model** (https://doi.org/10.5281/zenodo.21406134): the HBV model and the calibrated parameters (`calibrated_parameters.csv`). The scripts expect the repository to be cloned alongside this one (`../hbv-model/`). `benchmark_observed.py` and `co2_sensitivity_runs.py` run HBV and need the `hbv_model` package on the Python path (`export PYTHONPATH=../hbv-model`).
- **HBV output for CHESS-SCAPE** (https://doi.org/10.5281/zenodo.21861386): catchment-mean CHESS-SCAPE forcing and the HBV discharge and soil moisture for all RCPs and ensemble members. Place the files in `chess_scape_output/` at the repository root.
- **CAMELS-GB v2** (Coxon et al., 2026): catchment attributes (topographic, climatic, hydrologic and hydrogeology files, keyed on `gauge_id`) and daily observed time series. Available from the NERC Environmental Data Service.
- **CHESS-SCAPE** (Robinson et al., 2023): the climate projections, available from the Centre for Environmental Data Analysis. They are only needed to repeat the extraction of the catchment-mean forcing, which is part of the HBV output record above.

## Contents

Core modules

| File | Purpose |
|---|---|
| `catchment_filter.py` | Retains catchments with validation KGE ≥ 0.5 (621 of 671), read from `../hbv-model/calibrated_parameters.csv`. |
| `iha_categories.py` | The five-category classification of the hydrologic alteration factor (HA): high and moderate departure, low alteration, moderate and high concentration. |

Analysis scripts, in the order they are run

| Script | Purpose | Outputs |
|---|---|---|
| `calculate_ensemble_mean_q.py` | Ensemble-mean discharge per RCP (used only for the comparison with the member-wise analysis and for Supplementary Fig. S1). | `ensemble_mean_discharge/` |
| `chess_scape_iha_analysis.py` | The 32 IHA parameters for every water year of both periods, for one RCP and one ensemble member (`--ensemble`) or the ensemble mean (`--mean`), and the version 1 RVA tables (expected frequency 0.5). | `IHA_Results/annual/`, `IHA_Results/*.csv` |
| `iha_null_member_pairs.py` | Baseline IHA parameters of each member computed with the pulse thresholds of every other member, for the no-change reference; checks that the baseline is identical across RCPs. | `IHA_Results/annual/` |
| `rva_member_analysis.py` | RVA statistics per RCP, member, catchment and parameter (leave-one-out expected frequency, direction, median change, Cliff's delta), the no-change reference, and their summaries. | `IHA_Results/member_rva/` |
| `ensemble_mean_rva.py` | Departure and concentration from ensemble-mean discharge with the same method, for comparison with the member-wise analysis (Supplementary Fig. S18). | `IHA_Results/member_rva/` |
| `nochange_expected_frequency.py` | Mean HA without forced change for three choices of the expected frequency (0.5, in-sample, leave-one-out). | `IHA_Results/member_rva/` |
| `benchmark_observed.py` | Runs HBV with observed CAMELS-GB v2 forcing: KGE reproduction, low-flow metrics, and baseline IHA parameters for observed flows and each simulation. | `IHA_Results/benchmark/` |
| `benchmark_analysis.py` | RVA benchmark of the modelling chain (hydrological model, climate forcing, full chain). | `IHA_Results/benchmark/` |
| `lowflow_robustness.py` | Interannual low-flow agreement with observations, and projected change in catchment subsets with good low-flow performance. | `IHA_Results/benchmark/` |
| `extract_period_means.py` | Baseline and future means of precipitation, temperature, PET, discharge, actual evapotranspiration and soil moisture. | `IHA_Results/period_means/` |
| `attribute_associations.py` | Spearman and partial Spearman correlations of Cliff's delta with catchment attributes. | `IHA_Results/attributes/` |
| `co2_sensitivity_runs.py` | RCP8.5 HBV runs with CO2-adjusted PET. | `IHA_Results/co2_sensitivity*/` |
| `co2_sensitivity_analysis.py` | Alteration and water balance with the original and CO2-adjusted PET. | `IHA_Results/co2_sensitivity*/` |

`figures/` contains the scripts for all figures in the manuscript and the Supplementary Information; see `figures/README.md`.

## Reproducing the results

Run from the repository root, with `hbv-model` cloned alongside it and the HBV output in `chess_scape_output/`. Each script describes its inputs, outputs and options in its header.

```bash
# 1. Yearly IHA parameters for all RCPs, members and the ensemble mean
for rcp in 26 45 60 85; do
    python calculate_ensemble_mean_q.py --rcp rcp$rcp
    python chess_scape_iha_analysis.py --rcp $rcp --mean
    for ens in 01 04 06 15; do
        python chess_scape_iha_analysis.py --rcp $rcp --ensemble $ens
    done
done
python iha_null_member_pairs.py --rcp 85

# 2. Per-member RVA statistics and the no-change reference
python rva_member_analysis.py
python ensemble_mean_rva.py
python nochange_expected_frequency.py

# 3. Benchmark against observations (needs the hbv_model package)
export PYTHONPATH=../hbv-model
python benchmark_observed.py \
    --camels-pattern "/path/to/CAMELS-GB/camels_gb_v2_hydromet_daily_timeseries_{id}_*.csv" \
    --date-col date --precip-col precipitation_cehgear --temp-col temperature_chess \
    --pet-col pet_chess --q-col discharge_spec --run-start 1970-10-01
python benchmark_analysis.py
python lowflow_robustness.py

# 4. Changes in the forcing and model states
python extract_period_means.py

# 5. Catchment attributes
python attribute_associations.py --attributes-dir /path/to/CAMELS-GB/attributes

# 6. CO2 sensitivity (RCP8.5; needs the hbv_model package)
TOPO=/path/to/CAMELS-GB/attributes/camels_gb_v2_topographic_attributes.csv
python co2_sensitivity_runs.py --topo-csv $TOPO --u2 2
python co2_sensitivity_runs.py --topo-csv $TOPO --u2 4 --out-dir IHA_Results/co2_sensitivity_u2_4
python co2_sensitivity_analysis.py
python co2_sensitivity_analysis.py --co2-dir IHA_Results/co2_sensitivity_u2_4 --run-name co2pet_u2_4

# 7. Figures and tables
cd figures
export ATTR_DIR=/path/to/CAMELS-GB/attributes
python figure1_catchment_map.py
python benchmark_climate_figures.py
python alteration_figures.py
python supplementary_maps.py
python scenario_divergence.py
```

`co2_sensitivity_runs.py` reads RCP8.5 CO2 concentrations from `--co2-file` (`RCP85_MIDYEAR_CONCENTRATIONS.DAT` from the RCP database) if given, and otherwise interpolates the decadal values embedded in the script.

## Results

`IHA_Results/` holds the results used in the manuscript, so that a tagged release captures the exact code and results reported. Re-running the pipeline regenerates these files in place.

| Folder | Contents |
|---|---|
| `annual/` | `iha_annual_rcp<rcp>_<member>.csv.gz`: the 32 IHA parameters for every catchment and water year of the baseline and future periods (`<member>` is `01`, `04`, `06`, `15` or `ENSEMBLE`). `iha_annual_null_ref<A>_cmp<B>.csv.gz`: baseline IHA parameters of member B with the pulse thresholds of member A. `baseline_identity_check.csv`. |
| `member_rva/` | `rva_catchment_member.csv.gz` (one row per RCP, member, catchment and parameter), `rva_catchment_ensemble.csv.gz` (summarised across members), `rva_null_catchment.csv.gz` (no-change comparisons), `rva_summary_by_parameter.csv` (percentages of catchments per category, direction and magnitude, with the no-change levels), `rva_summary_ensemble_mean.csv` (the same from ensemble-mean discharge), `nochange_expected_frequency.csv`. |
| `benchmark/` | KGE reproduction, low-flow metrics, baseline IHA parameters of observed and simulated flows, the RVA benchmark and the low-flow robustness tables. |
| `period_means/` | Baseline and future means of the forcing and model states. |
| `attributes/` | Attribute intercorrelation and the rank correlations with Cliff's delta. |
| `co2_sensitivity/`, `co2_sensitivity_u2_4/` | CO2-adjusted runs for wind speeds of 2 and 4 m/s, with their summaries and a check that the control run reproduces the stored discharge. |
| top level | The per-run RVA tables in the version 1 format (`rva_results_*`, `iha_pre_*`, `iha_post_*`, `rva_summary_by_parameter_*`; expected frequency 0.5), still written by `chess_scape_iha_analysis.py`. The analysis in the manuscript uses the tables in `annual/` and `member_rva/`. |

`chess_scape_output/` and `ensemble_mean_discharge/` are not included (see `.gitignore`).

## Requirements

```bash
pip install numpy pandas scipy statsmodels
pip install matplotlib geopandas pyproj shapely   # figures
```

## Citation

If you use this code, please cite the archived release (see `CITATION.cff`).

## Licence

MIT (see `LICENSE`).
