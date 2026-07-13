# figures/

Figure-generation scripts for Patil et al., producing all figures in the main manuscript (Figures 1–5) and supplementary information (Figures S1–S12).

These scripts read the outputs of the main `iha-rva-analysis` pipeline (one level up), so run `chess_scape_iha_analysis.py` for all four RCPs and all five labels (`ENSEMBLE` + four members) before running anything here.

## Files

- `shared_utils.py` — shared constants, colour schemes, and helper functions (GB boundary/region loading, BNG coordinate conversion, ensemble agreement counting). **The five-category HA classification itself is imported directly from `../iha_categories.py`** rather than reimplemented here, so the figures are guaranteed to use exactly the same thresholds as the `rva_results_*.csv` files they plot.
- `figure1_catchment_map.py` — Figure 1: catchment locations coloured by validation KGE, with GB country/regional boundaries. Reads validation KGE from `hbv-model`'s `calibrated_parameters.csv` via `../catchment_filter.py` (same file, same column names, same 0.5 threshold as the main pipeline) rather than a separately maintained KGE file, so the retained/excluded counts shown here always match the 621/671 catchments actually used in the IHA/RVA analysis.
- `figure2_rcp_divergence.py` — Figure 2: RCP2.6 vs RCP8.5 discharge divergence check used to verify the scenario-neutral baseline period.
- `figures_3_4_5_main.py` — Figures 3, 4, and 5: IHA/RVA spatial alteration maps, ranked ensemble summary charts, and ensemble agreement maps, for four representative parameters. Reads `rva_results_{rcp}_ENSEMBLE.csv` and `rva_results_{rcp}_{member}.csv` from `../IHA_Results/`.
- `supplementary_figures_S1_S12.py` — Figures S1–S12: the complete spatial alteration maps (S1–S6) and ensemble agreement maps (S7–S12) for all 32 IHA parameters, organised by IHA group exactly as described in the Supplementary Information text.

## Requirements

```bash
pip install pandas numpy geopandas matplotlib pyproj shapely --break-system-packages
```

## Data expected

By default, scripts expect (edit the path constants near the top of each script if yours differ):

- `../../hbv-model/calibrated_parameters.csv` (via `catchment_filter.py`, for Figure 1's validation KGE)
- `camels_gb_v2_topographic_attributes.csv` (CAMELS-GB v2; gauge lat/lon for all spatial figures)
- `rcp26_hbv_ensemble_mean_discharge.csv`, `rcp85_hbv_ensemble_mean_discharge.csv` (Figure 2)
- `../IHA_Results/rva_results_{rcp}_ENSEMBLE.csv` and `../IHA_Results/rva_results_{rcp}_{member}.csv`, for `rcp` in `{rcp26, rcp45, rcp60, rcp85}` and `member` in `{01, 04, 06, 15}` (Figures 3–5, S1–S12) — these are written by `../chess_scape_iha_analysis.py`.

**Important:** the `rva_results_*.csv` files must be produced with `EXPECTED_FREQ = 0.5`, the WY1982–WY2010 baseline, and pulse thresholds computed once from the baseline period only (as `chess_scape_iha_analysis.py` currently does). You can verify you have the right files by checking that all `_ha` columns are bounded in [-1, +1]:

```python
import pandas as pd
df = pd.read_csv('../IHA_Results/rva_results_rcp85_ENSEMBLE.csv')
ha_cols = [c for c in df.columns if c.endswith('_ha')]
print(df[ha_cols].min().min(), df[ha_cols].max().max())  # should print -1.0 1.0
```

## Natural Earth boundary data

`figure1_catchment_map.py` and the Figures 3/4/5/S1–S12 scripts need Great Britain coastline and regional boundary data. The scripts will automatically download this from the Natural Earth GitHub mirror on first run and cache it locally (`_ne_10m_countries.geojson`, `_ne_10m_admin1.geojson`, `_gb_boundary.gpkg`, `_gb_regions.gpkg`). No manual download needed, but an internet connection is required the first time each script runs.

(The English-region name mapping in `get_gb_regions()` — e.g. Natural Earth's raw `'East'` region relabelled `'East Anglia'` — has been checked directly against the live Natural Earth admin-1 dataset field values, so all `region_map` keys are confirmed to match exactly.)

## A note on missing/insufficient-data catchments

`shared_utils.five_cat()` maps `NaN` HA values (catchments/parameters marked `"insufficient data"` in `rva_results_*.csv`, typically because too few valid years were available in the baseline or future period) to the `'low'` (grey) category purely so the spatial maps render cleanly rather than erroring or leaving blank pixels. This is a plotting convenience only — it is not a claim that these catchments show low alteration. If this affects more than a handful of catchments for any given parameter, the departure/concentration fractions shown in Figure 4 and computed by `compute_fractions()` will be marginally understated, since these catchments are counted in the denominator but not in either the departure or concentration category. Worth a quick cross-check against the `n_catchments` used in `../IHA_Results/rva_summary_by_parameter_*.csv` if exact match to the manuscript's reported fractions is needed.

## Running

Each script can be run standalone from within `figures/`:

```bash
python figure1_catchment_map.py
python figure2_rcp_divergence.py
python figures_3_4_5_main.py
python supplementary_figures_S1_S12.py
```

All outputs are saved as 300 dpi PNG files to the path set in `OUTPUTS`/`OUT_PATH` near the top of each script (defaults to `/mnt/user-data/outputs/` — change this to a local directory when running outside the Claude environment).
