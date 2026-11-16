# figures/

Scripts for all figures in the manuscript and the Supplementary Information, and the data for Table 1 and Supplementary Tables S1 and S2. They read the results in `../IHA_Results/`, so the analysis pipeline (see the main `README.md`) must be run first, or the committed results used.

## Files

| Script | Figures | Output files |
|---|---|---|
| `figure1_catchment_map.py` | Figure 1 | `figure1_catchment_map.png` |
| `benchmark_climate_figures.py` | Figures 2–4, Supplementary Fig. S20; Table 1 and Supplementary Tables S1–S2 | `fig_benchmark`, `fig_member_climate`, `fig_change_maps`, `figS_change_maps_all`; `tables/table1_members.csv`, `tables/benchmark_lowflow.csv`, `tables/benchmark_rva.csv` |
| `alteration_figures.py` | Figures 5–10, Supplementary Fig. S18 | `fig_direction_maps`, `fig_parameter_summary`, `fig_excess_departure`, `fig_effect_size`, `fig_agreement_maps`, `fig_attribute_scatter`, `fig_ensemble_mean_vs_members` |
| `supplementary_maps.py` | Supplementary Figs. S2–S17, S19 | `figS_direction_<group>`, `figS_agreement_<group>`, `figS_members_<parameter>`, `figS_member_spread` |
| `scenario_divergence.py` | Supplementary Fig. S1 | `figS_scenario_divergence.png` |
| `shared_utils.py` | Shared constants, colours and map helpers | |

`alteration_figures.py` defines the colours, labels and map layout used by `benchmark_climate_figures.py` and `supplementary_maps.py`.

## Inputs

- `../IHA_Results/` (committed with this repository).
- CAMELS-GB v2 attribute files (`camels_gb_v2_topographic_attributes.csv` for the gauge locations, and `camels_gb_v2_climatic_attributes.csv` for Figure 10), in the folder given by the `ATTR_DIR` environment variable (default: this folder).
- Figure 1: `../../hbv-model/calibrated_parameters.csv` (validation KGE).
- Supplementary Fig. S1: the ensemble-mean discharge for RCP2.6 and RCP8.5 in `../ensemble_mean_discharge/` (written by `../calculate_ensemble_mean_q.py`).

The Great Britain coastline and regional boundaries come from Natural Earth. They are downloaded on the first run and cached in this folder (`_ne_10m_countries.geojson`, `_ne_10m_admin1.geojson`, `_gb_boundary.gpkg`, `_gb_regions.gpkg`), so the first run needs an internet connection.

## Running

From this folder:

```bash
export ATTR_DIR=/path/to/CAMELS-GB/attributes
python figure1_catchment_map.py
python benchmark_climate_figures.py
python alteration_figures.py
python supplementary_maps.py
python scenario_divergence.py
```

Figures are saved as 300 dpi PNG files in the folder given by `FIG_OUT` (default: `output/`). With `FIG_PDF=1`, `alteration_figures.py` and `benchmark_climate_figures.py` also save vector PDF copies.

## Requirements

```bash
pip install numpy pandas matplotlib geopandas pyproj shapely
```
