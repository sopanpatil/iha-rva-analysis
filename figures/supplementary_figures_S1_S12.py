"""
supplementary_figures_S1_S12.py
================================
Supplementary Figures S1-S12: complete spatial alteration maps (S1-S6) and
ensemble agreement maps (S7-S12) for all 32 IHA parameters, organised by
IHA group, across all four RCP scenarios.

Inputs (edit paths as needed):
    rva_results_{rcp}_ENSEMBLE.csv        for rcp in [rcp26, rcp45, rcp60, rcp85]
    rva_results_{rcp}_{member}.csv        for member in [01, 04, 06, 15]
    camels_gb_v2_topographic_attributes.csv

Outputs:
    figureS1_group1_summer.png   ... figureS12_agreement_group5.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from shared_utils import (five_cat, CAT_COLOURS, CAT_ORDER, CAT_LABELS,
                           AGR_COLOURS, AGR_LABELS, RCPS, RCP_LABELS, MEMBERS,
                           XLIM, YLIM, MAP_KW, get_gb_boundary, load_topo_with_bng,
                           compute_agreement)

# chess_scape_iha_analysis.py writes rva_results_*.csv to IHA_Results/ at the
# repository root; this script lives in figures/, hence the one level up.
RVA_DIR = '../IHA_Results/'
UPLOADS = '/mnt/user-data/uploads/'
OUTPUTS = '/mnt/user-data/outputs/'

# Parameter groupings for the six spatial-map / six agreement-map figures
SUPP_GROUPS = [
    ('group1_summer', [
        ('month_08_May_ha', 'May', 'Group 1'), ('month_09_Jun_ha', 'June', 'Group 1'),
        ('month_10_Jul_ha', 'July', 'Group 1'), ('month_11_Aug_ha', 'August', 'Group 1'),
        ('month_12_Sep_ha', 'September', 'Group 1'), ('month_01_Oct_ha', 'October', 'Group 1')]),
    ('group1_winter', [
        ('month_02_Nov_ha', 'November', 'Group 1'), ('month_03_Dec_ha', 'December', 'Group 1'),
        ('month_04_Jan_ha', 'January', 'Group 1'), ('month_05_Feb_ha', 'February', 'Group 1'),
        ('month_06_Mar_ha', 'March', 'Group 1'), ('month_07_Apr_ha', 'April', 'Group 1')]),
    ('group2_minima', [
        ('min_1day_ha', '1-day minimum', 'Group 2'), ('min_3day_ha', '3-day minimum', 'Group 2'),
        ('min_7day_ha', '7-day minimum', 'Group 2'), ('min_30day_ha', '30-day minimum', 'Group 2'),
        ('min_90day_ha', '90-day minimum', 'Group 2'), ('base_flow_index_ha', 'Base flow index', 'Group 2')]),
    ('group2_maxima', [
        ('max_1day_ha', '1-day maximum', 'Group 2'), ('max_3day_ha', '3-day maximum', 'Group 2'),
        ('max_7day_ha', '7-day maximum', 'Group 2'), ('max_30day_ha', '30-day maximum', 'Group 2'),
        ('max_90day_ha', '90-day maximum', 'Group 2')]),
    ('groups3_4', [
        ('doy_min_ha', 'Day of minimum flow', 'Group 3'), ('doy_max_ha', 'Day of maximum flow', 'Group 3'),
        ('low_pulse_count_ha', 'Low pulse count', 'Group 4'),
        ('low_pulse_duration_mean_ha', 'Low pulse duration', 'Group 4'),
        ('high_pulse_count_ha', 'High pulse count', 'Group 4'),
        ('high_pulse_duration_mean_ha', 'High pulse duration', 'Group 4')]),
    ('group5', [
        ('mean_rise_rate_ha', 'Mean rise rate', 'Group 5'), ('mean_fall_rate_ha', 'Mean fall rate', 'Group 5'),
        ('number_of_reversals_ha', 'Number of reversals', 'Group 5')]),
]


def load_data():
    topo = load_topo_with_bng(UPLOADS + 'camels_gb_v2_topographic_attributes.csv')

    ens_dfs = {}
    for rcp in RCPS:
        df = pd.read_csv(f'{RVA_DIR}rva_results_{rcp}_ENSEMBLE.csv')
        ens_dfs[rcp] = df.merge(topo[['gauge_id', 'easting', 'northing']],
                                 left_on='catchment_id', right_on='gauge_id')

    mem_dfs = {}
    for rcp in RCPS:
        mem_dfs[rcp] = {}
        for mem in MEMBERS:
            df = pd.read_csv(f'{RVA_DIR}rva_results_{rcp}_{mem}.csv')
            df = df.merge(topo[['gauge_id', 'easting', 'northing']],
                           left_on='catchment_id', right_on='gauge_id')
            mem_dfs[rcp][mem] = df.set_index('catchment_id')

    return ens_dfs, mem_dfs


def make_supp_figure(params, out_path, gb, ens_dfs=None, mem_dfs=None, mode='alteration'):
    """mode = 'alteration' (uses ens_dfs) or 'agreement' (uses mem_dfs)."""
    nrows = len(params)
    fig_w = 6.85
    fig_h = nrows * 1.95 + 1.1
    fig, axes = plt.subplots(nrows, 4, figsize=(fig_w, fig_h), dpi=300)
    if nrows == 1:
        axes = axes.reshape(1, -1)

    for r, (ha_col, label, group) in enumerate(params):
        for c, rcp in enumerate(RCPS):
            ax = axes[r, c]
            gb.plot(ax=ax, **MAP_KW)

            if mode == 'alteration':
                df = ens_dfs[rcp]
                colours = df[ha_col].apply(five_cat).map(CAT_COLOURS)
                ax.scatter(df['easting'], df['northing'], c=colours, s=5.0,
                           zorder=2, linewidths=0, rasterized=True)
            else:
                agreement = compute_agreement(mem_dfs[rcp], ha_col, direction='departure')
                base = mem_dfs[rcp]['01'].loc[agreement.index, ['easting', 'northing']]
                colours = agreement.map(AGR_COLOURS)
                ax.scatter(base['easting'], base['northing'], c=colours, s=5.0,
                           zorder=2, linewidths=0, rasterized=True)

            ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect('equal'); ax.axis('off')
            if r == 0:
                ax.set_title(RCP_LABELS[c], fontsize=7.5, fontweight='bold', pad=3)
            if c == 0:
                ax.text(-0.14, 0.5, f'{group}\n{label}', transform=ax.transAxes,
                        fontsize=6.0, ha='right', va='center', multialignment='right',
                        linespacing=1.4)

    bottom_margin = 0.13 if nrows <= 3 else 0.09
    if mode == 'alteration':
        patches = [mpatches.Patch(facecolor=CAT_COLOURS[k], label=l) for k, l in zip(CAT_ORDER, CAT_LABELS)]
        ncol_leg = 2
    else:
        patches = [mpatches.Patch(facecolor=AGR_COLOURS[k], label=l) for k, l in zip(range(5), AGR_LABELS)]
        ncol_leg = 3

    fig.legend(handles=patches, ncol=ncol_leg, loc='lower center', bbox_to_anchor=(0.54, 0.01),
               bbox_transform=fig.transFigure, fontsize=5.8, frameon=True, edgecolor='#cccccc',
               framealpha=0.9, handlelength=1.2, handleheight=1.0, columnspacing=1.0,
               borderpad=0.6, labelspacing=0.4)
    plt.subplots_adjust(left=0.16, right=0.99, top=0.97, bottom=bottom_margin, wspace=0.04, hspace=0.06)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved {out_path}')


def main():
    gb = get_gb_boundary()
    ens_dfs, mem_dfs = load_data()

    # S1-S6: spatial alteration maps
    for i, (name, params) in enumerate(SUPP_GROUPS, start=1):
        make_supp_figure(params, f'{OUTPUTS}figureS{i}_{name}.png', gb,
                          ens_dfs=ens_dfs, mode='alteration')

    # S7-S12: ensemble agreement maps (same parameter groupings)
    for i, (name, params) in enumerate(SUPP_GROUPS, start=7):
        make_supp_figure(params, f'{OUTPUTS}figureS{i}_agreement_{name}.png', gb,
                          mem_dfs=mem_dfs, mode='agreement')


if __name__ == '__main__':
    main()
