"""
figures_3_4_5_main.py
=====================
Figures 3, 4, and 5: IHA/RVA spatial alteration maps, ranked ensemble
summary charts with member symbols, and ensemble agreement maps, using
four representative IHA parameters.

Inputs (edit paths as needed):
    rva_results_{rcp}_ENSEMBLE.csv        for rcp in [rcp26, rcp45, rcp60, rcp85]
    rva_results_{rcp}_{member}.csv        for member in [01, 04, 06, 15]
    camels_gb_v2_topographic_attributes.csv

Outputs:
    figure3_iha_maps.png
    figure4_ensemble_summary.png
    figure5_agreement_maps.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from shared_utils import (five_cat, CAT_COLOURS, CAT_ORDER, CAT_LABELS,
                           AGR_COLOURS, AGR_LABELS, RCPS, RCP_LABELS, MEMBERS,
                           XLIM, YLIM, MAP_KW, PARAM_COLS,
                           get_gb_boundary, load_topo_with_bng, compute_agreement)

UPLOADS = '/mnt/user-data/uploads/'
# chess_scape_iha_analysis.py writes rva_results_*.csv to IHA_Results/ at the
# repository root; this script lives in figures/, hence the one level up.
RVA_DIR = '../IHA_Results/'
OUTPUTS = '/mnt/user-data/outputs/'


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


def make_figure3(ens_dfs, gb):
    """Spatial alteration maps: 4 representative parameters x 4 RCPs."""
    params = [
        ('month_12_Sep_ha', 'September mean flow', 'Group 1'),
        ('min_7day_ha', '7-day minimum flow', 'Group 2'),
        ('doy_max_ha', 'Day of maximum flow', 'Group 3'),
        ('low_pulse_count_ha', 'Low pulse count', 'Group 4'),
    ]
    fig, axes = plt.subplots(4, 4, figsize=(6.85, 8.2), dpi=300)
    for r, (ha_col, label, group) in enumerate(params):
        for c, rcp in enumerate(RCPS):
            ax = axes[r, c]
            df = ens_dfs[rcp]
            gb.plot(ax=ax, **MAP_KW)
            colours = df[ha_col].apply(five_cat).map(CAT_COLOURS)
            ax.scatter(df['easting'], df['northing'], c=colours, s=5.5,
                       zorder=2, linewidths=0, rasterized=True)
            ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect('equal'); ax.axis('off')
            if r == 0:
                ax.set_title(RCP_LABELS[c], fontsize=7.5, fontweight='bold', pad=4)
            if c == 0:
                ax.text(-0.14, 0.5, f'{group}\n{label}', transform=ax.transAxes,
                        fontsize=6.2, ha='right', va='center', multialignment='right',
                        linespacing=1.4)
    patches = [mpatches.Patch(facecolor=CAT_COLOURS[k], label=l)
               for k, l in zip(CAT_ORDER, CAT_LABELS)]
    fig.legend(handles=patches, ncol=2, loc='lower center', bbox_to_anchor=(0.54, 0.00),
               fontsize=6, frameon=True, edgecolor='#cccccc', framealpha=0.9,
               handlelength=1.2, handleheight=1.0, columnspacing=1.0, borderpad=0.6,
               labelspacing=0.4)
    plt.subplots_adjust(left=0.16, right=0.99, top=0.96, bottom=0.13, wspace=0.04, hspace=0.06)
    fig.savefig(OUTPUTS + 'figure3_iha_maps.png', dpi=300, bbox_inches='tight')
    plt.close()
    print('Saved figure3_iha_maps.png')


def compute_fractions(df):
    rows = []
    for ha_col, label in PARAM_COLS:
        cats = df[ha_col].apply(five_cat)
        fracs = {c: (cats == c).mean() for c in CAT_ORDER}
        fracs['label'] = label
        fracs['total_departure'] = fracs['high_departure'] + fracs['moderate_departure']
        rows.append(fracs)
    return pd.DataFrame(rows).sort_values('total_departure', ascending=True)


def departure_fraction(df, ha_col):
    cats = df[ha_col].apply(five_cat)
    return ((cats == 'high_departure') | (cats == 'moderate_departure')).mean()


def make_figure4(ens_dfs, mem_dfs):
    """Ranked ensemble summary bar charts with per-member departure symbols."""
    mem_markers = {'01': 'o', '04': 's', '06': '^', '15': 'D'}
    mem_colour = '#222222'
    label_to_col = {label: col for col, label in PARAM_COLS}

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 9.8), dpi=300)
    axes_flat = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

    for ax, rcp, rcp_label in zip(axes_flat, RCPS, RCP_LABELS):
        df_fracs = compute_fractions(ens_dfs[rcp])
        labels = df_fracs['label'].tolist()
        n = len(labels); y = np.arange(n)

        left = np.zeros(n)
        for cat in CAT_ORDER:
            vals = df_fracs[cat].values
            ax.barh(y, vals, left=left, color=CAT_COLOURS[cat], height=0.75, zorder=2)
            left += vals

        for mem in MEMBERS:
            mem_df = mem_dfs[rcp][mem]
            dep_fracs = [departure_fraction(mem_df, label_to_col[lbl]) for lbl in labels]
            ax.scatter(dep_fracs, y, marker=mem_markers[mem], color=mem_colour, s=18,
                       zorder=4, linewidths=0.3, edgecolors='white')

        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.5)
        ax.set_xlim(0, 1); ax.set_xlabel('Fraction of catchments', fontsize=7)
        ax.set_title(rcp_label, fontsize=8, fontweight='bold', pad=4)
        ax.axvline(0.5, color='#555555', lw=0.7, ls='--', zorder=3)
        ax.tick_params(axis='x', labelsize=6.5); ax.tick_params(axis='y', pad=2)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', lw=0.3, alpha=0.4, color='#cccccc', zorder=0)
        ax.margins(y=0.01)

    cat_patches = [mpatches.Patch(facecolor=CAT_COLOURS[k], label=l)
                   for k, l in zip(CAT_ORDER, CAT_LABELS)]
    mem_labels = ['Member 01', 'Member 04', 'Member 06', 'Member 15']
    mem_handles = [mlines.Line2D([0], [0], marker=mem_markers[m], color='w',
                                  markerfacecolor=mem_colour, markersize=5, label=l,
                                  markeredgecolor='white', markeredgewidth=0.3)
                   for m, l in zip(MEMBERS, mem_labels)]
    fig.legend(handles=cat_patches + mem_handles, ncol=3, loc='lower center',
               bbox_to_anchor=(0.5, 0.00), fontsize=6, frameon=True, edgecolor='#cccccc',
               framealpha=0.9, handlelength=1.2, handleheight=1.0, columnspacing=1.0,
               borderpad=0.6, labelspacing=0.4)
    plt.subplots_adjust(left=0.26, right=0.98, top=0.97, bottom=0.10, wspace=0.42, hspace=0.12)
    fig.savefig(OUTPUTS + 'figure4_ensemble_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print('Saved figure4_ensemble_summary.png')


def make_figure5(mem_dfs, gb):
    """Ensemble agreement maps: 4 representative parameters x 4 RCPs."""
    params = [
        ('month_12_Sep_ha', 'September mean flow', 'Group 1'),
        ('min_7day_ha', '7-day minimum flow', 'Group 2'),
        ('base_flow_index_ha', 'Base flow index', 'Group 2'),
        ('low_pulse_count_ha', 'Low pulse count', 'Group 4'),
    ]
    fig, axes = plt.subplots(4, 4, figsize=(6.85, 8.4), dpi=300)
    for r, (ha_col, label, group) in enumerate(params):
        for c, rcp in enumerate(RCPS):
            ax = axes[r, c]
            agreement = compute_agreement(mem_dfs[rcp], ha_col, direction='departure')
            base = mem_dfs[rcp]['01'].loc[agreement.index, ['easting', 'northing']]
            gb.plot(ax=ax, **MAP_KW)
            colours = agreement.map(AGR_COLOURS)
            ax.scatter(base['easting'], base['northing'], c=colours, s=5.5,
                       zorder=2, linewidths=0, rasterized=True)
            ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect('equal'); ax.axis('off')
            if r == 0:
                ax.set_title(RCP_LABELS[c], fontsize=7.5, fontweight='bold', pad=4)
            if c == 0:
                ax.text(-0.14, 0.5, f'{group}\n{label}', transform=ax.transAxes,
                        fontsize=6.2, ha='right', va='center', multialignment='right',
                        linespacing=1.4)
    patches = [mpatches.Patch(facecolor=AGR_COLOURS[k], label=l)
               for k, l in zip(range(5), AGR_LABELS)]
    fig.legend(handles=patches, ncol=3, loc='lower center', bbox_to_anchor=(0.54, 0.00),
               fontsize=6, frameon=True, edgecolor='#cccccc', framealpha=0.9,
               handlelength=1.2, handleheight=1.0, columnspacing=1.0, borderpad=0.6,
               labelspacing=0.4)
    plt.subplots_adjust(left=0.16, right=0.99, top=0.96, bottom=0.08, wspace=0.04, hspace=0.06)
    fig.savefig(OUTPUTS + 'figure5_agreement_maps.png', dpi=300, bbox_inches='tight')
    plt.close()
    print('Saved figure5_agreement_maps.png')


def main():
    gb = get_gb_boundary()
    ens_dfs, mem_dfs = load_data()
    make_figure3(ens_dfs, gb)
    make_figure4(ens_dfs, mem_dfs)
    make_figure5(mem_dfs, gb)


if __name__ == '__main__':
    main()
