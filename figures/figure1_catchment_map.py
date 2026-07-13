"""
figure1_catchment_map.py
========================
Figure 1: Location and validation performance of the 671 CAMELS-GB v2 study
catchments, coloured by validation KGE, with excluded catchments shown in
grey and GB country/regional boundaries for geographic context.

KGE source
----------
Reads validation KGE directly from hbv-model's calibrated_parameters.csv
via catchment_filter.py (the same source and the same KGE >= 0.5 threshold
used by chess_scape_iha_analysis.py), rather than a separately maintained
file. This ensures the retained/excluded catchment counts shown here always
match the 621/671 catchments actually used in the IHA/RVA analysis.

Inputs (edit paths as needed):
    ../../hbv-model/calibrated_parameters.csv   (via catchment_filter.py, one
                                                  level up again since this
                                                  script lives in figures/)
    camels_gb_v2_topographic_attributes.csv

Output:
    figure1_catchment_map.png  (300 dpi)
"""

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
import pandas as pd

from shared_utils import get_gb_boundary, get_gb_regions, load_topo_with_bng, XLIM, YLIM

# Make the repository root (one level up from figures/) importable so we can
# reuse catchment_filter.py's KGE source path and threshold directly.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from catchment_filter import DEFAULT_THRESHOLD

# --- Paths (edit to match your local file locations) ---
# Default assumes hbv-model is a sibling of the iha-rva-analysis repository
# root, i.e. this script is at iha-rva-analysis/figures/figure1_catchment_map.py.
KGE_PATH  = os.path.join(os.path.dirname(__file__), '..', '..', 'hbv-model',
                          'calibrated_parameters.csv')
TOPO_PATH = '/mnt/user-data/uploads/camels_gb_v2_topographic_attributes.csv'
OUT_PATH  = '/mnt/user-data/outputs/figure1_catchment_map.png'


def main():
    kge  = pd.read_csv(KGE_PATH)
    topo = load_topo_with_bng(TOPO_PATH)
    df = topo.merge(kge, left_on='gauge_id', right_on='gauge_id')

    retained = df[df['validation_kge'] >= DEFAULT_THRESHOLD].copy()
    excluded = df[(df['validation_kge'] < DEFAULT_THRESHOLD) | df['validation_kge'].isna()].copy()

    gb = get_gb_boundary()
    regions = get_gb_regions()

    # Warm-neutral structural region fill; country boundaries thicker than
    # English regional boundaries. KGE colour scale (YlGnBu) is kept visually
    # independent of the region fill to avoid confusion between the two.
    countries = regions.copy()
    countries['country'] = countries['display_region'].apply(
        lambda x: 'Scotland' if x == 'Scotland'
        else 'Wales' if x == 'Wales' else 'England')
    country_bounds = countries.dissolve(by='country').reset_index()

    cmap = plt.cm.YlGnBu
    norm = mcolors.Normalize(vmin=DEFAULT_THRESHOLD, vmax=1.0)

    fig, ax = plt.subplots(figsize=(6.85, 8.2), dpi=300)

    regions.plot(ax=ax, color='#f2f2f2', edgecolor='#bbbbbb', linewidth=0.35, zorder=1)
    country_bounds.plot(ax=ax, color='none', edgecolor='#777777', linewidth=1.0, zorder=2)
    gb.plot(ax=ax, color='none', edgecolor='#555555', linewidth=0.7, zorder=3)

    ax.scatter(excluded['easting'], excluded['northing'],
               c='#cccccc', edgecolors='#999999', linewidths=0.3, s=12, zorder=4)
    ax.scatter(retained['easting'], retained['northing'],
               c=retained['validation_kge'], cmap=cmap, norm=norm,
               edgecolors='none', s=12, zorder=5)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect('equal')
    ax.axis('off')

    # Region labels: bold italic for countries, lighter italic for English regions
    label_pos = {
        'Scotland': (250000, 750000, 7.5, 'bold'),
        'Wales': (260000, 240000, 7.5, 'bold'),
        'North East\nEngland': (430000, 580000, 5.2, 'normal'),
        'North West\nEngland': (330000, 510000, 5.2, 'normal'),
        'Yorkshire and\nthe Humber': (470000, 490000, 5.2, 'normal'),
        'East\nMidlands': (490000, 375000, 5.2, 'normal'),
        'West\nMidlands': (385000, 330000, 5.2, 'normal'),
        'East Anglia': (565000, 295000, 5.2, 'normal'),
        'South East\nEngland': (510000, 175000, 5.2, 'normal'),
        'South West\nEngland': (215000, 110000, 5.2, 'normal'),
    }
    for label, (ex, ny, fs, fw) in label_pos.items():
        ax.text(ex, ny, label, ha='center', va='center', fontsize=fs, color='#555555',
                fontweight=fw, fontstyle='italic', multialignment='center', zorder=6)

    # Scale bar (200 km)
    bar_x0, bar_y, bar_len = -160000, 30000, 200000
    ax.plot([bar_x0, bar_x0+bar_len], [bar_y, bar_y], color='black', lw=1.5,
            solid_capstyle='butt')
    for x in [bar_x0, bar_x0+bar_len]:
        ax.plot([x, x], [bar_y-8000, bar_y+8000], color='black', lw=1.2)
    ax.text(bar_x0+bar_len/2, bar_y+18000, '200 km', ha='center', va='bottom', fontsize=6.5)

    ax.text(-160000, 130000, f'Grey dots: excluded catchments (n = {len(excluded)})',
            fontsize=6, ha='left', va='bottom', color='#666666')

    cbar_ax = fig.add_axes([0.15, 0.04, 0.55, 0.016])
    cb = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')
    cb.set_label('Validation KGE', fontsize=7, labelpad=3)
    cb.ax.tick_params(labelsize=6.5)
    cb.set_ticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    plt.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.07)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved {OUT_PATH}  (retained={len(retained)}, excluded={len(excluded)})')


if __name__ == '__main__':
    main()
