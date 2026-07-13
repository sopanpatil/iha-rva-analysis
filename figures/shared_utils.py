"""
shared_utils.py
================
Shared functions and constants used across all figure-generation scripts
for the CHESS-SCAPE IHA/RVA manuscript (Patil et al.).

Requires: pandas, numpy, geopandas, matplotlib, pyproj, shapely

Classification logic
---------------------
The five-category HA classification itself is NOT reimplemented here.
It is imported directly from ../iha_categories.py, which is the single
shared source of truth used by chess_scape_iha_analysis.py to produce
the rva_results_*.csv files these figures read, ensuring the figures
always use exactly the same thresholds as the underlying results. This
module translates iha_categories.py's short keys ("mod_departure", etc.)
into the more verbose keys used throughout the figure scripts
("moderate_departure", etc.) rather than re-deriving the thresholds.
"""

import os
import sys
import pandas as pd
import numpy as np
import geopandas as gpd
from pyproj import Transformer

# Make the repository root (one level up from figures/) importable so we
# can pull in the canonical classification module.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from iha_categories import classify_ha as _classify_ha_short

# ---------------------------------------------------------------------------
# Five-category HA classification (single source of truth: ../iha_categories.py)
# ---------------------------------------------------------------------------
# Translation from iha_categories.py's short keys to the verbose keys used
# throughout the figure scripts. Colours and boundary thresholds are NOT
# redefined here -- only the display keys/labels/colours are local to plotting.
_SHORT_TO_VERBOSE = {
    'low':                'low',
    'mod_departure':      'moderate_departure',
    'high_departure':     'high_departure',
    'mod_concentration':  'moderate_concentration',
    'high_concentration': 'high_concentration',
}

CAT_COLOURS = {
    'high_departure':          '#08519c',
    'moderate_departure':      '#6baed6',
    'low':                     '#cccccc',
    'moderate_concentration':  '#fd8d3c',
    'high_concentration':      '#a50f15',
}
CAT_ORDER = ['high_departure', 'moderate_departure', 'low',
             'moderate_concentration', 'high_concentration']
CAT_LABELS = ['High departure (HA < \u22120.67)',
              'Moderate departure (\u22120.67 \u2264 HA < \u22120.33)',
              'Low alteration (|HA| < 0.33)',
              'Moderate concentration (0.33 < HA \u2264 0.67)',
              'High concentration (HA > 0.67)']

# Ensemble agreement colour scale (0-4 members agreeing on departure)
AGR_COLOURS = {0: '#d9d9d9', 1: '#c6dbef', 2: '#6baed6', 3: '#2171b5', 4: '#08306b'}
AGR_LABELS = ['0 members (no departure)', '1 member', '2 members',
              '3 members', '4 members (full agreement)']

RCPS       = ['rcp26', 'rcp45', 'rcp60', 'rcp85']
RCP_LABELS = ['RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
MEMBERS    = ['01', '04', '06', '15']

# BNG map extent used consistently across all spatial figures
XLIM = (-200000, 700000)
YLIM = (0, 1250000)
MAP_KW = dict(color='#f5f5f5', edgecolor='#bbbbbb', linewidth=0.3, zorder=1)

# All 32 IHA parameters: (ha_column_name, display_label). These names match
# the columns actually written by chess_scape_iha_analysis.py's
# compute_rva_catchment() -- verified against that script directly.
PARAM_COLS = [
    ('month_01_Oct_ha', 'October'), ('month_02_Nov_ha', 'November'),
    ('month_03_Dec_ha', 'December'), ('month_04_Jan_ha', 'January'),
    ('month_05_Feb_ha', 'February'), ('month_06_Mar_ha', 'March'),
    ('month_07_Apr_ha', 'April'), ('month_08_May_ha', 'May'),
    ('month_09_Jun_ha', 'June'), ('month_10_Jul_ha', 'July'),
    ('month_11_Aug_ha', 'August'), ('month_12_Sep_ha', 'September'),
    ('min_1day_ha', '1-day minimum'), ('min_3day_ha', '3-day minimum'),
    ('min_7day_ha', '7-day minimum'), ('min_30day_ha', '30-day minimum'),
    ('min_90day_ha', '90-day minimum'), ('max_1day_ha', '1-day maximum'),
    ('max_3day_ha', '3-day maximum'), ('max_7day_ha', '7-day maximum'),
    ('max_30day_ha', '30-day maximum'), ('max_90day_ha', '90-day maximum'),
    ('base_flow_index_ha', 'Base flow index'), ('doy_min_ha', 'Day of minimum flow'),
    ('doy_max_ha', 'Day of maximum flow'), ('low_pulse_count_ha', 'Low pulse count'),
    ('low_pulse_duration_mean_ha', 'Low pulse duration'), ('high_pulse_count_ha', 'High pulse count'),
    ('high_pulse_duration_mean_ha', 'High pulse duration'), ('mean_rise_rate_ha', 'Mean rise rate'),
    ('mean_fall_rate_ha', 'Mean fall rate'), ('number_of_reversals_ha', 'Number of reversals'),
]


def five_cat(ha):
    """
    Classify a single HA value into one of the five alteration categories,
    using the canonical thresholds in ../iha_categories.py.

    NaN handling: a catchment/parameter with insufficient baseline or future
    data (rva_results_*.csv marks these "insufficient data", HA = NaN) is
    mapped to 'low' here purely so it renders as neutral grey on the spatial
    maps rather than raising an error or leaving a blank pixel. This is a
    plotting convenience only -- it does NOT mean these catchments are
    numerically classified as "low alteration" in the underlying results
    tables, and any reported departure/concentration fractions computed
    from these categories should be checked against the actual count of
    "insufficient data" catchments for the parameter in question (see
    rva_summary_by_parameter_*.csv), since silently folding them into
    "low" will understate both departure and concentration fractions if
    there are more than a handful of them for a given parameter.
    """
    if pd.isna(ha):
        return 'low'
    return _SHORT_TO_VERBOSE[_classify_ha_short(ha)]


def get_gb_boundary(cache_path='_gb_boundary.gpkg',
                     countries_geojson='_ne_10m_countries.geojson'):
    """
    Return a GeoDataFrame of the Great Britain coastline in British National
    Grid (EPSG:27700), with Northern Ireland and adjacent small islands
    removed. Downloads Natural Earth admin-0 countries data if not cached.
    """
    if os.path.exists(cache_path):
        return gpd.read_file(cache_path)

    if not os.path.exists(countries_geojson):
        import urllib.request
        url = ('https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
               'master/geojson/ne_10m_admin_0_countries.geojson')
        urllib.request.urlretrieve(url, countries_geojson)

    world = gpd.read_file(countries_geojson)
    uk = world[world['NAME'] == 'United Kingdom'].to_crs('EPSG:27700')
    polys = list(uk.geometry.iloc[0].geoms)

    transformer = Transformer.from_crs('EPSG:27700', 'EPSG:4326', always_xy=True)
    gb_polys = []
    for poly in polys:
        cx, cy = poly.centroid.x, poly.centroid.y
        lon, lat = transformer.transform(cx, cy)
        # Northern Ireland occupies roughly lon < -5.4, 53.8 < lat < 55.4
        if lon < -5.4 and 53.8 < lat < 55.4:
            continue
        gb_polys.append(poly)

    from shapely.geometry import MultiPolygon
    gb_gdf = gpd.GeoDataFrame(geometry=[MultiPolygon(gb_polys)], crs='EPSG:27700')
    gb_gdf.to_file(cache_path, driver='GPKG')
    return gb_gdf


def get_gb_regions(cache_path='_gb_regions.gpkg',
                    admin1_geojson='_ne_10m_admin1.geojson'):
    """
    Return a GeoDataFrame of GB regions in British National Grid (EPSG:27700):
    Scotland, Wales, and the nine official English regions (Natural Earth's
    'East' region relabelled 'East Anglia' to match common hydrological
    usage). Downloads Natural Earth admin-1 states/provinces data if not
    cached.

    Note: verified directly against the live Natural Earth admin-1 dataset
    (raw 'region' field) that the region_map keys below match exactly --
    all English/Scottish region names used as dict keys (e.g. 'East',
    'Eastern', 'South Western') are literal Natural Earth field values,
    not 'East of England' or similar longer forms.
    """
    if os.path.exists(cache_path):
        return gpd.read_file(cache_path)

    if not os.path.exists(admin1_geojson):
        import urllib.request
        url = ('https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
               'master/geojson/ne_10m_admin_1_states_provinces.geojson')
        urllib.request.urlretrieve(url, admin1_geojson)

    admin1 = gpd.read_file(admin1_geojson)
    gb_all = admin1[admin1['iso_a2'] == 'GB'].copy()
    gb_all = gb_all[~gb_all['region'].isin(['Northern Ireland'])]

    region_map = {
        'Eastern': 'Scotland', 'Highlands and Islands': 'Scotland',
        'North Eastern': 'Scotland', 'South Western': 'Scotland',
        'East Wales': 'Wales', 'West Wales and the Valleys': 'Wales',
        'North East': 'North East England', 'North West': 'North West England',
        'Yorkshire and the Humber': 'Yorkshire and\nthe Humber',
        'East Midlands': 'East Midlands', 'West Midlands': 'West Midlands',
        'East': 'East Anglia', 'South East': 'South East England',
        'Greater London': 'South East England', 'South West': 'South West England',
    }
    gb_all['display_region'] = gb_all['region'].map(region_map)
    regions = gb_all.dissolve(by='display_region').reset_index().to_crs('EPSG:27700')
    regions.to_file(cache_path, driver='GPKG')
    return regions


def load_topo_with_bng(topo_path):
    """Load CAMELS-GB v2 topographic attributes and add BNG easting/northing columns."""
    topo = pd.read_csv(topo_path)
    transformer = Transformer.from_crs('EPSG:4326', 'EPSG:27700', always_xy=True)
    topo['easting'], topo['northing'] = transformer.transform(
        topo['gauge_lon'].values, topo['gauge_lat'].values)
    return topo


def compute_agreement(mem_dfs_dict, ha_col, direction='departure'):
    """
    Compute the number of ensemble members (out of 4) agreeing on departure
    (HA < -0.33) or concentration (HA > 0.33) for a given parameter.

    mem_dfs_dict : dict of {member_label: DataFrame}, each indexed by catchment_id,
                   already merged with easting/northing columns.
    Returns a Series of agreement counts (0-4), aligned to the first member's index.
    """
    members = list(mem_dfs_dict.keys())
    base_index = mem_dfs_dict[members[0]].index
    agreement = np.zeros(len(base_index), dtype=int)
    for mem in members:
        ha = mem_dfs_dict[mem][ha_col].reindex(base_index)
        if direction == 'departure':
            mask = ha < -0.33
        else:
            mask = ha > 0.33
        agreement += mask.fillna(False).astype(int).values
    return pd.Series(agreement, index=base_index)
