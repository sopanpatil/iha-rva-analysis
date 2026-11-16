"""
supplementary_maps.py
=====================

Supplementary maps built from the per-member RVA statistics:

  figS_direction_<group>   alteration class and direction (median of members)
                           for every IHA parameter, four RCPs, by IHA group
  figS_agreement_<group>   number of members in departure, same layout
  figS_members_<param>     alteration class and direction for each ensemble
                           member (rows) and RCP (columns)
  figS_member_spread       standard deviation of HA across members

Run after rva_member_analysis.py, from figures/:

    python supplementary_maps.py
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

from alteration_figures import (AGR_COLOURS, AGR_LABELS, DIR_CLASSES, HERE, INK, LABEL,
                                MAP_KW, MAP_PARAMS, MEMBERS, RCP_LABELS, RCPS, RVA_DIR, TOPO_CSV,
                                XLIM, YLIM, direction_class, get_gb_boundary, save)

GROUPS = {
    "group1_may_oct": ["month_08_May", "month_09_Jun", "month_10_Jul", "month_11_Aug",
                       "month_12_Sep", "month_01_Oct"],
    "group1_nov_apr": ["month_02_Nov", "month_03_Dec", "month_04_Jan", "month_05_Feb",
                       "month_06_Mar", "month_07_Apr"],
    "group2_minima": ["min_1day", "min_3day", "min_7day", "min_30day", "min_90day",
                      "base_flow_index"],
    "group2_maxima": ["max_1day", "max_3day", "max_7day", "max_30day", "max_90day"],
    "group3_4_timing_pulses": ["doy_min", "doy_max", "low_pulse_count",
                               "low_pulse_duration_mean", "high_pulse_count",
                               "high_pulse_duration_mean"],
    "group5_rates": ["mean_rise_rate", "mean_fall_rate", "number_of_reversals"],
}
DIR_COLOURS = {k: c for k, _, c in DIR_CLASSES}
DIR_HANDLES = [Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="white", label=lab)
               for _, lab, c in DIR_CLASSES]
AGR_HANDLES = [Line2D([], [], marker="o", ls="", ms=5, mfc=AGR_COLOURS[k], mec="white", label=lab)
               for k, lab in enumerate(AGR_LABELS)]


def map_grid(rows, cols, draw, row_label, col_label, handles, name, legend_cols=3):
    gb = get_gb_boundary(os.path.join(HERE, "_gb_boundary.gpkg"),
                         os.path.join(HERE, "_ne_10m_countries.geojson"))
    fig, axes = plt.subplots(len(rows), len(cols), figsize=(7.2, 2.4 * len(rows)), squeeze=False)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            ax = axes[i, j]
            gb.plot(ax=ax, **MAP_KW)
            draw(ax, r, c)
            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.set_aspect("equal")
            ax.axis("off")
            if i == 0:
                ax.set_title(col_label(c), fontsize=9, color=INK)
            if j == 0:
                ax.text(-0.02, 0.5, row_label(r), transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.5, color=INK)
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=legend_cols, fontsize=7,
                   bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(wspace=0.02, hspace=0.04, bottom=0.25 / len(rows))
    return fig


def scatter_classes(ax, s, colour_col):
    s = s.assign(z=(s.cls != "low").astype(int)).sort_values("z")
    ax.scatter(s.x, s.y, c=s[colour_col], s=4.5, linewidths=0.15, edgecolors="white", zorder=2)


def main():
    topo = pd.read_csv(TOPO_CSV).set_index("gauge_id")
    coords = topo[["gauge_easting", "gauge_northing"]].rename(
        columns={"gauge_easting": "x", "gauge_northing": "y"})
    ens = pd.read_csv(os.path.join(RVA_DIR, "rva_catchment_ensemble.csv.gz"))
    ens["cls"] = ens.apply(direction_class, axis=1)
    ens["colour"] = ens.cls.map(DIR_COLOURS)
    ens = ens.merge(coords, left_on="catchment_id", right_index=True)
    idx = ens.set_index(["parameter", "rcp"]).sort_index()

    for group, params in GROUPS.items():
        fig = map_grid(params, RCPS, lambda ax, p, rcp: scatter_classes(ax, idx.loc[(p, rcp)], "colour"),
                       lambda p: LABEL[p], lambda r: RCP_LABELS[RCPS.index(r)], DIR_HANDLES,
                       f"figS_direction_{group}")
        save(fig, f"figS_direction_{group}")

        def draw_agr(ax, p, rcp):
            s = idx.loc[(p, rcp)].sort_values("n_departure")
            ax.scatter(s.x, s.y, c=s.n_departure.map(AGR_COLOURS), s=4.5, linewidths=0.15,
                       edgecolors="white", zorder=2)
        fig = map_grid(params, RCPS, draw_agr, lambda p: LABEL[p],
                       lambda r: RCP_LABELS[RCPS.index(r)], AGR_HANDLES, f"figS_agreement_{group}")
        save(fig, f"figS_agreement_{group}")

    member = pd.read_csv(os.path.join(RVA_DIR, "rva_catchment_member.csv.gz"), dtype={"member": str})
    member = member[member.parameter.isin(MAP_PARAMS)].merge(coords, left_on="catchment_id",
                                                              right_index=True)
    member["ha_median"] = member.ha
    member["frac_below_median"] = member.frac_below
    member["frac_above_median"] = member.frac_above
    member["cls"] = member.apply(direction_class, axis=1)
    member["colour"] = member.cls.map(DIR_COLOURS)
    midx = member.set_index(["parameter", "member", "rcp"]).sort_index()
    for p in MAP_PARAMS:
        fig = map_grid(MEMBERS, RCPS, lambda ax, m, rcp: scatter_classes(ax, midx.loc[(p, m, rcp)], "colour"),
                       lambda m: f"Member {m}", lambda r: RCP_LABELS[RCPS.index(r)], DIR_HANDLES,
                       f"figS_members_{p}")
        fig.suptitle(LABEL[p], fontsize=10, color=INK, y=0.93)
        save(fig, f"figS_members_{p}")

    bounds = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 2]
    cmap = ListedColormap(["#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"])
    norm = BoundaryNorm(bounds, cmap.N)

    def draw_sd(ax, p, rcp):
        s = idx.loc[(p, rcp)].sort_values("ha_sd")
        ax.scatter(s.x, s.y, c=s.ha_sd, cmap=cmap, norm=norm, s=4.5, linewidths=0.15,
                   edgecolors="white", zorder=2)
    labels = ["< 0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5", "≥ 0.5"]
    handles = [Line2D([], [], marker="o", ls="", ms=5, mfc=cmap(k), mec="white", label=lab)
               for k, lab in enumerate(labels)]
    fig = map_grid(MAP_PARAMS, RCPS, draw_sd, lambda p: LABEL[p],
                   lambda r: RCP_LABELS[RCPS.index(r)], handles, "figS_member_spread", legend_cols=6)
    fig.text(0.5, 0.035, "Standard deviation of HA across the four members", ha="center",
             fontsize=8, color=INK)
    save(fig, "figS_member_spread")


if __name__ == "__main__":
    main()
