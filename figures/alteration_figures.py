"""
alteration_figures.py
=====================

Figures built from the per-member RVA statistics (../rva_member_analysis.py):

  fig_parameter_summary   percentage of catchments departing from the baseline
                          range per parameter and RCP (mean of the four members,
                          split by direction), with the member range and the
                          no-change level
  fig_excess_departure    departure minus the no-change level, parameters x RCPs
  fig_effect_size         Cliff's delta across catchments for RCP2.6 and RCP8.5,
                          against the no-change 5-95% range
  fig_direction_maps      maps of the median-across-members HA class, coloured by
                          direction of departure, for four parameters x four RCPs
  fig_agreement_maps      number of members in departure, same layout
  fig_attribute_scatter   Cliff's delta (RCP8.5) against aridity, snow-influenced
                          catchments highlighted
  fig_ensemble_mean_vs_members
                          departure from ensemble-mean discharge against the
                          mean of member-wise departure (supplementary)

Inputs: ../IHA_Results/member_rva/* (../rva_member_analysis.py, ../ensemble_mean_rva.py),
../IHA_Results/attributes/ (attribute_associations.py), and the CAMELS-GB v2
topographic and climatic attribute files in ATTR_DIR (default: figures/).
Outputs: PNG files (300 dpi) in FIG_OUT (default: figures/output/); with FIG_PDF=1,
also vector PDF copies.

    python alteration_figures.py
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shared_utils import (AGR_COLOURS, AGR_LABELS, MAP_KW, PARAM_COLS, RCP_LABELS,
                          XLIM, YLIM, get_gb_boundary)

HERE = os.path.dirname(os.path.abspath(__file__))
RVA_DIR = os.path.join(HERE, "..", "IHA_Results", "member_rva")
ATTR_DIR = os.environ.get("ATTR_DIR", HERE)
TOPO_CSV = os.path.join(ATTR_DIR, "camels_gb_v2_topographic_attributes.csv")
OUT_DIR = os.environ.get("FIG_OUT", os.path.join(HERE, "output"))

RCPS = [26, 45, 60, 85]
MEMBERS = ["01", "04", "06", "15"]
LABEL = {c[:-3]: lab for c, lab in PARAM_COLS}
MAP_PARAMS = ["month_12_Sep", "min_7day", "base_flow_index", "max_1day"]

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BELOW, ABOVE = "#2a78d6", "#eb6834"
MEMBER_STYLE = {"01": ("#2a78d6", "o"), "04": ("#eb6834", "s"),
                "06": ("#1baf7a", "^"), "15": ("#eda100", "D")}
DIR_CLASSES = [
    ("high_below", "High departure, lower", "#08519c"),
    ("mod_below", "Moderate departure, lower", "#6baed6"),
    ("low", "Low alteration", "#d9d9d9"),
    ("concentration", "Concentration", "#9e9ac8"),
    ("mod_above", "Moderate departure, higher", "#fdae6b"),
    ("high_above", "High departure, higher", "#d94801"),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
})


def save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, f"{name}.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    if os.environ.get("FIG_PDF"):  # vector copy
        with plt.rc_context({"pdf.fonttype": 42}):
            fig.savefig(path[:-4] + ".pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {path}")


def load():
    summary = pd.read_csv(os.path.join(RVA_DIR, "rva_summary_by_parameter.csv"))
    member = pd.read_csv(os.path.join(RVA_DIR, "rva_catchment_member.csv.gz"),
                         dtype={"member": str})
    ens = pd.read_csv(os.path.join(RVA_DIR, "rva_catchment_ensemble.csv.gz"))
    null = pd.read_csv(os.path.join(RVA_DIR, "rva_null_catchment.csv.gz"))
    return summary, member, ens, null


def parameter_order(summary):
    s = summary[summary.rcp == 85].sort_values("pct_departure", ascending=False)
    return list(s.parameter)


# ---------------------------------------------------------------------------
def fig_parameter_summary(summary, member):
    order = parameter_order(summary)
    dep = member.category.isin(["mod_departure", "high_departure"])
    per_member = (member.assign(dep=dep).groupby(["rcp", "parameter", "member"]).dep.mean() * 100)

    fig, axes = plt.subplots(1, 4, figsize=(7.2, 6.2), sharey=True)
    y = np.arange(len(order))
    for ax, rcp, lab in zip(axes, RCPS, RCP_LABELS):
        s = summary[summary.rcp == rcp].set_index("parameter").loc[order]
        ax.barh(y, s.pct_departure_below, height=0.7, color=BELOW, edgecolor="white", linewidth=0.5)
        ax.barh(y, s.pct_departure_above, left=s.pct_departure_below, height=0.7, color=ABOVE,
                edgecolor="white", linewidth=0.5)
        for m in MEMBERS:
            mk = MEMBER_STYLE[m][1]
            v = per_member.loc[rcp].xs(m, level="member").loc[order]
            ax.plot(v, y, mk, ms=2.8, mfc="white", mec=INK, mew=0.5, zorder=3)
        ax.plot(s.pct_departure_nochange, y, "|", ms=8, mew=1.6, color="#7a2a8c", zorder=4)
        ax.set_title(lab, fontsize=9, color=INK)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.grid(axis="x", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([LABEL[p] for p in order])
    axes[0].invert_yaxis()
    fig.supxlabel("Catchments in departure (%); bars = mean of the four members", fontsize=8,
                  y=0.065)
    handles = [Patch(color=BELOW, label="Departure, shifted lower"),
               Patch(color=ABOVE, label="Departure, shifted higher"),
               Line2D([], [], color="#7a2a8c", marker="|", ls="", ms=8, mew=1.6,
                      label="No-change level")]
    handles += [Line2D([], [], marker=MEMBER_STYLE[m][1], ls="", ms=4, mfc="white",
                       mec=INK, mew=0.6, label=f"Member {m}") for m in MEMBERS]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, -0.02), fontsize=7)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, "fig_parameter_summary")


def fig_excess_departure(summary):
    order = parameter_order(summary)
    w = summary.pivot(index="parameter", columns="rcp", values="pct_departure").loc[order]
    n = summary.drop_duplicates("parameter").set_index("parameter").pct_departure_nochange.loc[order]
    ex = w.sub(n, axis=0)
    fig, ax = plt.subplots(figsize=(3.6, 6.2))
    im = ax.imshow(ex.values, cmap="Blues", vmin=0, vmax=80, aspect="auto")
    for i in range(ex.shape[0]):
        for j in range(ex.shape[1]):
            v = ex.values[i, j]
            ax.text(j, i, f"{int(round(v))}", ha="center", va="center", fontsize=6.5,
                    color="white" if v > 45 else INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(RCP_LABELS)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([LABEL[p] for p in order])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("Departure above the no-change level\n(percentage points)")
    cb.outline.set_visible(False)
    save(fig, "fig_excess_departure")


def fig_effect_size(summary, ens, null):
    order = parameter_order(summary)
    lo = null.groupby("parameter").cliffs_delta.quantile(0.05)
    hi = null.groupby("parameter").cliffs_delta.quantile(0.95)
    fig, ax = plt.subplots(figsize=(3.8, 6.2))
    y = np.arange(len(order))
    ax.barh(y, (hi - lo).loc[order], left=lo.loc[order], height=0.8, color="#ecebe7",
            label="No-change 5–95% range")
    for k, (rcp, col, lab) in enumerate([(26, "#9ec5f0", "RCP2.6"), (85, "#1c5aa6", "RCP8.5")]):
        d = ens[ens.rcp == rcp]
        data = [d.loc[d.parameter == p, "cliffs_delta_median"].dropna().values for p in order]
        ax.boxplot(data, positions=y + (k - 0.5) * 0.32, orientation="horizontal", widths=0.26,
                   patch_artist=True, showfliers=False, whis=(5, 95),
                   boxprops=dict(facecolor=col, edgecolor=col, linewidth=0.6),
                   medianprops=dict(color="white", linewidth=1),
                   whiskerprops=dict(color=col, linewidth=0.8), capprops=dict(linewidth=0))
        ax.plot([], [], "s", color=col, label=lab)
    ax.axvline(0, color=INK2, linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[p] for p in order])
    ax.invert_yaxis()
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Cliff's delta, future vs baseline\n(median of members; box 25–75%, whiskers 5–95%)")
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", bbox_to_anchor=(0, -0.2), ncol=3, fontsize=7)
    save(fig, "fig_effect_size")


def direction_class(row):
    h = row.ha_median
    if h > 0.33:
        return "concentration"
    if h >= -0.33:
        return "low"
    side = "below" if row.frac_below_median >= row.frac_above_median else "above"
    return f"{'high' if h < -0.67 else 'mod'}_{side}"


def _map_grid(draw, title_fn, legend_handles, name):
    gb = get_gb_boundary(os.path.join(HERE, "_gb_boundary.gpkg"),
                         os.path.join(HERE, "_ne_10m_countries.geojson"))
    fig, axes = plt.subplots(len(MAP_PARAMS), 4, figsize=(7.2, 9.6))
    for i, p in enumerate(MAP_PARAMS):
        for j, rcp in enumerate(RCPS):
            ax = axes[i, j]
            gb.plot(ax=ax, **MAP_KW)
            draw(ax, p, rcp)
            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.set_aspect("equal")
            ax.axis("off")
            if i == 0:
                ax.set_title(RCP_LABELS[j], fontsize=9, color=INK)
            if j == 0:
                ax.text(-0.02, 0.5, title_fn(p), transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.5, color=INK)
    fig.legend(handles=legend_handles, loc="lower center", ncol=3, fontsize=7,
               bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(wspace=0.02, hspace=0.04, bottom=0.06)
    save(fig, name)


def fig_direction_maps(ens, coords):
    d = ens[ens.parameter.isin(MAP_PARAMS)].copy()
    d["cls"] = d.apply(direction_class, axis=1)
    colours = {k: c for k, _, c in DIR_CLASSES}
    d = d.merge(coords, left_on="catchment_id", right_index=True)

    def draw(ax, p, rcp):
        s = d[(d.parameter == p) & (d.rcp == rcp)].copy()
        s["z"] = s.cls.map(lambda c: 0 if c == "low" else 1)
        s = s.sort_values("z")
        ax.scatter(s.x, s.y, c=s.cls.map(colours), s=4.5, linewidths=0.15, edgecolors="white", zorder=2)

    handles = [Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="white", label=lab)
               for _, lab, c in DIR_CLASSES]
    _map_grid(draw, lambda p: LABEL[p], handles, "fig_direction_maps")


def fig_agreement_maps(ens, coords):
    d = ens[ens.parameter.isin(MAP_PARAMS)].merge(coords, left_on="catchment_id", right_index=True)

    def draw(ax, p, rcp):
        s = d[(d.parameter == p) & (d.rcp == rcp)].sort_values("n_departure")
        ax.scatter(s.x, s.y, c=s.n_departure.map(AGR_COLOURS), s=4.5, linewidths=0.15,
                   edgecolors="white", zorder=2)

    handles = [Line2D([], [], marker="o", ls="", ms=5, mfc=AGR_COLOURS[k], mec="white", label=lab)
               for k, lab in enumerate(AGR_LABELS)]
    _map_grid(draw, lambda p: LABEL[p], handles, "fig_agreement_maps")


def fig_ensemble_mean_vs_members(summary):
    ens_mean = pd.read_csv(os.path.join(RVA_DIR, "rva_summary_ensemble_mean.csv"))
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    for rcp, lab, col in zip(RCPS, RCP_LABELS, ["#9ec5f0", "#5a9be0", "#2a78d6", "#0d3d7a"]):
        v1 = ens_mean[ens_mean.rcp == rcp].set_index("parameter").pct_departure
        s = summary[summary.rcp == rcp].set_index("parameter").pct_departure
        ax.plot(s.loc[v1.index], v1, "o", ms=3.5, mfc=col, mec="white", mew=0.3, label=lab)
    ax.plot([0, 100], [0, 100], color=INK2, lw=0.6, ls="--")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal")
    ax.set_xlabel("Departure, mean of member-wise analyses (%)")
    ax.set_ylabel("Departure, ensemble-mean discharge (%)")
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "fig_ensemble_mean_vs_members")


def fig_attribute_scatter(ens, attrs_dir):
    """Cliff's delta (RCP8.5) against aridity, snow-influenced catchments highlighted."""
    files = ["camels_gb_v2_climatic_attributes.csv"]
    clim = pd.concat([pd.read_csv(os.path.join(attrs_dir, f)).set_index("gauge_id") for f in files], axis=1)
    rho = pd.read_csv(os.path.join(HERE, "..", "IHA_Results", "attributes",
                                   "attribute_rank_correlations.csv"))
    params = ["min_7day", "month_12_Sep", "max_1day", "month_06_Mar"]
    snow_col, other_col = "#6250d6", "#b4b2ab"
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.3), sharey=True)
    for ax, p in zip(axes, params):
        d = ens[(ens.rcp == 85) & (ens.parameter == p)].set_index("catchment_id")
        d = d.join(clim[["aridity", "frac_snow"]], how="inner")
        snow = d.frac_snow > 0.05
        ax.scatter(d.aridity[~snow], d.cliffs_delta_median[~snow], s=5, color=other_col,
                   linewidths=0, alpha=0.8)
        ax.scatter(d.aridity[snow], d.cliffs_delta_median[snow], s=9, color=snow_col,
                   edgecolors="white", linewidths=0.3)
        r = rho[(rho.rcp == 85) & (rho.parameter == p)].set_index("attribute")
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(LABEL[p], fontsize=8.5, color=INK, pad=14)
        ax.text(0.5, 1.02, f"ρ aridity {r.loc['aridity', 'rho']:.2f}; "
                           f"partial ρ snow {r.loc['frac_snow', 'partial_rho']:.2f}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=6.5, color=INK2)
        ax.set_xlabel("Aridity (PET/P)")
        ax.set_ylim(-1, 1)
        ax.grid(color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Cliff's delta, RCP8.5\n(median of members)")
    handles = [Line2D([], [], marker="o", ls="", ms=4, color=other_col, label="Snow fraction ≤ 0.05"),
               Line2D([], [], marker="o", ls="", ms=4.5, mfc=snow_col, mec="white",
                      label="Snow fraction > 0.05")]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.08), fontsize=7)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "fig_attribute_scatter")


def main():
    summary, member, ens, null = load()
    topo = pd.read_csv(TOPO_CSV).set_index("gauge_id")
    coords = topo[["gauge_easting", "gauge_northing"]].rename(
        columns={"gauge_easting": "x", "gauge_northing": "y"})
    fig_parameter_summary(summary, member)
    fig_excess_departure(summary)
    fig_effect_size(summary, ens, null)
    fig_direction_maps(ens, coords)
    fig_agreement_maps(ens, coords)
    fig_ensemble_mean_vs_members(summary)
    fig_attribute_scatter(ens, ATTR_DIR)


if __name__ == "__main__":
    main()
