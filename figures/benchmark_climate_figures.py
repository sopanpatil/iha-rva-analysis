"""
benchmark_climate_figures.py
============================

Figures and tables from the benchmark and period-mean outputs:

  fig_benchmark          (a) cumulative distributions of KGE and log-KGE for the
                         calibration and validation periods; (b) Q95 bias and
                         low-flow volume bias; (c) percentage of catchments with
                         low alteration per IHA parameter for each step of the
                         modelling chain, against the no-change level
  fig_member_climate     (a, b) change in temperature against change in summer and
                         winter precipitation by member and RCP; (c) low-flow and
                         high-flow departure against warming
  fig_change_maps        RCP8.5 change (median of members) in summer precipitation,
                         annual PET, annual AET and summer soil moisture
  figS_change_maps_all   as fig_change_maps plus temperature, for all four RCPs

Tables written as CSV files to FIG_OUT/tables/:
  table1_members.csv     changes in temperature and precipitation per member
                         (RCP2.6 / RCP8.5), averaged over the catchments
  benchmark_lowflow.csv  HBV performance with observed forcing (Table S1)
  benchmark_rva.csv      benchmark of the modelling chain (Table S2)

Inputs: ../IHA_Results/benchmark/, ../IHA_Results/period_means/,
../IHA_Results/member_rva/ and the CAMELS-GB v2 topographic attributes file
in ATTR_DIR (default: figures/).

    python benchmark_climate_figures.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

from alteration_figures import (GRID, HERE, INK, INK2, LABEL, MAP_KW, MEMBER_STYLE, MEMBERS,
                                OUT_DIR, RCP_LABELS, RCPS, RVA_DIR, TOPO_CSV, XLIM, YLIM,
                                get_gb_boundary, save)

ROOT = os.path.join(HERE, "..")
BENCH = os.path.join(ROOT, "IHA_Results", "benchmark")
PMEANS = os.path.join(ROOT, "IHA_Results", "period_means", "period_means.csv.gz")
TABLE_DIR = os.path.join(OUT_DIR, "tables")
LOW = ["min_1day", "min_3day", "min_7day", "min_30day", "min_90day", "base_flow_index",
       "month_10_Jul", "month_11_Aug", "month_12_Sep"]
HIGH = ["max_1day", "max_3day", "max_7day", "max_30day", "max_90day"]
RCP_SHADE = {26: "#c6dbef", 45: "#6baed6", 60: "#2171b5", 85: "#08306b"}


def period_changes():
    """Per rcp, member, catchment: changes in each variable/season (abs for tas, % otherwise)."""
    d = pd.read_csv(PMEANS, dtype={"member": str})
    d["member"] = d.member.str.zfill(2)
    w = d.pivot_table(index=["rcp", "member", "catchment_id", "variable", "season"],
                      columns="period", values="value").reset_index()
    w["change"] = np.where(w.variable == "tas", w.future - w.baseline,
                           100 * (w.future / w.baseline - 1))
    return d, w


def gb_means(d):
    """Catchment-averaged baseline/future values per rcp, member, variable, season."""
    g = d.groupby(["rcp", "member", "variable", "season", "period"]).value.mean().unstack("period")
    g["change"] = np.where(g.index.get_level_values("variable") == "tas", g.future - g.baseline,
                           100 * (g.future / g.baseline - 1))
    return g.change


# ---------------------------------------------------------------------------
def fig_benchmark():
    lf = pd.read_csv(os.path.join(BENCH, "lowflow_metrics.csv"))
    summ = pd.read_csv(os.path.join(BENCH, "benchmark_rva_summary.csv"))
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.6], hspace=0.32, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    for period, ls in [("calibration", "-"), ("validation", "--")]:
        for metric, col, lab in [("kge", "#2a78d6", "KGE"), ("kge_log", "#eb6834", "KGE (log flows)")]:
            v = np.sort(lf.loc[lf.period == period, metric].dropna().values)
            ax.plot(v, np.arange(1, len(v) + 1) / len(v), ls, color=col, lw=1.4,
                    label=f"{lab}, {period}")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Efficiency")
    ax.set_ylabel("Fraction of catchments")
    ax.grid(color=GRID, lw=0.5)
    ax.legend(fontsize=6, loc="upper left")
    ax.set_title("(a)", loc="left", fontsize=9, color=INK)

    ax = fig.add_subplot(gs[0, 1])
    val = lf[lf.period == "validation"]
    for metric, col, lab in [("q95_bias_pct", "#2a78d6", "Q95 bias"), ("flv_pct", "#eb6834", "FLV")]:
        v = np.sort(val[metric].dropna().clip(-100, 100).values)
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, lw=1.4, label=lab)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlim(-100, 100)
    ax.set_xlabel("Bias, validation period (%)")
    ax.grid(color=GRID, lw=0.5)
    ax.legend(fontsize=6.5, loc="upper left")
    ax.set_title("(b)", loc="left", fontsize=9, color=INK)

    ax = fig.add_subplot(gs[1, :])
    hm = summ[summ.comparison == "hydrological_model"].set_index("parameter")
    order = list(hm.sort_values("pct_low").index)
    y = np.arange(len(order))
    styles = [("hydrological_model", "#2a78d6", "o", "Observed vs HBV (observed forcing)"),
              ("climate_forcing", "#eb6834", "s", "HBV (observed forcing) vs HBV (CHESS-SCAPE)"),
              ("full_chain", "#1baf7a", "D", "Observed vs HBV (CHESS-SCAPE)")]
    for k, (comp, col, mk, lab) in enumerate(styles):
        s = summ[summ.comparison == comp].set_index("parameter").loc[order]
        yy = y + (k - 1) * 0.25
        if comp != "hydrological_model":
            ax.hlines(yy, s.pct_low_min, s.pct_low_max, color=col, lw=0.8)
        ax.plot(s.pct_low, yy, mk, ms=3.6, color=col, mec="white", mew=0.3, label=lab)
    nc = hm.loc[order, "pct_low_nochange"]
    ax.plot(nc, y, "|", ms=9, mew=1.6, color="#7a2a8c", label="No-change level")
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[p] for p in order], fontsize=7)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Catchments with low alteration, baseline period (%)\n"
                  "(CHESS-SCAPE comparisons: mean of members, line = range)")
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=2)
    ax.set_title("(c)", loc="left", fontsize=9, color=INK)
    save(fig, "fig_benchmark")


def fig_member_climate(d):
    ch = gb_means(d)
    summ = pd.read_csv(os.path.join(RVA_DIR, "rva_catchment_member.csv.gz"), dtype={"member": str})
    summ["dep"] = summ.category.isin(["mod_departure", "high_departure"])
    dep = summ.groupby(["rcp", "member", "parameter"]).dep.mean() * 100
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for m in MEMBERS:
        _, mk = MEMBER_STYLE[m]
        for rcp in RCPS:
            dt = ch.loc[(rcp, m, "tas", "annual")]
            kw = dict(marker=mk, ms=6, mfc=RCP_SHADE[rcp], mec=INK, mew=0.5, ls="")
            axes[0].plot(dt, ch.loc[(rcp, m, "pr", "JJA")], **kw)
            axes[1].plot(dt, ch.loc[(rcp, m, "pr", "DJF")], **kw)
            low = dep.loc[(rcp, m)].loc[LOW].mean()
            high = dep.loc[(rcp, m)].loc[HIGH].mean()
            axes[2].plot(dt, low, **kw)
            axes[2].plot(dt, high, marker=mk, ms=5, mfc="white", mec=RCP_SHADE[rcp], mew=1, ls="")
    for ax, yl, t in zip(axes, ["Summer (JJA) precipitation change (%)",
                                "Winter (DJF) precipitation change (%)",
                                "Catchments in departure (%)"], ["(a)", "(b)", "(c)"]):
        ax.set_xlabel("Warming (°C)")
        ax.set_ylabel(yl, fontsize=7.5)
        ax.grid(color=GRID, lw=0.5)
        ax.set_title(t, loc="left", fontsize=9, color=INK)
    for ax in axes[:2]:
        ax.axhline(0, color=INK2, lw=0.6)
    axes[2].set_ylim(0, 100)
    axes[2].text(0.03, 0.97, "filled: low-flow parameters\nopen: annual maxima", transform=axes[2].transAxes,
                 va="top", fontsize=6.5, color=INK2)
    handles = [Line2D([], [], marker=MEMBER_STYLE[m][1], ls="", ms=5, mfc="white", mec=INK,
                      label=f"Member {m}") for m in MEMBERS]
    handles += [Line2D([], [], marker="o", ls="", ms=6, mfc=RCP_SHADE[r], mec=INK, mew=0.5,
                       label=lab) for r, lab in zip(RCPS, RCP_LABELS)]
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=6.5, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, "fig_member_climate")


MAP_VARS = [("pr", "JJA", "Summer precipitation (%)", "BrBG", 40),
            ("pet", "annual", "Annual PET (%)", "Oranges", None),
            ("aet", "annual", "Annual AET (%)", "BrBG", 20),
            ("sm", "JJA", "Summer soil moisture (%)", "BrBG", 40)]


def _map(ax, gb, coords, values, cmap, vlim):
    gb.plot(ax=ax, **MAP_KW)
    v = coords.join(values.rename("v"), how="inner")
    if vlim is None:
        sc = ax.scatter(v.x, v.y, c=v.v, cmap=cmap, s=5, vmin=0, vmax=np.nanpercentile(v.v, 98),
                        linewidths=0.1, edgecolors="white")
    else:
        sc = ax.scatter(v.x, v.y, c=v.v, cmap=cmap, s=5, norm=TwoSlopeNorm(0, -vlim, vlim),
                        linewidths=0.1, edgecolors="white")
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    return sc


def fig_change_maps(w, coords, gb):
    from matplotlib.ticker import MaxNLocator
    fig = plt.figure(figsize=(7.2, 3.4))
    width, gap = 0.235, 0.02
    for k, (var, season, title, cmap, vlim) in enumerate(MAP_VARS):
        left = 0.02 + k * (width + gap)
        ax = fig.add_axes([left, 0.17, width, 0.72])
        s = w[(w.rcp == 85) & (w.variable == var) & (w.season == season)]
        med = s.groupby("catchment_id").change.median()
        sc = _map(ax, gb, coords, med, cmap, vlim)
        cax = fig.add_axes([left + 0.025, 0.1, width - 0.05, 0.03])
        cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cb.ax.xaxis.set_major_locator(MaxNLocator(nbins=4, symmetric=vlim is not None))
        cb.ax.tick_params(labelsize=6.5)
        cb.outline.set_visible(False)
        fig.text(left + width / 2, 0.89, title, ha="center", va="bottom", fontsize=8, color=INK)
    save(fig, "fig_change_maps")


def figS_change_maps_all(w, coords, gb):
    rows = [("tas", "annual", "Annual temperature (°C)", "Reds", None)] + MAP_VARS
    fig, axes = plt.subplots(len(rows), 4, figsize=(7.2, 2.3 * len(rows)))
    for i, (var, season, title, cmap, vlim) in enumerate(rows):
        s = w[(w.variable == var) & (w.season == season)]
        lim = None if vlim is None else vlim
        for j, rcp in enumerate(RCPS):
            med = s[s.rcp == rcp].groupby("catchment_id").change.median()
            sc = _map(axes[i, j], gb, coords, med, cmap, lim)
            if i == 0:
                axes[i, j].set_title(RCP_LABELS[j], fontsize=9, color=INK)
        cb = fig.colorbar(sc, ax=axes[i, :], orientation="vertical", fraction=0.02, pad=0.01)
        cb.set_label(title, fontsize=7)
        cb.ax.tick_params(labelsize=6.5)
    save(fig, "figS_change_maps_all")


# ---------------------------------------------------------------------------
def table1(d):
    ch = gb_means(d)
    rows = []
    for m in MEMBERS:
        r = {"member": m}
        for key, var, season, fmt in [("dT", "tas", "annual", "{:.1f}"), ("dP_annual", "pr", "annual", "{:+.0f}"),
                                      ("dP_DJF", "pr", "DJF", "{:+.0f}"), ("dP_JJA", "pr", "JJA", "{:+.0f}")]:
            r[key] = " / ".join(fmt.format(ch.loc[(rcp, m, var, season)]) for rcp in (26, 85))
        rows.append(r)
    os.makedirs(TABLE_DIR, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(TABLE_DIR, "table1_members.csv"), index=False)


def si_tables():
    lf = pd.read_csv(os.path.join(BENCH, "lowflow_metrics.csv"))
    rows = []
    names = [("kge", "KGE"), ("kge_log", "KGE, log-transformed flows"), ("q95_bias_pct", "Q95 bias (%)"),
             ("q70_bias_pct", "Q70 bias (%)"), ("flv_pct", "Low-flow volume bias, FLV (%)"),
             ("recession_ratio_obs", "Recession ratio, observed"),
             ("recession_ratio_sim", "Recession ratio, simulated")]
    for col, lab in names:
        r = [lab]
        for period in ["calibration", "validation", "baseline"]:
            v = lf.loc[lf.period == period, col].dropna()
            r.append(f"{v.median():.2f} ({v.quantile(0.25):.2f} to {v.quantile(0.75):.2f})")
        rows.append(r)
    os.makedirs(TABLE_DIR, exist_ok=True)
    pd.DataFrame(rows, columns=["Metric", "Calibration", "Validation", "Baseline"]).to_csv(
        os.path.join(TABLE_DIR, "benchmark_lowflow.csv"), index=False)
    summ = pd.read_csv(os.path.join(BENCH, "benchmark_rva_summary.csv")).set_index(["comparison", "parameter"])
    rows = []
    for p in [p for p in LABEL]:
        r = [LABEL[p]]
        for comp in ["hydrological_model", "climate_forcing", "full_chain"]:
            x = summ.loc[(comp, p)]
            r.append(f"{x.pct_low:.0f}")
            r.append(f"{x.cliffs_delta:+.2f}")
        r.append(f"{summ.loc[('hydrological_model', p), 'pct_low_nochange']:.0f}")
        rows.append(r)
    header = ["Parameter", "Model: low (%)", "Model: delta", "Forcing: low (%)", "Forcing: delta",
              "Chain: low (%)", "Chain: delta", "No change: low (%)"]
    pd.DataFrame(rows, columns=header).to_csv(os.path.join(TABLE_DIR, "benchmark_rva.csv"), index=False)


def main():
    topo = pd.read_csv(TOPO_CSV).set_index("gauge_id")
    coords = topo[["gauge_easting", "gauge_northing"]].rename(columns={"gauge_easting": "x",
                                                                      "gauge_northing": "y"})
    gb = get_gb_boundary(os.path.join(HERE, "_gb_boundary.gpkg"),
                         os.path.join(HERE, "_ne_10m_countries.geojson"))
    fig_benchmark()
    d, w = period_changes()
    fig_member_climate(d)
    fig_change_maps(w, coords, gb)
    figS_change_maps_all(w, coords, gb)
    table1(d)
    si_tables()


if __name__ == "__main__":
    main()
