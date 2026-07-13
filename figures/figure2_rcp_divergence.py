"""
figure2_rcp_divergence.py
=========================
Figure 2: Annual median absolute percentage difference in ensemble mean
HBV-simulated discharge between RCP2.6 and RCP8.5, on a water-year basis,
used to verify the scenario-neutral baseline period (WY1982-WY2010).

Inputs (edit paths as needed):
    rcp26_hbv_ensemble_mean_discharge.csv
    rcp85_hbv_ensemble_mean_discharge.csv

Output:
    figure2_rcp_divergence.png  (300 dpi)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

RCP26_PATH = '/mnt/user-data/uploads/rcp26_hbv_ensemble_mean_discharge.csv'
RCP85_PATH = '/mnt/user-data/uploads/rcp85_hbv_ensemble_mean_discharge.csv'
OUT_PATH   = '/mnt/user-data/outputs/figure2_rcp_divergence.png'

BASELINE_END_WY = 2010   # end of scenario-neutral baseline (verified zero-divergence)
FUTURE_START_WY = 2051
FUTURE_END_WY   = 2080


def wy_stats(df26, df85, wy):
    """Median/IQR of |% difference| between RCP2.6 and RCP8.5 for one water year."""
    dates = df26.index
    months_this_wy = [f'{wy}-0{m}' if m < 10 else f'{wy}-{m}' for m in range(1, 10)]
    months_prev_wy = [f'{wy-1}-{m}' for m in (10, 11, 12)]
    mask = pd.Series(False, index=dates)
    for m in months_this_wy + months_prev_wy:
        mask |= dates.str.startswith(m)
    m26 = df26[mask].mean(axis=0)
    m85 = df85[mask].mean(axis=0)
    pct = 100 * (m85 - m26).abs() / m26.replace(0, np.nan)
    return pct.median(), pct.quantile(0.25), pct.quantile(0.75)


def main():
    rcp26 = pd.read_csv(RCP26_PATH, index_col=0)
    rcp85 = pd.read_csv(RCP85_PATH, index_col=0)

    wys, meds, q25s, q75s = [], [], [], []
    for wy in range(1982, 2081):
        med, q25, q75 = wy_stats(rcp26, rcp85, wy)
        wys.append(wy); meds.append(med); q25s.append(q25); q75s.append(q75)
    wys = np.array(wys); meds = np.array(meds); q25s = np.array(q25s); q75s = np.array(q75s)

    COL_LINE, COL_BAND = '#1a6faf', '#a8c8e8'
    fig, ax = plt.subplots(figsize=(6.85, 3.6), dpi=300)

    ax.axvspan(1982, BASELINE_END_WY + 1, color='#2c7bb6', alpha=0.07, zorder=0, lw=0)
    ax.axvspan(FUTURE_START_WY, FUTURE_END_WY + 1, color='#d7191c', alpha=0.07, zorder=0, lw=0)

    # Step plot: each water year is a discrete annual value, not a continuous
    # quantity, so a step function (not linear interpolation) is the honest
    # representation, especially right at the baseline/divergence boundary.
    ax.fill_between(wys, q25s, q75s, step='post', color=COL_BAND, alpha=0.6, zorder=1, lw=0)
    ax.plot(wys, meds, drawstyle='steps-post', color=COL_LINE, lw=1.4, zorder=2)

    ax.axvline(BASELINE_END_WY + 1, color='#333333', lw=0.9, ls='--', zorder=4)

    ax.set_xlim(1982, 2080)
    ax.set_ylim(-1, 52)

    ax.text(1996.5, 49, f'Baseline period\n(WY1982\u2013WY{BASELINE_END_WY})',
            ha='center', va='top', fontsize=6.5, color='#1a5276', fontstyle='italic')
    ax.text(2066, 49, f'Future period\n(WY{FUTURE_START_WY}\u2013WY{FUTURE_END_WY})',
            ha='center', va='top', fontsize=6.5, color='#922b21', fontstyle='italic')

    ax.set_xlabel('Water year', fontsize=8)
    ax.set_ylabel('Median absolute discharge\ndifference (%)', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', lw=0.4, alpha=0.4, color='#cccccc')

    band_p = mpatches.Patch(color=COL_BAND, alpha=0.6, label='Interquartile range across catchments')
    line_p = plt.Line2D([0], [0], color=COL_LINE, lw=1.4, label='Annual median')
    ax.legend(handles=[line_p, band_p], fontsize=6.5, framealpha=0.92, edgecolor='#cccccc',
              loc='upper right', bbox_to_anchor=(1.0, -0.18), ncol=2, borderaxespad=0)

    plt.subplots_adjust(left=0.11, right=0.97, top=0.97, bottom=0.22)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved {OUT_PATH}')


if __name__ == '__main__':
    main()
