__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import os
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import networkx as nx
from collections import deque
from Swing.util.Evaluator import Evaluator
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams, gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from scipy.stats import fisher_exact, linregress, ttest_rel, mannwhitneyu, ttest_ind, pearsonr, zscore, wilcoxon
import palettable

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'

def get_swing_list(pkl):
    x = pd.read_pickle(pkl)
    edges = x.edge_list
    full_edges = list(set([e for sub in edges for e in sub]))
    ranking = np.array([[sub.index(e) for e in full_edges] for sub in edges])
    edge_rankings = pd.DataFrame(ranking.T, index=full_edges)
    return edge_rankings

def get_true_edges(gold_filename):
    evaluator = Evaluator(gold_filename, '\t')
    edges = evaluator.gs_flat.tolist()
    return edges, evaluator

if __name__ == '__main__':
    swing_rankings = get_swing_list('gardner_RF_w7_results.pkl')
    mean_swing = swing_rankings.mean(axis=1).sort_values().to_frame().reset_index()
    # mean_swing = swing_rankings[15].sort_values().to_frame().reset_index()
    mean_swing.columns = ['regulator-target', 'rank_importance']
    mean_swing['rank'] = np.arange(1, len(mean_swing)+1)


    params = {'axes.labelsize': 16, 'font.size': 14, 'legend.fontsize': 14,
              'xtick.labelsize': 14, 'ytick.labelsize': 14, 'text.usetex': False}
    rcParams.update(params)

    colors1 = palettable.colorbrewer.qualitative.Set1_8.mpl_colors
    colors2 = palettable.colorbrewer.qualitative.Set2_8.mpl_colors
    colors3 = palettable.colorbrewer.qualitative.Set3_8.mpl_colors
    split_c2 = [colors2[i] for i in [0, 2]]
    paired = palettable.colorbrewer.qualitative.Paired_10.mpl_colors

    gs = "../../data/invitro/gardner_goldstandard.tsv"
    true_edges, ee = get_true_edges(gs)

    _, _, swing_roc = ee.calc_roc(mean_swing)
    swing_roc = swing_roc.values[-1]

    _, _, swing_pr = ee.calc_pr(mean_swing)
    swing_pr = swing_pr.values[-1]

    data = get_experiment_list(gs.replace('goldstandard', 'timeseries'), 14, 1.0)
    gene_list = data[0].columns.values.tolist()
    xcorr_array = xcorr_experiments(data)
    edge_lags = calc_edge_lag(xcorr_array, gene_list, 1, 0.3, timestep=1)
    edge_lags = edge_lags.set_index('Edge').iloc[:, 2]
    true_lags = edge_lags[edge_lags.index.isin(true_edges)]

    roc_list = pd.read_pickle('gardner_rf_roc.pkl')
    pr_list = pd.read_pickle('gardner_rf_pr.pkl')
    rankings = pd.read_pickle('gardner_rf_rankings.pkl')
    # s_roc_list = pd.read_pickle('gardner_swing_rf_roc.pkl')
    # s_pr_list = pd.read_pickle('gardner_swing_rf_pr.pkl')
    # s_rankings = pd.read_pickle('gardner_swing_rf_rankings.pkl')
    avg_rank = pd.DataFrame()
    # s_avg_rank = pd.DataFrame()
    avg_rank['regulator-target'] = rankings[0].sort_values('regulator-target').values[:, 0]
    # s_avg_rank['regulator-target'] = s_rankings[0].sort_values('regulator-target').values[:, 0]
    rank = np.array([rankings[ii].sort_values('regulator-target')['Rank'].values for ii in range(len(rankings))]).T
    avg_rank['rank'] = np.mean(rank, axis=1)
    # s_rank = np.array([s_rankings[ii].sort_values('regulator-target')['Rank'].values for ii in range(len(s_rankings))]).T
    # s_avg_rank['rank'] = np.mean(s_rank, axis=1)

    s_avg_rank = mean_swing.copy()
    avg_rank.set_index('regulator-target', drop=True, inplace=True)
    s_avg_rank.set_index('regulator-target', drop=True, inplace=True)
    avg_rank['s_rank'] = s_avg_rank['rank']
    avg_rank.sort_values('rank', inplace=True)
    avg_rank['rank'] = range(len(avg_rank))
    avg_rank.sort_values('s_rank', inplace=True)
    avg_rank['s_rank'] = range(len(avg_rank))

    avg_rank = pd.concat([edge_lags, avg_rank], axis=1, join='inner')
    te_ranks = avg_rank[avg_rank.index.isin(true_edges)]
    # avg_rank.reset_index(inplace=True)
    # auroc = ee.calc_roc(avg_rank.sort_values('rank'))[2].values[-1]
    # s_auroc = ee.calc_roc(avg_rank.sort_values('s_rank'))[2].values[-1]
    # avg_rank.set_index('regulator-target', drop=True, inplace=True)

    # Promoted
    te_p = te_ranks[te_ranks['s_rank'] < te_ranks['rank']]

    # Demoted
    te_d = te_ranks[(te_ranks['s_rank'] > te_ranks['rank'])]

    #Same
    te_s = te_ranks[(te_ranks['s_rank'] == te_ranks['rank'])]
    fe_ranks = avg_rank[~avg_rank.index.isin(true_edges)]

    fig = plt.figure(figsize=(10, 9))
    left_int = 2
    right_int = 2
    cols = left_int + right_int
    ax1 = plt.subplot2grid((2, cols), (0, 0), rowspan=2, colspan=left_int)
    ax2 = plt.subplot2grid((2, cols), (0, left_int), colspan=right_int)
    ax3 = plt.subplot2grid((2, cols), (1, left_int), colspan=right_int)
    lw = 3
    shift_down = ['ruvA', 'polB', 'recA']
    shift_up = ['umuDC']
    shift_pad = 0.2

    for row in fe_ranks.iterrows():
        edge = row[0]
        lag = row[1]["Lag"]
        start = row[1].values[1]
        end = row[1].values[2]
        verts = [(0, start),
                 (0.5, start),
                 (0.5, end),
                 (1, end)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', lw=2, ec='0.75', zorder=0)
        ax1.add_patch(patch)

    for row in te_d.iterrows():
        edge = row[0]
        lag = row[1]["Lag"]
        start = row[1].values[1]
        end = row[1].values[2]
        verts = [(0, start),
                 (0.5, start),
                 (0.5, end),
                 (1, end)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

        path = Path(verts, codes)
        ec = paired[5]
        patch = patches.PathPatch(path, facecolor='none', lw=lw, ec=ec)
        ax1.add_patch(patch)
        va = 'center'
        if edge[1] in shift_down:
            va = 'top'
            start -= shift_pad
        elif edge[1] in shift_up:
            va = 'bottom'
            start += shift_pad
        if lag > 0:
            edge_str = edge[0] + u"\u2192" + edge[1] + "\nLag = " + str(int(lag))
        else:
            edge_str = edge[0] + u"\u2192" + edge[1]
        ax1.text(-0.02, start, edge_str, ha='right', va=va)

    for row in te_s.iterrows():
        edge = row[0]
        lag = row[1]["Lag"]
        start = row[1].values[1]
        end = row[1].values[2]
        verts = [(0, start),
                 (0.5, start),
                 (0.5, end),
                 (1, end)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

        path = Path(verts, codes)
        ec = 'k'
        patch = patches.PathPatch(path, facecolor='none', lw=lw, ec=ec, zorder=0)
        ax1.add_patch(patch)
        va = 'center'
        if edge[1] in shift_down:
            va = 'top'
            start -= shift_pad
        elif edge[1] in shift_up:
            va = 'bottom'
            start += shift_pad
        if lag > 0:
            edge_str = edge[0] + u"\u2192" + edge[1] + "\nLag = " + str(int(lag))
        else:
            edge_str = edge[0] + u"\u2192" + edge[1]
        ax1.text(-0.02, start, edge_str, ha='right', va=va)

    for row in te_p.iterrows():
        edge = row[0]
        lag = row[1]["Lag"]
        start = row[1].values[1]
        end = row[1].values[2]
        verts = [(0, start),
                 (0.5, start),
                 (0.5, end),
                 (1, end)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

        path = Path(verts, codes)
        ec = paired[1]
        if edge == ('lexA', 'umuDC'):
            ec = paired[3]
        patch = patches.PathPatch(path, facecolor='none', lw=lw, ec=ec)
        ax1.add_patch(patch)
        if lag > 0:
            edge_str = edge[0] + u"\u2192" + edge[1] + "\nLag = " + str(int(lag))
        else:
            edge_str = edge[0] + u"\u2192" + edge[1]
        va = 'center'
        if edge[1] in shift_down:
            va = 'top'
            start -= shift_pad
        elif edge[1] in shift_up:
            va = 'bottom'
            start += shift_pad
        ax1.text(-0.02, start, edge_str, ha='right', va=va)

    ax1.set_ylim([57, -1])
    ax1.plot([0, 0], [-1, 57], '-', c='k', lw=1)
    # ax1.spines['right'].set_visible = False
    ax1.spines['top'].set_bounds(0, 1)
    ax1.spines['bottom'].set_bounds(0, 1)
    ax1.spines['left'].set_position(('data', 0))
    ax1.set_xlim([-0.6, 1])
    ax1.yaxis.tick_right()
    ax1.set_ylabel('Rank', rotation=0, ha='center')
    ax1.yaxis.set_label_coords(1.005, 1.005)
    ax1.tick_params(axis='x', which='both', bottom='off', top='off')
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['RF', 'SWING-RF'])

    # In experiment 4, G2 and nothing else upstream is perturbed, so the relation between G2, G1 is clearer
    # The apparent lag is 2, so shift it that much
    # Plot unlagged time series and correlation
    time = data[0].index.values
    lexa = zscore(data[0]['lexA'].values, ddof=1)
    umudc = zscore(data[0]['umuDC'].values, ddof=1)

    ax2.plot(time, lexa, '.-', label='lexA', c=paired[3], lw=3, ms=15)
    ax2.plot(time, umudc, '.--', dashes=(3, 2), label='umuDC', c='k', lw=3, ms=15, mfc='w', mew=2)
    ax2.plot(time[:-1], umudc[1:], '.-', label='umuDC-shifted', c='k', lw=3, ms=15)
    arrow_nums = [1, 2, 3, 7, 10, 12, 13]
    pad = 0.2
    for idx, point in enumerate(time):
        if idx > 0 and idx in arrow_nums:
            ax2.annotate("", xy=(time[idx-1]+pad, umudc[idx]+.001),
                         xytext=(point-pad, umudc[idx]+.001),
                         arrowprops=dict(facecolor='k', headlength=8, headwidth=7, width=1.45))
    # ax2.set_yticks(np.arange(0.0, 0.6, 0.1))
    # ax2.set_xticks(np.arange(0.0, 1200, 200))
    ax2.set_ylabel('Normalized expression')
    ax2.set_xlabel('Time')
    ax2.legend(loc='best', numpoints=2)

    reg_fit = linregress(lexa, umudc)
    shift_fit = linregress(lexa[:-1], umudc[1:])

    ax3.plot(lexa, umudc, '.', c='k', mfc='w', mew=2, ms=15, label='umuDC')
    ax3.plot(lexa[:-1], umudc[1:], '.', c='k', ms=15, label='umuDC-shifted')
    xvals = np.array([np.min(lexa), np.max(lexa)])
    ax3.plot(xvals, xvals * reg_fit.slope + reg_fit.intercept,
             c='k', ls='--', dashes=(5, 5), lw=1, zorder=0, label=('$\mathregular{R^2 = %0.3f}$' % reg_fit.rvalue))
    xvals = np.array([np.min(lexa[:-1]), np.max(lexa[:-1])])
    ax3.plot(xvals, xvals * shift_fit.slope + shift_fit.intercept,
             c='k', zorder=1, label=('$\mathregular{R^2 = %0.3f}$' % shift_fit.rvalue))
    # ax3.set_yticks(np.arange(0.0, 0.4, 0.1))
    # ax3.set_xticks(np.arange(0.05, 0.5, 0.1))
    ax3.set_xlabel('lexA normalized expression')
    ax3.set_ylabel('Normalized expression')
    ax3.legend(loc='best', numpoints=1, ncol=2, handletextpad=0.5, handlelength=1, columnspacing=0.5)

    plt.tight_layout()
    plt.savefig('gardner_sos_edge_promotion_w7.pdf', fmt='pdf')
