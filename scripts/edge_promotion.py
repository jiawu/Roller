__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import os
import pandas as pd
import numpy as np
import networkx as nx
from collections import deque
from Swing.util.Evaluator import Evaluator
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Polygon
import brewer2mpl
from scipy.stats import fisher_exact, linregress, ttest_rel, mannwhitneyu, ttest_ind, pearsonr


def is_square(n):
    """
    Determine if a number is a perfect square
    :param n: int or float
        The number to check
    :return: Boolean
        Return True if the number is a perfect square
    """
    return np.sqrt(n).is_integer()


def get_factors(n):
    """
    Calculate the factors of a number
    :param n: int
        The number to be factored
    :return: list
        A sorted list of the unique factors from smallest to largest
    """
    factor_list = np.array([[i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0]).flatten().astype(int)
    return sorted(factor_list.tolist())


def calc_subplot_dimensions(x):
    """
    Calculate the dimensions for a matplotlib subplot object.
    :param x: int
        Number of plots that need to be made
    :return: rows, columns
        The number of rows and columns that should be in the subplot
    """
    if x <= 3:
        rows = x
        columns = 1
    else:
        factor_list = get_factors(x)
        while len(factor_list) <= 2 and not is_square(x):
            x += 1
            factor_list = get_factors(x)
        if is_square(x):
            rows = int(np.sqrt(x))
            columns = int(np.sqrt(x))

        else:
            rows = factor_list[int(len(factor_list)/2-1)]
            columns = factor_list[int(len(factor_list)/2)]

    return rows, columns


def get_true_edges(gold_filename):
    evaluator = Evaluator(gold_filename, '\t')
    edges = evaluator.gs_flat.tolist()
    return edges, evaluator


def get_edge_lags(data_filename):
    df = pd.read_csv(data_filename, sep="\t")
    gene_list = df.columns.values[1:].tolist()
    experiment_list = get_experiment_list(data_filename, 21, 10)
    xcorr_array = xcorr_experiments(experiment_list)
    lags = calc_edge_lag(xcorr_array, gene_list, 0.1, 0.5, timestep=1)
    return lags, df


def get_network_changes(pickle_filename, edge_str='regulator-target', base_str='rank_importance_RF-td_21',
                        shortener_str='rank_importance_', replace=''):
    results_df = pd.read_pickle(pickle_filename)
    edges = results_df[edge_str].values
    baseline = results_df[base_str].values

    rank_df = pd.DataFrame()
    rank_df[edge_str] = edges
    rank_df[('Base_%s' %replace)] = baseline
    for column in results_df.columns:
        if column != edge_str and column!=base_str:
            short_name = column.replace(shortener_str, replace)
            rank_df[short_name] = results_df[column].values
    rank_df.set_index(['regulator-target'], inplace=True)
    diff_df = (rank_df.T.iloc[0]-rank_df.T).T
    parameters = set(rank_df.columns[1:].values)
    return diff_df, rank_df, parameters


def get_network_data(goldstandard, timeseries, ignore_self=True):
    # Get true network
    true_edges, evaluator = get_true_edges(goldstandard)
    dg = nx.DiGraph()
    dg.add_edges_from(true_edges)

    #Network statistics - deprecated
    #degree = nx.degree_centrality(dg)
    #b_cent = pd.DataFrame.from_dict({k: [v] for k, v in nx.edge_betweenness_centrality(dg).items()}, 'index')
    #b_cent.columns = ['Bcent']

    #Calculate edge lags
    edge_lags, data = get_edge_lags(timeseries)
    if ignore_self:
        edge_lags = edge_lags[edge_lags['Parent'] != edge_lags['Child']]
    edge_df = pd.DataFrame(edge_lags['Lag'].values, index=edge_lags['Edge'].values, columns=['Lag'])

    return true_edges, edge_df, data, dg, evaluator


def get_signed_edges(signed):
    df = pd.read_csv(signed, sep='\t', header=None)
    df['regulator-target'] = list(zip(df[0], df[1]))
    df.set_index(['regulator-target'], inplace=True)
    df.drop([0, 1], axis=1, inplace=True)
    df.columns=['sign']
    return df


def calc_scores(ranking_df, evaluator):
    filtered_ranks = ranking_df.copy()
    filtered_ranks.reset_index(level=0, inplace=True)
    roc = []
    aupr = []
    for c in filtered_ranks.columns[1:]:
        filtered_ranks.sort_values(filtered_ranks.columns[1], inplace=True)
        roc.append(evaluator.calc_roc(filtered_ranks.iloc[:, :2])[2].values[-1])
        aupr.append(evaluator.calc_pr(filtered_ranks.iloc[:, :2])[2].values[-1])
        filtered_ranks.drop(c, axis=1, inplace=True)
    return roc, aupr


def calc_promotion(change_df, columns):
    t_promoted = np.sum(change_df.loc[:, columns].values > 0, axis=0)
    t_demoted = np.sum(change_df.loc[:, columns].values < 0, axis=0)
    t_same = np.sum(change_df.loc[:, columns].values == 0, axis=0)

    t_lagged = change_df[change_df['Lag'] > 0]
    l_promoted = np.sum(t_lagged.loc[:, columns].values > 0, axis=0)
    l_demoted = np.sum(t_lagged.loc[:, columns].values < 0, axis=0)
    l_same = np.sum(t_lagged.loc[:, columns].values == 0, axis=0)
    rows = ['true+', 'true-', 'true=', 'lag+', 'lag-', 'lag=']
    return pd.DataFrame([t_promoted, t_demoted, t_same, l_promoted, l_demoted, l_same], index=rows, columns=columns).T


def get_net_stats(dg):
    g = dg.to_undirected()
    assort = nx.degree_pearson_correlation_coefficient(dg)
    if np.isnan(assort):
        assort = 0
    clust = nx.average_clustering(g)
    trans= nx.transitivity(g)
    try:
        rad = nx.radius(g)
    except nx.NetworkXError:
        rad = 0
    try:
        diam = nx.diameter(g)
    except nx.NetworkXError:
        diam = 0
    return [assort, clust, trans, rad, diam]


def get_c_table(summary_df):
    """
    C table format: [[lagged_promoted, lagged_not_promoted],
                 [not_lagged_promoted, not_lagged_not_promoted]]
    """

    c_table = np.array([[summary_df['lag+'], summary_df['lag-']+summary_df['lag=']],
                        [summary_df['true+']-summary_df['lag+'],
                         summary_df['true-']+summary_df['true=']-summary_df['lag-']-summary_df['lag=']]])
    c_table = np.array([c_table[:, :, ii] for ii in range(c_table.shape[2])])
    return c_table


def get_enrichment(c_array, conditions):
    pvals = pd.DataFrame([fisher_exact(c_tab)[1] for c_tab in c_array], index=conditions, columns=['pval'])
    return pvals


def make_dictionary(methods, replace_dict, models, num_nets=20, directory_path="../data/gnw_insilico/network_data/"):
    lag_range = {'ml_0': [0, 1], 'ml_1': [0, 2], 'ml_2': [0, 3], 'ml_3': [1, 2], 'ml_4': [1, 4], 'ml_5': [2, 3]}

    s_dict = {}
    for ii, method in enumerate(methods):
        s_dict[method] = {"aupr": pd.DataFrame(), "auroc": pd.DataFrame(),
                          "te_change": pd.DataFrame(), "te_rank": pd.DataFrame()}
        for model in models:
            s_dict[method][model] = {"aupr": pd.DataFrame(), "auroc": pd.DataFrame(),
                                     "te_change": pd.DataFrame(), "te_rank": pd.DataFrame()}
            for net in range(1, num_nets + 1):
                s_dict[method][model][net] = {"aupr": pd.DataFrame(), "auroc": pd.DataFrame(),
                                              "te_change": pd.DataFrame(), "te_rank": pd.DataFrame()}
                short = 'rank_importance_%s' % method
                pickle_file = "%s_net%i_%s_promotion.pkl" % (model, net, method.lower())
                base_str = ('rank_importance_%s-td_21' % method)

                roc_df = pd.DataFrame()
                pr_df = pd.DataFrame()
                # Get the network information
                gold_file = directory_path+"%s/%s-%i_goldstandard.tsv" % (model, model, net)
                signed_file = gold_file.replace('.tsv', '_signed.tsv')
                data_file = directory_path+"%s/%s-%i_timeseries.tsv" % (model, model, net)
                true_edges, edge_df, data, dg, evaluator = get_network_data(gold_file, data_file)
                signed_edges = get_signed_edges(signed_file)

                change, ranks, params = get_network_changes(pickle_file, base_str=base_str,
                                                            shortener_str=short, replace=replace_dict[method])

                change_df = change.reindex_axis(sorted(change.columns), axis=1)
                ranks_df = ranks.reindex_axis(sorted(ranks.columns), axis=1)
                conditions = change_df.columns.values
                s_dict[method][model][net]['rank'] = ranks_df
                s_dict[method][model][net]['rank_change'] = change_df

                # Calculate the auroc and aupr for each parameter set of the network
                roc_df[model + str(net)], pr_df[model + str(net)] = calc_scores(ranks_df, evaluator)
                roc_df.index = conditions
                pr_df.index = conditions
                s_dict[method][model][net]['auroc'] = roc_df.T
                s_dict[method][model][net]['aupr'] = pr_df.T

                # Compile results
                full_change = pd.concat([edge_df, change_df], axis=1, join='inner')
                full_rank = pd.concat([edge_df, ranks_df], axis=1, join='inner')
                te_rank = pd.concat([signed_edges, full_rank[full_rank.index.isin(true_edges)]],
                                    axis=1, join='inner')
                te_change = pd.concat([signed_edges, full_change[full_change.index.isin(true_edges)]],
                                      axis=1, join='inner')
                promotions = calc_promotion(te_change, conditions)
                contingency = get_c_table(promotions)
                s_dict[method][model][net]['conditions'] = conditions
                s_dict[method][model][net]['change'] = full_change
                s_dict[method][model][net]['rank'] = full_rank
                s_dict[method][model][net]['te_change'] = te_change
                s_dict[method][model][net]['te_rank'] = te_rank
                s_dict[method][model][net]['promotion'] = promotions
                s_dict[method][model][net]['contingency'] = contingency
                s_dict[method][model][net]['enrich_pvals'] = get_enrichment(contingency, conditions)
                s_dict[method][model][net]['stats'] = pd.DataFrame(get_net_stats(dg), index=['assort', 'clust',
                                                                                             'trans', 'rad', 'diam'])
                # Summarize it for each model organism
                s_dict[method][model]['aupr'] = pd.concat([s_dict[method][model]['aupr'], pr_df.T], join='inner')
                s_dict[method][model]['auroc'] = pd.concat([s_dict[method][model]['auroc'], roc_df.T], join='inner')
                s_dict[method][model]['te_change'] = pd.concat([s_dict[method][model]['te_change'], te_change],
                                                               join='outer')
                s_dict[method][model]['te_rank'] = pd.concat([s_dict[method][model]['te_rank'], te_rank], join='outer')
            auroc = s_dict[method][model]['auroc']
            aupr = s_dict[method][model]['aupr']
            s_dict[method][model]['auroc_diff'] = pd.DataFrame((auroc.values.T - auroc.values[:, 0]).T,
                                                               index=auroc.index, columns=auroc.columns)
            s_dict[method][model]['aupr_diff'] = pd.DataFrame((aupr.values.T - aupr.values[:, 0]).T,
                                                              index=aupr.index, columns=aupr.columns)
            s_dict[method][model]['promotion'] = calc_promotion(s_dict[method][model]['te_change'], conditions)
            s_dict[method][model]['contingency'] = get_c_table(s_dict[method][model]['promotion'])
            s_dict[method][model]['enrich_pvals'] = get_enrichment(s_dict[method][model]['contingency'], conditions)

            # Summarize it for each method
            s_dict[method]['aupr'] = pd.concat([s_dict[method]['aupr'], s_dict[method][model]['aupr']], join='inner')
            s_dict[method]['auroc'] = pd.concat([s_dict[method]['auroc'], s_dict[method][model]['auroc']], join='inner')
            s_dict[method]['te_change'] = pd.concat([s_dict[method]['te_change'],
                                                     s_dict[method][model]['te_change']], join='outer')
            s_dict[method]['te_rank'] = pd.concat([s_dict[method]['te_rank'], s_dict[method][model]['te_rank']],
                                                  join='outer')

        auroc = s_dict[method]['auroc']
        aupr = s_dict[method]['aupr']
        s_dict[method]['auroc_diff'] = pd.DataFrame((auroc.values.T - auroc.values[:, 0]).T,
                                                    index=auroc.index, columns=auroc.columns)
        s_dict[method]['aupr_diff'] = pd.DataFrame((aupr.values.T - aupr.values[:, 0]).T,
                                                   index=aupr.index, columns=aupr.columns)
        s_dict[method]['promotion'] = calc_promotion(s_dict[method]['te_change'], s_dict[method]['aupr'].columns)
        s_dict[method]['contingency'] = get_c_table(s_dict[method]['promotion'])
        s_dict[method]['enrich_pvals'] = get_enrichment(s_dict[method]['contingency'], conditions)
    return s_dict


def stars(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "-"


def draw_boxes(axes, box_obj, color_list, ww=1.5):
    fliers_exist = True if len(box_obj['fliers']) > 0 else False

    # Replace the boxes on the boxplots with nice patches
    for ii, box in enumerate(box_obj['boxes']):
        # Clear the existing box
        box.set_linewidth(0)

        # Add the new box at the coordinates for the existing box
        box_poly = Polygon(box.get_xydata(), facecolor=color_list[ii], linewidth=0)
        axes.add_patch(box_poly)

        # Remove the whiskers, there are 2 of them
        for jj in range(2):
            box_obj['whiskers'][ii * 2 + jj].set_color(color_list[ii])
            box_obj['whiskers'][ii * 2 + jj].set_linewidth(ww)

        # top and bottom fliers
        # (set allows us to set many parameters at once)
        if fliers_exist:
            box_obj['fliers'][ii].set(marker='x', ms=5, mew=1, mec=color_list[ii])

        box_obj['medians'][ii].set_color('black')
        box_obj['medians'][ii].set_linewidth(2)
        box_obj['medians'][ii].set_solid_capstyle('butt')

    # and 4 caps to remove
    for c in box_obj['caps']:
        c.set_linewidth(0)


def plot_scores(axes, ctrl, swing, pval, pos, palette):
    x_array = np.tile(pos, (len(ctrl), 1))
    y_array = np.array([ctrl, swing])
    y_max = np.max(y_array, axis=None)
    y_min = np.min(y_array, axis=None)

    # Plot the paired lines
    axes.plot(x_array.T, y_array, '.-', c='k', alpha=0.4, zorder=10)

    # Add the boxplots
    bp = axes.boxplot([ctrl, swing], positions=pos, widths=0.25)
    s = stars(pval)
    if p_value < 0.05:
        axes.annotate("", xy=(pos[0], y_max + .01), xycoords='data', xytext=(pos[1], y_max + .01), textcoords='data',
                      arrowprops=dict(arrowstyle="-", ec='k', connectionstyle="bar,fraction=0.1"))
        axes.text((pos[0]+pos[1])/2, y_max+.03, s, horizontalalignment='center', verticalalignment='center')

    draw_boxes(axes, bp, palette)
    axes.get_xaxis().tick_bottom()
    axes.get_yaxis().tick_left()
    axes.tick_params(axis='x', direction='out')
    axes.tick_params(axis='y', length=0)
    # axes.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axes.set_axisbelow(True)
    return y_min, y_max


def plot_rank_change(axes, ranking, ctrl_str, swing_str, color, cutoff=90):
    nl = ranking[ranking["Lag"] == 0]
    l = ranking[ranking["Lag"] != 0]
    axes.plot(nl[ctrl_str][nl[ctrl_str] < cutoff], nl[swing_str][nl[ctrl_str] < cutoff],
              'x', c='k', alpha=0.7, label='Not Lagged', ms=5, mew=2, zorder=2)
    axes.plot(l[ctrl_str][l[ctrl_str] < cutoff], l[swing_str][l[ctrl_str] < cutoff],
              '.', c=color, label='Lagged', ms=10, zorder=1)
    axes.plot([0, cutoff], [0, 90], color='k', lw=1, ls='-', zorder=0)


def plot_diff_distribution(axes, promotion, swing_str, pos, palette, width=0.25, ww=3):
    nl = promotion[swing_str][promotion["Lag"] == 0]
    l = promotion[swing_str][promotion["Lag"] != 0]
    bp = axes.boxplot([nl, l], positions=pos, showfliers=False, widths=width)

    y_max = 0
    y_min = 0
    for whisker in bp['whiskers']:
        coords = whisker._xy[:, 1]
        y_max = max(np.max(coords), y_max)
        y_min = min(np.min(coords), y_min)

    draw_boxes(axes, bp, palette, ww=ww)
    pval = mannwhitneyu(nl, l).pvalue
    print(np.mean(l) - np.mean(nl), pval)
    s = stars(pval)
    if pval < 0.05:
        axes.annotate("", xy=(pos[0], y_max + 1), xycoords='data', xytext=(pos[1], y_max + 1), textcoords='data',
                      arrowprops=dict(arrowstyle="-", ec='k', connectionstyle="bar,fraction=0.05"))
        axes.text((pos[0]+pos[1])/2, y_max+1.1, s, horizontalalignment='center')

    # axes.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axes.set_axisbelow(True)
    return y_min, y_max


# def draw_bezier(axes, )

#NOTES FROM JIA
"""
Some columns are labeled as RF-td_2, RF-td10, RF-ml_1, etc.

td_X refers to a window size of X. For all these networks, the rest of the parameters are the same:
min_lag is 1, max lag is 3, except RF-td 21, which is just regular random forest.

ml_X refers to a set of min_lags/max_lags for each X.
if X = 0; min_lag = 0, max_lag = 1
elif X = 1; min_lag = 0, max_lag = 2
elif X = 2; min_lag = 0, max_lag = 3
elif X = 3; min_lag = 1, max_lag = 2
elif X = 4; min_lag = 1, max_lag = 4
elif X = 5; min_lag = 2, max_lag = 3
"""


if __name__ == "__main__":
    # Values for displaying or saving figures are True, 'show' or False
    savefig1 = False
    savefig2 = 'show'
    savefigs_group2 = False
    params = {'axes.labelsize': 16, 'axes.labelweight': 'bold', 'font.size': 14, 'legend.fontsize': 14,
              'xtick.labelsize': 14, 'ytick.labelsize': 14,'text.usetex': False}

    colors1 = list(brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors)
    colors2 = list(brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors)
    colors3 = list(brewer2mpl.get_map('Set3', 'qualitative', 8).mpl_colors)
    split_c2 = [colors2[i] for i in [0, 2]]
    paired = list(brewer2mpl.get_map('Paired', 'qualitative', 10).mpl_colors)
    rcParams.update(params)
    m = ['RF', 'Dionesus']
    rd = {'Dionesus': 'D', 'RF': 'RF'}
    mod = ['Yeast', 'Ecoli']
    try:
        print("Loading_pickle")
        summary = pd.read_pickle('./param_sweep_summary.pickle')
    except FileNotFoundError:
        print("Pickle doesn't exist. Summarizing data now")
        summary = make_dictionary(m, rd, mod)
        pd.to_pickle(summary, './param_sweep_summary.pickle')

    # Overall difference
    swing_rf, base_rf = summary['RF']['auroc']['RF-ml_4'], summary['RF']['auroc']['Base_RF']
    swing_d, base_d = summary['Dionesus']['auroc']['D-ml_4'], summary['Dionesus']['auroc']['Base_D']
    print('AUROC increase, and pvalue')
    print('RF', np.mean(swing_rf)-np.mean(base_rf), (np.mean(swing_rf)-np.mean(base_rf))/np.mean(base_rf),
          ttest_rel(swing_rf, base_rf).pvalue)
    print('Dionesus', np.mean(swing_d) - np.mean(base_d), (np.mean(swing_d) - np.mean(base_d))/np.mean(base_d),
          ttest_rel(swing_d, base_d).pvalue)

    print('\nYeast, AUROC increase, and pvalue')
    swing_d, base_d = summary['Dionesus']['Yeast']['auroc']['D-ml_4'], summary['Dionesus']['Yeast']['auroc']['Base_D']
    print('Dionesus', np.mean(swing_d) - np.mean(base_d), (np.mean(swing_d) - np.mean(base_d)) / np.mean(base_d),
          ttest_rel(swing_d, base_d).pvalue)


    swing_rf, base_rf = summary['RF']['aupr']['RF-ml_4'], summary['RF']['aupr']['Base_RF']
    swing_d, base_d = summary['Dionesus']['aupr']['D-ml_4'], summary['Dionesus']['aupr']['Base_D']
    print('\nAUPR increase and pvalue')
    print('RF', np.mean(swing_rf) - np.mean(base_rf), (np.mean(swing_rf) - np.mean(base_rf)) / np.mean(base_rf),
          ttest_rel(swing_rf, base_rf).pvalue)
    print('Dionesus', np.mean(swing_d) - np.mean(base_d), (np.mean(swing_d) - np.mean(base_d)) / np.mean(base_d),
      ttest_rel(swing_d, base_d).pvalue)

    print('Promotion')
    c_index = 5
    rf_c = summary['RF']['contingency'][c_index]
    d_c = summary['Dionesus']['contingency'][c_index]
    print(rf_c[1, 0]/(rf_c[1,0]+rf_c[1, 1]), rf_c[0, 0]/(rf_c[0,0]+rf_c[0, 1]), summary['RF']['enrich_pvals'].iloc[c_index])
    print(d_c[1, 0] / (d_c[1, 0] + d_c[1, 1]), d_c[0, 0] / (d_c[0, 0] + d_c[0, 1]),
          summary['Dionesus']['enrich_pvals'].iloc[c_index])

    scores = ['auroc', 'aupr']
    bprops = dict(linewidth=3)
    mprops = dict(linewidth=3, color='r')
    wprops = dict(linewidth=3, linestyle='--')

    """
    ####################################################################################################################
    ####################################################################################################################
    Score changes figure
    ####################################################################################################################
    ####################################################################################################################
    """
    # Parameter ml_4 includes the expected lag and has high enrichment. To keep things simpler I will use this
    f = plt.figure(figsize=(10, 9))
    positions = np.reshape(np.arange(2*(len(m)+len(mod))), (len(m)+len(mod), 2))
    network_size = 10
    box_pairs = [[paired[ii * 2 + 6], paired[ii * 2 + 1 + 6]] for ii in range(len(mod))]

    for kk, score in enumerate(scores):
        if not savefig1:
            plt.close(f)
            break
        ax = f.add_subplot(2, 1, kk+1)
        # Initialize parameters
        ax_label, ax_patch, num_te, ax_min, ax_max = [], [], 0, 1, 0
        for ii, model in enumerate(mod):
            for jj, method in enumerate(m):
                # Retrieve/Calculate plotting data
                control = summary[method][model][score].iloc[:, 0]
                test = summary[method][model][score][(rd[method]+'-ml_4')]
                p_value = ttest_rel(control, test).pvalue
                rank = summary[method][model]['te_rank']
                change = summary[method][model]['te_change']

                # Plot results for network scores
                f_min, f_max = plot_scores(ax, control, test, p_value, positions[ii*2+jj], box_pairs[ii])
                ax_min, ax_max = min(ax_min, f_min), max(ax_max, f_max)

                # Legend patches
                if jj == 0:
                    for cc in box_pairs[ii]:
                        c_patch = patches.Patch(color=cc, label=model)
                        ax_patch.append(c_patch)
                ax_label.append(method)
                ax_label.append('SWING\n' + method)

            num_te += len(change)

        # Plot expected null models
        if score == 'auroc':
            ax.plot([np.min(positions) - 1, np.max(positions) + 1], [0.5, 0.5], '--', c='k', zorder=0)
        elif score == 'aupr':
            null_aupr = num_te/len(control)/(network_size**2-network_size)/len(m)
            ax.plot([np.min(positions) - 1, np.max(positions) + 1], [null_aupr, null_aupr], '--', c='k', zorder=0)

        #Format axes
        ax.set_xlim([np.min(positions)-0.5, np.max(positions)+0.5])
        ax.set_ylim([ax_min-0.05, ax_max+0.07])
        if kk == 1:
            ax.set_xticks(positions.flatten())
            ax.set_xticklabels(ax_label, rotation=45, ha='center')
        else:
            ax.set_xticks([])
        ax.set_ylabel(score.upper())
        f.tight_layout()
        f_filename = '../manuscript/Figures/gnw_score_improvement.pdf'

        if savefig1 is True and kk > 0:
            f.savefig(f_filename, fmt='pdf')
        elif savefig1 == 'show' and kk > 0:
            plt.show()
    """
    ####################################################################################################################
    ####################################################################################################################
    Edge promotion figures
    ####################################################################################################################
    ####################################################################################################################
    """
    top = 0.98
    bottom = 0.07
    hspace = 0.1
    g = plt.figure(figsize=(11, 7))
    gs1 = gridspec.GridSpec(len(mod), len(m))
    gs1.update(left=0.06, right=0.6, wspace=0.1, hspace=hspace, top=top, bottom=bottom)
    gs2 = gridspec.GridSpec(len(mod), 1)
    gs2.update(left=gs1.right+0.09, right=0.99, hspace=hspace, top=top, bottom=bottom)

    rank_change_lim = 90
    for ii, model in enumerate(mod):
        if not savefig2:
            plt.close(g)
            break
        hax_label = []
        hax_min, hax_max = 1, 0
        for jj, method in enumerate(m):
            # Retrieve/Calculate plotting data
            rank = summary[method][model]['te_rank']
            change = summary[method][model]['te_change']

            axnum = ii * 2 + jj
            if ii == 0:
                axnum2 = 3
            else:
                axnum2 = 6

            # Add subplot
            gax = plt.subplot(gs1[axnum])
            hax_label.append(method)
            hax_label.append('SWING\n' + method)
            if axnum2:
                hax = plt.subplot(gs2[ii])
                hax.set_ylabel('True Edge Rank Change')
                # Plot results for rank change distributions
                h_pos = positions[:2].flatten()
                h_min, h_max = plot_diff_distribution(hax, change, rd[method] + '-ml_4', positions[jj], split_c2)
                # Format rank difference distribution
                hax.set_xlim([np.min(h_pos) - 0.5, np.max(h_pos)+0.5])
                hax_min, hax_max = min(hax_min, h_min), max(hax_max, h_max)
                hax.set_ylim([hax_min - 0.3 * np.abs(hax_min), hax_max + 0.3 * np.abs(hax_max)])
                hax.set_xticks(h_pos)
                hax.set_xticklabels(hax_label)
                hax.set_ylabel('True Edge Rank Change')

                if ii == 0:
                    nl_patch = patches.Patch(color=split_c2[0], label='Not lagged')
                    l_patch = patches.Patch(color=split_c2[1], label='Lagged')
                    hax.legend(handles=[nl_patch, l_patch], loc='best')

            if axnum == 0:
                gax.set_ylabel('SWING Rank')
            elif axnum == 2:
                gax.set_ylabel('SWING Rank')
                gax.set_xlabel(method + ' Rank')
            elif axnum == 3:
                gax.set_xlabel(method + ' Rank')

            # Plot results for rank changes
            plot_rank_change(gax, rank, 'Base_'+rd[method], rd[method]+'-ml_4', box_pairs[ii][0], rank_change_lim)
            gax.set_yticks(np.arange(0, 100, 30))
            gax.set_xticks(np.arange(0, rank_change_lim+10, 30))
            if axnum == 0 or axnum == 2:
                gax.legend(loc='best', numpoints=1, handletextpad=-0.25)
            if ii == 0:
                gax.set_xticklabels([])
                hax.set_xticklabels([])
            if jj == 1:
                gax.set_yticklabels([])

        g_filename = '../manuscript/Figures/lagged_edge_analysis.pdf'
        if savefig2 is True and ii > 0:
            g.savefig(g_filename, fmt='pdf')
        elif savefig2 == 'show' and ii > 0:
            plt.show()

    """
    ####################################################################################################################
    ####################################################################################################################
    Specific network figure
    ####################################################################################################################
    ####################################################################################################################
    """

    # Look at one specific network. Yeast12 has the most increase in AUROC)
    gs = '../data/gnw_insilico/network_data/Yeast/Yeast-12_goldstandard.tsv'
    true_edges, edge_df, _, dg, ee = get_network_data(gs, gs.replace('goldstandard', 'timeseries'))
    data = get_experiment_list(gs.replace('goldstandard', 'timeseries'), 21, 10)

    ranks = summary['RF']['Yeast'][12]['rank'].iloc[:, [0, 1, 6]]
    te_ranks = ranks[ranks.index.isin(true_edges)]
    te_p = te_ranks[te_ranks['Base_RF'] > te_ranks['RF-ml_4']]
    te_d = te_ranks[~(te_ranks['Base_RF'] > te_ranks['RF-ml_4'])]
    fe_ranks = ranks[~ranks.index.isin(true_edges)]

    fig = plt.figure(figsize=(10, 9))
    left_int = 4
    right_int = 3
    cols = left_int + right_int
    ax1 = plt.subplot2grid((2, cols), (0, 0), rowspan=2, colspan=left_int)
    ax2 = plt.subplot2grid((2, cols), (0, left_int), colspan=right_int)
    ax3 = plt.subplot2grid((2, cols), (1, left_int), colspan=right_int)
    lw = 3

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
        patch = patches.PathPatch(path, facecolor='none', lw=2, ec='0.75')
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
        edge_str = edge[0] + u"\u2192" + edge[1] + " Lag = " + str(int(lag))
        ax1.text(-0.02, start, edge_str, ha='right', va='center')

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
        if edge == ('G2', 'G1'):
            ec = paired[3]
        patch = patches.PathPatch(path, facecolor='none', lw=lw, ec=ec)
        ax1.add_patch(patch)
        edge_str = edge[0] + u"\u2192" + edge[1] + " Lag = " + str(int(lag))
        ax1.text(-0.02, start, edge_str, ha='right', va='center')

    ax1.set_ylim([91, -1])
    ax1.plot([0, 0], [-1, 91], '-', c='k', lw=1)
    ax1.yaxis.tick_right()
    ax1.spines['top'].set_bounds(0, 1)
    ax1.spines['bottom'].set_bounds(0, 1)
    ax1.spines['left'].set_position(('data', 0))
    ax1.set_xlim([-0.5, 1])
    ax1.set_ylabel('Rank', rotation=0, ha='center')
    ax1.yaxis.set_label_coords(1.005, 1.005)
    ax1.tick_params(axis='x', which='both', bottom='off', top='off')
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['RF', 'SWING\nRF'])

    # In experiment 4, G2 and nothing else upstream is perturbed, so the relation between G2, G1 is clearer
    # The apparent lag is 2, so shift it that much
    # Plot unlagged time series and correlation
    ax2.plot(data[4].index.values, data[4]['G2'].values, '.-', label='G2', c=paired[3], lw=3, ms=15)
    ax2.plot(data[4].index.values, data[4]['G1'].values, '.--', dashes=(5, 5),
             label='G1', c='k', lw=3, ms=15, mfc='w', mew=2)
    ax2.plot(data[4].index.values[:-2], data[4]['G1'].values[2:], '.-', label='G1-shifted', c='k', lw=3, ms=15)
    arrow_nums = [4, 7, 16, 18]
    pad = 15
    for idx, point in enumerate(data[4].index.values):
        if idx > 1 and idx in arrow_nums:
            ax2.annotate("", xy=(data[4].index.values[idx-2]+pad, data[4]['G1'].values[idx]),
                         xytext=(point-pad, data[4]['G1'].values[idx]),
                         arrowprops=dict(facecolor='k', headlength=8, headwidth=7, width=1.45))
    ax2.set_yticks(np.arange(0.0, 0.6, 0.1))
    ax2.set_xticks(np.arange(0.0, 1200, 200))
    ax2.set_ylabel('Normalized expression')
    ax2.set_xlabel('Time')
    ax2.legend(loc='best')

    reg_fit = linregress(data[4]['G2'], data[4]['G1'])
    shift_fit = linregress(data[4]['G2'].values[:-2], data[4]['G1'].values[2:])

    ax3.plot(data[4]['G2'].values, data[4]['G1'].values, '.', c='k', mfc='w', mew=2, ms=15, label='$\mathregular{G_1}$')
    ax3.plot(data[4]['G2'].values[:-2], data[4]['G1'].values[2:], '.', c='k', ms=15,
            label='$\mathregular{G_{1} shifted}$')
    xvals = np.array([np.min(data[4]['G2'].values), np.max(data[4]['G2'].values)])
    ax3.plot(xvals, xvals * reg_fit.slope + reg_fit.intercept,
            c='k', ls='--', dashes=(5, 5), lw=1, zorder=0, label=('$\mathregular{R^2 = %0.3f}$' % reg_fit.rvalue))
    xvals = np.array([np.min(data[4]['G2'].values[:-2]), np.max(data[4]['G2'].values[:-2])])
    ax3.plot(xvals, xvals * shift_fit.slope + shift_fit.intercept,
            c='k', zorder=1, label=('$\mathregular{R^2 = %0.3f}$' % shift_fit.rvalue))
    ax3.set_yticks(np.arange(0.0, 0.4, 0.1))
    ax3.set_xticks(np.arange(0.05, 0.5, 0.1))
    ax3.set_xlabel('G2 normalized expression')
    ax3.set_ylabel('Normalized expression')
    ax3.legend(loc='best', numpoints=1, ncol=2, handletextpad=0.5, handlelength=1, columnspacing=0.5)

    plt.tight_layout()
    if savefigs_group2 is True:
        plt.savefig('../manuscript/Figures/RF_yeast12_edge_promotion_analysis.pdf', fmt='pdf')
    elif savefigs_group2 == 'show':
        plt.show()
    else:
        plt.close()
