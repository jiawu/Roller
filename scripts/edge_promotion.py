__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import os
import pandas as pd
import numpy as np
import networkx as nx
from Swing.util.Evaluator import Evaluator
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag
from nxpd import draw
from nxpd import nxpdParams
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib.patches import Polygon
import brewer2mpl
import seaborn as sns
from scipy.stats import fisher_exact, linregress, ttest_rel, mannwhitneyu, ttest_ind


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
    elif (p < 0.001):
        return "***"
    elif (p < 0.01):
        return "**"
    elif (p < 0.05):
        return "*"
    else:
        return "-"

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
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'text.usetex': False,
        'figure.figsize': [2.5, 4.5]
    }
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors
    rcParams.update(params)
    m = ['Dionesus', 'RF']
    rd = {'Dionesus': 'D', 'RF': 'RF'}
    mod = ['Ecoli', 'Yeast']
    try:
        print("Loading_pickle")
        summary = pd.read_pickle('./param_sweep_summary.pickle')
    except FileNotFoundError:
        print("Pickle doesn't exist. Summarizing data now")
        summary = make_dictionary(m, rd, mod)
        pd.to_pickle(summary, './param_sweep_summary.pickle')

    scores = ['aupr', 'auroc']
    bprops = dict(linewidth=3)
    mprops = dict(linewidth=3, color='r')
    wprops = dict(linewidth=3, linestyle='--')
    # for method in m:
    #     diffs = pd.DataFrame()
    #     for model in mod:
    #         for score in scores:
    #             #             print('\n', method, model, score)
    #             avg_diff = pd.DataFrame(np.mean(summary[method][model][(score + '_diff')], axis=0),
    #                                     columns=[model + "_" + score])
    #             avg_diff['_'.join([model, score, 'pval'])] = [ttest_rel(summary[method][model][score].iloc[:, 0],
    #                                                                     summary[method][model][score].iloc[:, ii]).pvalue
    #                                                           for ii in
    #                                                           range(len(summary[method][model][score].columns))]
    #             x_array = np.array([[0] * len(summary[method][model][score]), [1] * len(summary[method][model][score])])
    #             f = plt.figure(figsize=(10, 10))
    #             for ii, idx in enumerate(avg_diff.index):
    #                 if ii == 0:
    #                     continue
    #                 y_array = np.array(
    #                     [summary[method][model][score].iloc[:, 0], summary[method][model][score].iloc[:, ii]])
    #                 ax = f.add_subplot(3, 3, ii)
    #                 ax.plot(x_array, y_array, '.-', c='k', alpha=0.4, zorder=2)
    #                 if score == 'aupr':
    #                     ax.plot([-0.50, 1.5], [.17, .17], c='g', lw=3, zorder=1)
    #                 else:
    #                     ax.plot([-0.50, 1.5], [0.5, 0.5], c='g', lw=3, zorder=1)
    #                 p_value = avg_diff.iloc[ii, 1]
    #                 s = stars(p_value)
    #                 bp = ax.boxplot(y_array.T, positions=[0, 1])
    #                 y_max = np.max(np.concatenate((y_array[0], y_array[1])))
    #                 y_min = np.min(np.concatenate((y_array[0], y_array[1])))
    #                 ax.annotate("", xy=(0, y_max + .01), xycoords='data',
    #                             xytext=(1, y_max + .01), textcoords='data',
    #                             arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
    #                                             connectionstyle="bar,fraction=0.1"))
    #                 ax.text(0.5, y_max + abs(y_max - y_min) * 0.1, stars(p_value),
    #                         horizontalalignment='center',
    #                         verticalalignment='center')
    #                 #                 for i in range(0, len(bp['boxes'])):
    #                 #                     bp['boxes'][i].set_color(colors[i])
    #                 #                     # we have two whiskers!
    #                 #                     bp['whiskers'][i*2].set_color(colors[i])
    #                 #                     bp['whiskers'][i*2 + 1].set_color(colors[i])
    #                 #                     bp['whiskers'][i*2].set_linewidth(2)
    #                 #                     bp['whiskers'][i*2 + 1].set_linewidth(2)
    #                 #                     # top and bottom fliers
    #                 #                     # (set allows us to set many parameters at once)
    #                 #                     bp['fliers'][i].set(markerfacecolor=colors[i],
    #                 #                                    marker='o', alpha=0.75, markersize=6,
    #                 #                                    markeredgecolor='none')
    #                 #                     bp['fliers'][i].set(markerfacecolor=colors[i],
    #                 #                                    marker='o', alpha=0.75, markersize=6,
    #                 #                                    markeredgecolor='none')
    #                 #                     bp['medians'][i].set_color('black')
    #                 #                     bp['medians'][i].set_linewidth(3)
    #                 #                     # and 4 caps to remove
    #                 #                     for c in bp['caps']:
    #                 #                         c.set_linewidth(0)
    #                 #                 for i in range(len(bp['boxes'])):
    #                 #                     box = bp['boxes'][i]
    #                 #                     box.set_linewidth(0)
    #                 #                     boxX = []
    #                 #                     boxY = []
    #                 #                     for j in range(5):
    #                 #                         boxX.append(box.get_xdata()[j])
    #                 #                         boxY.append(box.get_ydata()[j])
    #                 #                         boxCoords = zip(boxX,boxY)
    #                 #                         print(boxCoords)
    #                 #                         boxPolygon = Polygon(boxCoords, facecolor = colors[i], linewidth=0)
    #                 #                         ax.add_patch(boxPolygon)
    #                 ax.spines['top'].set_visible(False)
    #                 ax.spines['right'].set_visible(False)
    #                 ax.spines['left'].set_visible(False)
    #                 ax.get_xaxis().tick_bottom()
    #                 ax.get_yaxis().tick_left()
    #                 ax.tick_params(axis='x', direction='out')
    #                 ax.tick_params(axis='y', length=0)
    #                 ax.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    #                 ax.set_axisbelow(True)
    #                 ax.set_xlim([-0.5, 1.5])
    #                 #                 ax.set_title(' '.join([idx, model, score]))
    #             plt.tight_layout()
    #             plt.close()
    #             # plt.show()
    #
    #             diffs = pd.concat([diffs, avg_diff], axis=1)

    nl = summary['RF']['te_change'][(summary['RF']['te_change']["Lag"] == 0)]
    l = summary['RF']['te_change'][(summary['RF']['te_change']["Lag"] > 0)]
    for lag in set(summary['RF']['te_change']["Lag"]):
        data = summary['RF']['te_change']['RF-ml_4'][summary['RF']['te_change']['Lag'] == lag]
        if len(data) >= 3:
            vp = plt.violinplot(data, positions=[lag], showmedians=False, showextrema=False)
            for pc in vp['bodies']:
                pc.set_facecolor('c')
                pc.set_edgecolor('w')
                pc.set_alpha(1)
            plt.plot([lag-0.1, lag+0.1], [np.median(data), np.median(data)], c='k')
    plt.xlim([-1, 6])
    plt.ylim([-90, 90])
    plt.show()
    sys.exit()
    distros = [np.mean(l['RF-ml_4'][l['Lag']==lag]) for lag in set(l["Lag"])]
    plt.plot(np.arange(1, 9), distros)
    # plt.boxplot(distros, positions=np.arange(8))
    plt.show()
    sys.exit()
    print(np.median(l.iloc[:, 2:], axis=0))
    print(np.median(l.iloc[:, 2:], axis=0) - np.median(nl.iloc[:, 2:], axis=0))
    cond = 'RF-ml_4'
    plt.plot([0.5, 2.5], [0, 0], c='k', lw=1)
    vp = plt.violinplot([nl[cond].values, l[cond].values], showextrema=False, widths=0.75, points=200)
    for pc in vp['bodies']:
        pc.set_facecolor('c')
        pc.set_edgecolor('w')
    plt.plot([0.75, 1.25], [np.median(nl[cond].values), np.median(nl[cond].values)], c='k', lw=3)
    plt.plot([1.75, 2.25], [np.median(l[cond].values), np.median(l[cond].values)], c='k', lw=3)
    ax = plt.gca()
    ax.set_axis_bgcolor('w')
    plt.show()
    print(mannwhitneyu(nl[cond].values, l[cond].values))

    plt.figure(figsize=(5, 10))
    ax = sns.swarmplot(data=[nl[cond].values, l[cond].values], size=5)
    ax.plot([-0.10, 0.10], [np.median(nl[cond].values), np.median(nl[cond].values)], c='k', lw=3, zorder=10)
    ax.plot([0.90, 1.10], [np.median(l[cond].values), np.median(l[cond].values)], c='k', lw=3, zorder=10)
    ax.set_axis_bgcolor('w')
    plt.show()




