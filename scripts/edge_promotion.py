__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import pandas as pd
import numpy as np
import networkx as nx
from Swing.util.Evaluator import Evaluator
from nxpd import draw
from nxpd import nxpdParams
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import fisher_exact, linregress

from lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag


def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


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
    return edges


def get_edge_lags(data_filename):
    df = pd.read_csv(data_filename, sep="\t")
    gene_list = df.columns.values[1:].tolist()
    experiment_list = get_experiment_list(data_filename, 21, 10)
    xcorr_array = xcorr_experiments(experiment_list)
    lags = calc_edge_lag(xcorr_array, gene_list, 0.1, 0.5, timestep=1)
    return lags


def get_network_changes(pickle_filename, edge_str='regulator-target',
                        base_str='rank_importance_RF-td_21', shortener_str='rank_importance_'):

    results_df = pd.read_pickle(pickle_filename)
    edges = results_df[edge_str].values
    baseline = results_df[base_str].values

    diff_df = pd.DataFrame()
    diff_df[edge_str] = edges
    diff_df['base_rank'] = baseline
    rank_df = pd.DataFrame()
    rank_df[edge_str] = edges
    rank_df['base_rank'] = baseline
    for column in results_df.columns:
        if column != edge_str and column != base_str:
            short_name = column.replace(shortener_str, "")
            diff_df[short_name] = baseline - results_df[column].values
            rank_df[short_name] = results_df[column].values
    return diff_df, rank_df


def calc_stats(df, edges, lags):
    true_df = df[df['regulator-target'].isin(edges)]
    false_df = df[~df['regulator-target'].isin(edges)]


    t_promoted, t_demoted, t_same = count_change(np.sign(true_df.iloc[:, 2:].values))
    f_promoted, f_demoted, f_same = count_change(np.sign(false_df.iloc[:, 2:].values))


def count_change(x, axis=0, normalized=True):
    frac = 1
    if normalized:
        frac = x.shape[axis]
    # x is 2-d np array of signed values
    p = np.sum(x == 1, axis=axis)/frac
    d = np.sum(x == -1, axis=axis)/frac
    s = np.sum(x == 0, axis=axis)/frac
    return p, d, s

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
    num_nets = 20
    models = ['Ecoli', 'Yeast']
    true_edge_df = pd.DataFrame()
    roc_df = pd.DataFrame()
    pr_df = pd.DataFrame()
    for model in models:
        for net in range(1, num_nets+1):
            gold_file = "../data/gnw_insilico/network_data/%s/%s-%i_goldstandard.tsv" % (model, model, net)
            data_file = "../data/gnw_insilico/network_data/%s/%s-%i_timeseries.tsv" % (model, model, net)
            pickle_file = "%s_net%i_promotion.pkl" % (model, net)
            pd.set_option('display.width', 500)
            true_edges = get_true_edges(gold_file)
            ee = Evaluator(gold_file, sep='\t')
            edge_lags = get_edge_lags(data_file)

            # Remove self edges
            #todo: functionalize this better
            edge_lags = edge_lags[edge_lags['Parent'] != edge_lags['Child']]
            edge_df = pd.DataFrame(edge_lags['Lag'].values, index=edge_lags['Edge'].values, columns=['Lag'])
            change_df, ranks_df = get_network_changes(pickle_file)
            conditions = ranks_df.columns[1:].values
            roc = []
            aupr = []
            for c in conditions:
                ranks_df.sort_values(ranks_df.columns[1], inplace=True)
                roc.append(ee.calc_roc(ranks_df.iloc[:, :2])[2].values[-1])
                aupr.append(ee.calc_pr(ranks_df.iloc[:, :2])[2].values[-1])
                ranks_df.drop(c, axis=1, inplace=True)
            roc_df[model + str(net)] = roc
            pr_df[model + str(net)] = aupr
            full_df = pd.concat([edge_df, change_df.set_index(['regulator-target'])], axis=1, join='inner')
            true_edge_df = true_edge_df.append(full_df[full_df.index.isin(true_edges)])
    roc_df.index = conditions
    roc_df = roc_df.T
    baseline = np.reshape(np.repeat(roc_df.values[:, 0], len(conditions)), (len(roc_df), len(conditions)))
    roc_delta = roc_df-baseline
    pr_df.index = conditions
    pr_df = pr_df.T
    baseline = np.reshape(np.repeat(pr_df.values[:, 0], len(conditions)), (len(pr_df), len(conditions)))
    pr_delta = pr_df - baseline
    roc_mean = np.mean(roc_delta, axis=0)
    roc_std = np.std(roc_delta, axis=0)
    pr_mean = np.mean(pr_delta, axis=0)
    pr_std = np.std(pr_delta, axis=0)
    plt.errorbar(roc_mean, pr_mean, roc_std, pr_std, '.')
    plt.show()
    print(roc_mean)
    sys.exit()
    for ii, c in enumerate(conditions):
        plt.plot(roc_df.iloc[:, ii], pr_df.iloc[:, ii], '.', label=c)
    x = roc_df.values.flatten()
    y = pr_df.values.flatten()
    fit = linregress(x, y)
    # plt.plot(x, y, '.')
    plt.plot(x, fit[0]*x+fit[1], c='0.75', label=('R2 = %0.3f' % (fit[2])))
    plt.xlabel('AUROC')
    plt.ylabel('AUPR')
    plt.legend(loc='best')
    plt.show()
    # figure1 = plt.figure()
    # heatmap = figure1.add_subplot(1, 1, 1)
    # my_axis = heatmap.imshow(roc_df)
    # my_axi = my_axis.get_axes()
    # clean_axis(my_axi)
    # heatmap.set_yticks(np.arange(roc_df.shape[0]))
    # heatmap.yaxis.set_ticks_position('left')
    # # heatmap.set_yticklabels(row_labels)
    # heatmap.set_xticks(np.arange(roc_df.shape[1]))
    # heatmap.xaxis.set_ticks_position('top')
    # # xlabelsL = heatmap.set_xticklabels(col_labels)
    # plt.show()
    sys.exit()

    t_promoted = np.sum(true_edge_df.iloc[:, 2:].values > 0, axis=0)
    t_demoted = np.sum(true_edge_df.iloc[:, 2:].values < 0, axis=0)
    t_same = np.sum(true_edge_df.iloc[:, 2:].values == 0, axis=0)
    print(t_promoted/len(true_edge_df), '\n')

    t_lagged = true_edge_df[true_edge_df['Lag'] > 0]
    l_promoted = np.sum(t_lagged.iloc[:, 2:].values > 0, axis=0)
    l_demoted = np.sum(t_lagged.iloc[:, 2:].values < 0, axis=0)
    l_same = np.sum(t_lagged.iloc[:, 2:].values == 0, axis=0)
    print(l_promoted / len(t_lagged))
    c_table = np.array([[l_promoted, t_promoted-l_promoted], [l_demoted+l_same, t_demoted+t_same-l_demoted-l_same]])
    conditions = true_edge_df.columns[2:].values
    for run in range(c_table.shape[2]):
        print(conditions[run])
        print(fisher_exact(c_table[:, :, run]))

    lag_set = sorted(list(set(true_edge_df['Lag'].values)))
    conditions = true_edge_df.columns[1:].values
    f = plt.figure()
    n_rows, n_cols = calc_subplot_dimensions(len(conditions))
    for ii, run in enumerate(conditions):
        ax = f.add_subplot(n_rows, n_cols, ii+1)
        # plot_data = [true_edge_df[true_edge_df['Lag'] == lag][[run]].values for lag in lag_set]
        ax.plot([0, 90], [0, 0], '-', c='k', lw=3)
        base_zero = true_edge_df[true_edge_df['Lag'] == 0][['base_rank']].values
        new_zero = true_edge_df[true_edge_df['Lag'] == 0][[run]].values
        base = true_edge_df[true_edge_df['Lag'] > 0][['base_rank']].values
        new = true_edge_df[true_edge_df['Lag'] > 0][[run]].values
        ax.plot(base_zero, new_zero, '.', c='b')
        ax.plot(base, new, '.', c='r')
        ax.set_title(run)
        ax.set_ylim([-90, 90])
    plt.tight_layout()
    plt.show()

    # Wilcoxon test for each network to see if it changed? - may not be appropriate
    # Hypergeometric test for promoted edges enriched for lagged edges?
    # Calculate for each experimental condition
    # Get true edge lag





"""
Old plotting code

demoted = new_df[new_df['diff'] <= 0]
promoted = new_df[new_df['diff'] > 0]
p_edges = set(promoted['regulator-target'].values)
d_edges = set(demoted['regulator-target'].values)
print(len(p_edges.intersection(set(true_edges))), len(d_edges.intersection(set(true_edges))))
# print(len(d_edges.intersection(set(true_edges))))
# print(len(true_edges))
for row in demoted.iterrows():
    print(row[1])
    sys.exit()
    c = 'r'
    plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
for row in promoted.iterrows():
    c = 'b'
    plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
plt.xticks([0, 1], ['base_rank', column.replace('rank_importance_', "")])
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
"""






