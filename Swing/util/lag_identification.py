__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

# todo: Clean this up! Make it into a real module

import os, sys, itertools
import networkx as nx
import pandas as pd
from statsmodels.tsa.stattools import ccf
import matplotlib.pyplot as plt
import numpy as np
from Swing.util.Evaluator import Evaluator
from collections import Counter
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'


def get_experiment_list(filename, timepoints=21, perturbs=5):
    # load files
    timecourse = pd.read_csv(filename, sep="\t")
    # divide into list of dataframes
    experiments = []
    for i in range(0, timepoints * perturbs - timepoints + 1, timepoints):
        experiments.append(timecourse.ix[i:i + timepoints - 1])
    # reformat
    for idx, exp in enumerate(experiments):
        exp = exp.set_index('Time')
        experiments[idx] = exp
    return (experiments)


def xcorr_experiments(experiments, gene_axis=1):
    """
    Cross correlate the g
    :param experiments: list
        list of dataframes
    :param gene_axis: int
        axis corresponding to each gene. 0 for rows, 1 for columns
    :return:
    """
    return np.array([cc_experiment(experiment.values.T) if gene_axis == 1 else cc_experiment(experiment.values)
                     for experiment in experiments])


def cc_experiment(experiment):
    """
    For one experiment.
    x should be n rows (genes) by m columns (timepoints)
    :param experiment:
    :return:
    """
    ccf_array = np.zeros((experiment.shape[0], experiment.shape[0], experiment.shape[1]))
    for ii, static in enumerate(experiment):
        for jj, moving in enumerate(experiment):
            if ii == jj:
                unbiased = True
            else:
                unbiased = False
            ccf_array[ii][jj] = ccf(static, moving, unbiased=unbiased)
    return ccf_array


def calc_edge_lag(xcorr, genes, sc_frac=0.1, min_ccf=0.5, timestep=1):
    """

    :param xcorr: 4d array
        4 axes in order: experiments, parent, child, time
    :param genes: list
    :return:
    """
    e, p, c, t = xcorr.shape
    edges = itertools.product(genes, genes)
    lag_estimate = np.zeros((p,c))
    sc_thresh = sc_frac * t

    for edge in edges:
        # Ignore self edges
        if edge[0] == edge[1]:
            continue
        p_idx = genes.index(edge[0])
        c_idx = genes.index(edge[1])

        # The ccf keeps the parent static and moves the child. Therefore the reversed xcorr would show the true lag
        reverse = xcorr[:, c_idx, p_idx]
        filtered = filter_ccfs(reverse, sc_thresh, min_ccf)
        if filtered.shape[0] > 0:
            # f, axarr = plt.subplots(1,2)
            # axarr[0].plot(reverse.T)
            # axarr[1].plot(filtered.T)
            # plt.show()
            lag_estimate[p_idx, c_idx] = np.ceil(float(np.mean(np.argmax(np.abs(filtered), axis=1))))*timestep
            # print(edge, np.argmax(filtered, axis=0), np.mean(np.argmax(filtered, axis=0)))
    col, row = np.meshgrid(range(len(genes)), range(len(genes)))
    edge_lag = pd.DataFrame()
    edge_lag['Parent'] = np.array(genes)[row.flatten()]
    edge_lag['Child'] = np.array(genes)[col.flatten()]
    edge_lag['Lag'] = lag_estimate.flatten()
    edge_lag['Edge'] = list(zip(edge_lag['Parent'], edge_lag['Child']))
    return edge_lag


def round_to(x, base, type='ceil'):
    if type == 'ceil':
        r = np.ceil(x/base)*base
    elif type == 'round':
        r = round(x/base, 0)*base
    elif type == 'floor':
        r = np.floor(x/base)*base

    return r


def filter_ccfs(ccfs, sc_thresh, min_ccf):
    """
    Remove noisy ccfs from irrelevant experiments
    :param ccfs: 2d array
    :param sc_thresh: int
        number of sign changes expected
    :param min_ccf: float
        cutoff value for a ccf to be above the noise threshold
    :return:
    """
    if sc_thresh is None:
        sc_thresh = np.inf
    asign = np.sign(ccfs)
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[:, 0] = 0
    # (np.sum(signchange, axis=1) <= sc_thresh) &
    filtered_ccf = ccfs[(np.sum(signchange, axis=1) <= sc_thresh) & (np.max(np.abs(ccfs), axis=1) > min_ccf), :]
    return filtered_ccf

if __name__ == '__main__':
    data_folder = "../data/gnw_insilico/network_data/Yeast100/"
    nets = 20
    thresh = 0.3
    insilico_dict = {ii: {} for ii in range(1, nets + 1)}
    lags = []
    for net in insilico_dict:
        data_file = data_folder + "Yeast100-%i_dream4_timeseries.tsv" % (net)
        gold_file = data_folder + "Yeast100-%i_goldstandard.tsv" % (net)
        perturb_file = data_folder + "Yeast100-%i_dream4__timeseries_perturbations.tsv" % (net)

        # Get the true edges
        evaluator = Evaluator(gold_file, '\t')
        true_edges = evaluator.gs_flat.tolist()

        # Calculate the xcorr for each gene pair
        df = pd.read_csv(data_file, sep="\t")
        gene_list = df.columns.values[1:].tolist()
        experiment_list = get_experiment_list(data_file, 21, 10)
        xcorr_array = xcorr_experiments(experiment_list)
        edge_lags = calc_edge_lag(xcorr_array, gene_list, 0.1, 0.3, timestep=50)

        true_lags = edge_lags[edge_lags['Edge'].isin(true_edges)]
        lags += true_lags['Lag'].values.tolist()

    lags = [round_to(lag, 50) for lag in lags]
    z = np.nan_to_num(np.array(lags))
    data = pd.DataFrame(np.asarray(list(Counter(z).items())))
    data.sort_values(0, inplace=True)
    # plt.plot(data[0].values, data[1]/len(lags), 'o-', lw=5, ms=10, markerfacecolor='w', mew=3, mec='b')
    plt.bar(range(len(data[0].values)), data[1]/len(lags), width=0.7, color='k')
    # plt.hist(lags, bins=71, cumulative=True, normed=1)
    plt.xlabel('Apparent Lag', fontsize=20)
    plt.ylabel('% of Edges', fontsize=20)
    # labels = [ii if ii==0 else str(data[0].values.astype(int)[ii-1])+"-"+str(val) for ii, val in enumerate(data[0].values.astype(int))]
    labels = [val for val in data[0].values.astype(int)]
    plt.xticks(np.arange(len(data[0]))+0.35, labels)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tight_layout()
    plt.show()
    # plt.savefig('../manuscript/Figures/true_lag_DREAM4.pdf', fmt='pdf')
