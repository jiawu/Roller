__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

# todo: Clean this up! Make it into a real module

import os, sys, itertools
import networkx as nx
import pandas as pd
from statsmodels.tsa.stattools import ccf
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
from scipy.stats import pearsonr
from sklearn.metrics.cluster import entropy
from sklearn.metrics.cluster import expected_mutual_information
from sklearn.utils.validation import check_array
from math import log
from Swing.util.Evaluator import Evaluator

def get_experiment_list(filename, timepoints=None, perturbs=None, time_col='Time'):
    # load files
    timecourse = pd.read_csv(filename, sep="\t")
    if timepoints is None:
        timepoints = len(set(timecourse[time_col]))

    if perturbs is None:
        perturbs = len(timecourse)/timepoints

    if perturbs.is_integer():
        perturbs = int(perturbs)
    else:
        raise ValueError("Uneven number of timepoints between perturbation experiments")

    # divide into list of dataframes
    experiments = []
    for i in range(0, timepoints * perturbs - timepoints + 1, timepoints):
        experiments.append(timecourse.ix[i:i + timepoints - 1])
    # reformat
    for idx, exp in enumerate(experiments):
        exp = exp.set_index(time_col)
        experiments[idx] = exp
    return experiments


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


def get_xcorr_indices(diff_ts, lag, tolerance):
    pair_list = []
    # get all pairs
    targets = np.array(np.where((diff_ts >= lag-tolerance ) & (diff_ts <= lag+tolerance)))
    n_ind = targets.shape[1]
    pair_list = [tuple(targets[:,x]) for x in range(n_ind)]
    # only keep tuples where the parent index is greater than the child
    if lag != 0:
        pair_list = [ x for x in pair_list if x[1] < x[2]]
    p_pair_list = [(x[0],x[1]) for x in pair_list]
    c_pair_list = [(x[0],x[2]) for x in pair_list]

    return(p_pair_list,c_pair_list)

def get_pairwise_xcorr(parent,child,experiment,time_map,lag,tolerance,rc):      
    ts_shape = time_map.shape[1]-1
    ts = time_map.iloc[:,:ts_shape]
    ts = ts.values

    all_ps_values = np.zeros(rc)
    all_cs_values = np.zeros(rc)

    # make an array of differences

    diff_ts = np.abs(ts[:,:,None] - ts[:,None,:])
    # get all indices with the same difference
    ps_values = np.zeros(rc)
    cs_values = np.zeros(rc)
    ps = [x[parent].values for x in experiment]
    cs = [x[child].values for x in experiment]
    all_ps_values = np.vstack(ps)
    all_cs_values = np.vstack(cs)

    p_idx,c_idx = get_xcorr_indices(diff_ts, lag, tolerance)
    

    ps_values = [all_ps_values[x] for x in p_idx]
    cs_values = [all_cs_values[x] for x in c_idx]

    rsq, pval = pearsonr(ps_values,cs_values)
    #c_xy = np.histogram2d(ps_values, cs_values, 10)[0]
    #n_samples = len(ps_values)
    #mi = my_mutual_info_score(n_samples, x_val = ps_values, y_val = cs_values, labels_true=None, labels_pred=None, contingency=c_xy)
    #print(mi, parent, child, lag)

    return(rsq,pval)

def my_mutual_info_score(n_samples, x_val, y_val, labels_true, labels_pred, contingency=None):
    """Mutual Information between two clusterings.
    The Mutual Information is a measure of the similarity between two labels of
    the same data. Where :math:`P(i)` is the probability of a random sample
    occurring in cluster :math:`U_i` and :math:`P'(j)` is the probability of a
    random sample occurring in cluster :math:`V_j`, the Mutual Information
    between clusterings :math:`U` and :math:`V` is given as:
    .. math::
        MI(U,V)=\sum_{i=1}^R \sum_{j=1}^C P(i,j)\log\\frac{P(i,j)}{P(i)P'(j)}
    This is equal to the Kullback-Leibler divergence of the joint distribution
    with the product distribution of the marginals.
    This metric is independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score value in any way.
    This metric is furthermore symmetric: switching ``label_true`` with
    ``label_pred`` will return the same score value. This can be useful to
    measure the agreement of two independent label assignments strategies
    on the same dataset when the real ground truth is not known.
    Read more in the :ref:`User Guide <mutual_info_score>`.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        A clustering of the data into disjoint subsets.
    labels_pred : array, shape = [n_samples]
        A clustering of the data into disjoint subsets.
    contingency : {None, array, sparse matrix},
                  shape = [n_classes_true, n_classes_pred]
        A contingency matrix given by the :func:`contingency_matrix` function.
        If value is ``None``, it will be computed, otherwise the given value is
        used, with ``labels_true`` and ``labels_pred`` ignored.
    Returns
    -------
    mi : float
       Mutual information, a non-negative value
    See also
    --------
    adjusted_mutual_info_score: Adjusted against chance Mutual Information
    normalized_mutual_info_score: Normalized Mutual Information
    """
    if contingency is None:
        labels_true, labels_pred = check_clusterings(labels_true, labels_pred)
        contingency = contingency_matrix(labels_true, labels_pred, sparse=True)
    else:
        contingency = check_array(contingency,
                                  accept_sparse=['csr', 'csc', 'coo'],
                                  dtype=[int, np.int32, np.int64])

    if isinstance(contingency, np.ndarray):
        # For an array
        nzx, nzy = np.nonzero(contingency)
        nz_val = contingency[nzx, nzy]
    elif sp.issparse(contingency):
        # For a sparse matrix
        nzx, nzy, nz_val = sp.find(contingency)
    else:
        raise ValueError("Unsupported type for 'contingency': %s" %
                         type(contingency))

    contingency_sum = contingency.sum()
    pi = np.ravel(contingency.sum(axis=1))
    pj = np.ravel(contingency.sum(axis=0))
    log_contingency_nm = np.log(nz_val)
    contingency_nm = nz_val / contingency_sum
    # Don't need to calculate the full outer product, just for non-zeroes
    outer = pi.take(nzx) * pj.take(nzy)
    log_outer = -np.log(outer) + log(pi.sum()) + log(pj.sum())
    mi = (contingency_nm * (log_contingency_nm - log(contingency_sum)) +
          contingency_nm * log_outer)
    mi = mi.sum()
    emi = expected_mutual_information(contingency, n_samples)
    # Calculate entropy for each labeling
    h_true, h_pred = entropy(x_val), entropy(y_val)
    ami = (mi - emi) / (max(h_true, h_pred) - emi)
    return ami

def calc_edge_lag2(experiments,genes, signed_edge_list=None, tolerance = 8, rc = (23,6), mode=None):
    
    # load the interval file
    edges = signed_edge_list['regulator-target']
    

    #initialize dataframe to return
    col, row = np.meshgrid(range(len(genes)), range(len(genes)))
    edge_lag = pd.DataFrame()
    edge_lag['parent'] = np.array(genes)[row.flatten()]
    edge_lag['child'] = np.array(genes)[col.flatten()]
    edge_lag['Edge'] = list(zip(edge_lag['parent'], edge_lag['child']))

    lag_results = []
    if mode is 'marbach':
        time_map = pd.read_csv('../../data/invitro/marbach_timesteps.tsv', sep='\t')
        rc = (23,6)
        lags = [0,5,10,20,30,40]
        tolerance = 3
    else:
        time_map = pd.read_csv('../../data/invitro/omranian_timesteps.tsv', sep='\t')
        lags = [0,10,20,30,60,90]

    time_steps = time_map['Timestep'].tolist()

    for edge in edges:
        # Ignore self edges
        if edge[0] == edge[1]:
            continue

        tolerance = 8
        c_list = []
        for lag in lags:
            r,p = get_pairwise_xcorr(edge[0],edge[1],experiments,time_map,lag,tolerance,rc)
            c_list.append((lag,r,p))

        sign = signed_edge_list[signed_edge_list['regulator-target'] == edge]['signs'].tolist()

        best_lag = min(c_list, key = lambda x: x[2])
        if best_lag[2] > 0.05/len(edges):
            true_lag = np.nan
        else:
            true_lag = best_lag[0]
        lag_results.append({'Edge':edge, 'Lag':true_lag, 'Sign': sign, 'Lag_list': c_list})

    lag_results = pd.DataFrame(lag_results)

    edge_lag = pd.merge(edge_lag, lag_results, how='outer', on='Edge')

    lag_results['parent'] = [x[0] for x in lag_results['Edge'].tolist()]
    lag_results['child'] = [x[1] for x in lag_results['Edge'].tolist()]

    return(lag_results, edge_lag)

 
def calc_edge_lag(xcorr, genes, sc_frac=0.1, min_ccf=0.5, timestep=1, signed_edge_list=None, flat=True, return_raw=False):
    """

    :param xcorr: 4d array
        4 axes in order: experiments, parent, child, time
    :param genes: list
    :param sc_frac: float
        related filtering. see filter_ccfs
    :param min_ccf: float
        minimum cross correlation needed to call a lag
    :param timestep: int
    :param signed: dataframe
        can be a list of signed edges or none (default)
        maximize either negative or positive correlation depending on prior information
    :param flat: boolean
        true: return the mean lag for each edge
        false: return the list of all lags (for each exp) for each edge
    :return:
    """
    e, p, c, t = xcorr.shape
    if signed_edge_list is not None:
        edges = signed_edge_list['regulator-target']
    else:
        edges = itertools.product(genes, genes)
    lag_estimate = np.zeros((p,c))
    sc_thresh = sc_frac * t

    #initialize dataframe to return
    col, row = np.meshgrid(range(len(genes)), range(len(genes)))
    edge_lag = pd.DataFrame()
    edge_lag['Parent'] = np.array(genes)[row.flatten()]
    edge_lag['Child'] = np.array(genes)[col.flatten()]
    edge_lag['Edge'] = list(zip(edge_lag['Parent'], edge_lag['Child']))

    lag_results = []

    for edge in edges:
        # Ignore self edges
        if edge[0] == edge[1]:
            continue
        p_idx = genes.index(edge[0])
        c_idx = genes.index(edge[1])
        if signed_edge_list is not None:
            sign = signed_edge_list[signed_edge_list['regulator-target'] == edge]['signs'].tolist()[0]

        # The ccf keeps the parent static and moves the child. Therefore the reversed xcorr would show the true lag
        reverse = xcorr[:, c_idx, p_idx]
        filtered = filter_ccfs(reverse, sc_thresh, min_ccf)
        if filtered.shape[0] > 0:
            # f, axarr = plt.subplots(1,2)
            # axarr[0].plot(reverse.T)
            # axarr[1].plot(filtered.T)
            # plt.show()

            # default setting
            if flat:
                if signed_edge_list is None:
                    lag_estimate[p_idx, c_idx] = float(np.mean(np.argmax(np.abs(filtered), axis=1)))*timestep
                elif sign == '+':
                    lag_estimate[p_idx, c_idx] = float(np.mean(np.argmax(filtered, axis=1)))*timestep
                elif sign == '-':
                    lag_estimate[p_idx, c_idx] = float(np.mean(np.argmin(filtered, axis=1)))*timestep
                elif sign == '+-':
                    lag_estimate[p_idx, c_idx] = float(np.mean(np.argmax(np.abs(filtered), axis=1)))*timestep
                edge_lag['Lag'] = lag_estimate.flatten()

            elif not flat:
                
                if sign == '+':
                    lag = [float(x) for x in np.argmax(filtered, axis=1)]*timestep
                elif sign == '-':
                    lag = [float(x) for x in np.argmin(filtered, axis=1)]*timestep
                elif sign == '+-':
                    lag = [float(x) for x in np.argmax(np.abs(filtered), axis=1)]*timestep
                lag_results.append({'Edge':edge, 'Lag':lag, 'Raw_CCF': lag})

    if not flat:
        lag_results = pd.DataFrame(lag_results)
        edge_lag = pd.merge(edge_lag, lag_results, how='outer', on='Edge')
            # print(edge, np.argmax(filtered, axis=0), np.mean(np.argmax(filtered, axis=0)))

    # If all cross correlations were filtered out, no estimate will be made
    # Assign all lags to 0
    if 'Lag' not in edge_lag.columns:
        edge_lag['Lag'] = 0

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
    
    ### Do not cross correlate with a lag greater than 1/2 of the dataset when the timeseries is short
    ### throw out these cross correlations in filtered time-series
    max_lag = ccfs.shape[1]

    # if max_lag < 10:
    #     max_lag = int(np.ceil(ccfs.shape[1]/2.0))

    filtered_ccf = ccfs[(np.sum(signchange, axis=1) <= sc_thresh) & (np.max(np.abs(ccfs), axis=1) > min_ccf),
                        :max_lag + 1]

    return filtered_ccf


