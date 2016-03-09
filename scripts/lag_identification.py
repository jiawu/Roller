__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

#todo: Clean this up! Make it into a real module

import os, sys
import networkx as nx
import pandas as pd
from scipy import stats
from statsmodels.tsa.stattools import acf, ccf
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
from Swing import Swing
from Swing.util.Evaluator import Evaluator
from collections import Counter
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['font.sans-serif'] = 'Arial'


def get_experiment_list(filename, timepoints=21, perturbs=5):
    # load files
    timecourse = pd.read_csv(filename, sep="\t")
    # divide into list of dataframes
    experiments = []
    for i in range(0,timepoints*perturbs-timepoints+1, timepoints):
        experiments.append(timecourse.ix[i:i+timepoints-1])
    #reformat
    for idx,exp in enumerate(experiments):
        exp = exp.set_index('Time')
        experiments[idx]=exp
    return(experiments)

def xcorr_experiments(experiments, gene_axis=1):
    #### NOTE #####
    #### SCIPY may have a way to do this efficiently for 2d arrays. Take a look when possible ##############
    return np.array([cc_experiment(experiment.values.T) if gene_axis==1 else cc_experiment(experiment.values)
            for experiment in experiments])

def cc_experiment(x):
    """
    For one experiment.
    x should be n rows (genes) by m columns (timepoints)
    """
    ccf_array = np.zeros((x.shape[0], x.shape[0], x.shape[1]))
    for ii, static in enumerate(x):
        for jj, moving in enumerate(x):
            if ii==jj:
                unbiased=True
            else:
                unbiased=False
            ccf_array[ii][jj] = ccf(static, moving, unbiased=unbiased)
    return ccf_array

if __name__ == '__main__':
    data_folder = "../data/dream4/"
    nets = 5
    thresh=0.5
    insilico_dict = {ii:{} for ii in range(1, nets+1)}
    lags = []
    for net in insilico_dict:
        data_file = data_folder + "insilico_size10_%i_timeseries.tsv"%(net)
        gold_file = data_folder + "insilico_size10_%i_goldstandard.tsv"%(net)
        perturb_file = data_folder + "insilico_size10_%i_timeseries_perturbations.tsv"%(net)

        # Calculate the xcorr for each gene pair
        df = pd.read_csv(data_file, sep="\t")
        gene_list = df.columns.values[1:].tolist()
        experiments=get_experiment_list(data_file)
        xcorr_list = xcorr_experiments(experiments)

        # Get the true edges
        evaluator = Evaluator(gold_file, '\t')
        true_edges = evaluator.gs_flat.tolist()
        dg = nx.DiGraph()
        dg.add_edges_from(true_edges)
    #     print(nx.shortest_path(dg, 'G9', 'G7'))

        # Get the true perturbations
        true_perturbs = pd.read_csv(perturb_file, sep="\t")
        edge_max_ccf= pd.DataFrame(np.max(np.max(np.abs(xcorr_list), axis=3), axis=0), index=gene_list, columns=gene_list)
        a, b = np.meshgrid(range(len(gene_list)), range(len(gene_list)))
        ll = pd.DataFrame()
        gene_list = np.array(gene_list)
        ll['Parent'] = gene_list[a.flatten()]
        ll['Child'] = gene_list[b.flatten()]
        ll['Max_ccf'] = edge_max_ccf.values.flatten()
        ll['is_lag'] = ll['Max_ccf']>=thresh
        ll['Edge'] = list(zip(ll['Parent'], ll['Child']))
        ll['True_Edge'] = ll['Edge'].isin(true_edges)
    #     print(ll[ll['is_lag']==False], '\n')
        false_neg = ll[(ll['is_lag']==False) & (ll['True_Edge']==True)]['Edge'].values
        genes = list(gene_list)
        x = np.ma.masked_where((np.max(np.abs(xcorr_list), axis=3)<thresh),
                               np.argmax(np.abs(xcorr_list), axis=3))
        for edge in true_edges:
            p_idx = list(gene_list).index(edge[0])
            c_idx = list(gene_list).index(edge[1])
            print(edge, np.ceil(float(np.ma.mean(x[:, c_idx,p_idx]))))
            lags.append(np.ceil(float(np.ma.mean(x[:, c_idx,p_idx]))))
        for falsey in false_neg:
            p_idx = genes.index(falsey[0])
            c_idx = genes.index(falsey[1])
            forward = xcorr_list[:, p_idx, c_idx].T
            reverse = xcorr_list[:, c_idx, p_idx].T
            diff = (forward-reverse)
            diff = diff/np.max(np.abs(diff))/np.array(range(1,len(diff)+1))[:, None]
            max_abs = np.max(np.abs(reverse))
    #         plt.figure()
    #         plt.plot(reverse, 'o-')
    #         plt.title(falsey)
    z = np.nan_to_num(np.array(lags))
    data = pd.DataFrame(np.asarray(list(Counter(z).items())))
    data.sort_values(0, inplace=True)
    plt.plot(data[0].values, data[1]/len(lags), 'o-', lw=5, ms=10, markerfacecolor='w', mew=3, mec='b')
    plt.xlabel('Peak Lag', fontsize=20)
    plt.ylabel('% of True Edges', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.show()
    # plt.savefig('../manuscript/Figures/true_lag_DREAM4.pdf', fmt='pdf')