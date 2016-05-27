import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

import pandas as pd
import pdb

import Swing.util.lag_identification as lag_id
from Swing.util.Evaluator import Evaluator
import pickle
import numpy as np
import math

from timeit import default_timer as timer


def get_true_edges(gold_filename):
    evaluator = Evaluator(gold_filename, '\t')
    edges = evaluator.gs_flat.tolist()
    return edges, evaluator


def get_true_lags(exp_path, timepoints, perturbs, dset = 'omranian'):
    exp_inter = lag_id.get_experiment_list(exp_path, timepoints = timepoints, perturbs=perturbs)
    
    start = timer()
    exp_xcor = lag_id.xcorr_experiments(exp_inter, 1)
    end = timer()
    print(end - start)      

    print("Calculated x correlation matrices")
    gene_list = list(exp_inter[0].columns.values)
    print("Calculated edge lag from x correlation matrices")
    if 'omranian' in dset:
        signed_edge_list = pd.read_csv('../data/invitro/omranian_signed_parsed_goldstandard.tsv',sep='\t',header=None)
        goldstandard = '../data/invitro/omranian_parsed_goldstandard.tsv'

    elif 'marbach' in dset:
        signed_edge_list = pd.read_csv('../data/invitro/marbach_signed_parsed_goldstandard.tsv',sep='\t',header=None)
        goldstandard = '../data/invitro/marbach_parsed_goldstandard.tsv'
                
    signed_edge_list.columns=['regulator', 'target', 'signs']
    signed_edge_list['regulator-target'] = tuple(zip(signed_edge_list['regulator'], signed_edge_list['target']))

    lags=lag_id.calc_edge_lag(exp_xcor, gene_list, 0.5, 0.6, timestep=1, signed_edge_list = signed_edge_list, flat = False)

    lags = lags[lags['Parent'] != lags['Child']]
    edge_df = pd.DataFrame(lags['Lag'].values, index=lags['Edge'].values, columns=['Lag'])
    
    true_edges, evaluator = get_true_edges(goldstandard)

    all_lag_list = edge_df['Lag'].tolist()
    edge_df['lag_stderr'] = [np.std(x)/math.sqrt(len(x)) if type(x) is list else 'nan' for x in all_lag_list]
    edge_df['lag_std'] = [np.std(x) for x in all_lag_list]
    edge_df['lag_mean'] = [np.mean(x) for x in all_lag_list]

    
    only_true = edge_df[edge_df.index.isin(true_edges)]
    
    lag_list = only_true['Lag'].tolist()
    only_true['lag_stderr'] = [np.std(x)/math.sqrt(len(x)) if type(x) is list else 'nan' for x in lag_list]
    only_true['lag_std'] = [np.std(x) for x in lag_list]
    only_true['lag_mean'] = [np.mean(x) for x in lag_list]
    return(only_true, edge_df)

def get_mean_list(lag_df):
    return(lag_df[~lag_df['lag_mean'].isnull()]['lag_mean'].values)

def get_kde(mean_list):
    density = gaussian_kde(mean_list)
    density.covariance_factor = lambda : .24
    density._compute_covariance()
    return(density)

def main():
    my_paths = ['../data/invitro/omranian_parsed_timeseries.tsv','../data/invitro/omranian_parsed_heatstress_timeseries.tsv','../data/invitro/omranian_parsed_coldstress_timeseries.tsv', '../data/invitro/omranian_parsed_control_timeseries.tsv']

    xs = np.linspace(0,5,200)
    colors = ['r','g','y','b']
    for idx, path in enumerate(my_paths):
        if idx is 0:
            lag_df = get_true_lags(path, 6, 9)
        else:
            lag_df = get_true_lags(path, 6, 3)
        ms = get_mean_list(lag_df)
        density = get_kde(ms)
        plt.plot(xs,density(xs), lw=2.0, color = colors[idx])


    plt.savefig('multiple_lag_kernel_density.png')

    plt.hist(ms, bins=5)
    plt.savefig('multiple_lag_hist.png')

if __name__ == '__main__':
    main()

  #combine/parse omranian

  # Time Gene Names

  #map omranian to strong regulondb gold standard

  #use swing on gold standard

  #check if edges are lagged in omranian

  #convert network into modules

  #get lagged modules (percentage of each module that is lagged)

  #get functional enrichment analysis of modules


