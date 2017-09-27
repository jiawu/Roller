import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag, round_to
from Swing.util.Evaluator import Evaluator

sampling = [10, 30, 50, 100, 200, 333, 500]
models = ['Ecoli', 'Yeast']
nodes = 10

for s in sampling:
    if s == 'dream4':
        srate = 1
    else:
        srate = s
    print('Sampling rate: '+str(srate))
    lags = []
    thresh = 0.3
    for m in models:
        nets = 20
        data_folder = "../../data/gnw_insilico/high_sampling/%s%i/" % (m, nodes)
        insilico_dict = {ii: {} for ii in range(1, nets + 1)}
        for net in insilico_dict:
            data_file = data_folder + "%s-%i_%s_timeseries.tsv" % (m, net, str(s))
            gold_file = data_folder + "%s-%i_goldstandard.tsv" % (m, net)
            perturb_file = data_folder + "%s-%i_timeseries_perturbations.tsv" % (m, net)

            # Get the true edges
            evaluator = Evaluator(gold_file, '\t')
            true_edges = evaluator.gs_flat.tolist()

            # Calculate the xcorr for each gene pair
            df = pd.read_csv(data_file, sep="\t")
            gene_list = df.columns.values[1:].tolist()
            experiment_list = get_experiment_list(data_file)
            xcorr_array = xcorr_experiments(experiment_list)
            edge_lags = calc_edge_lag(xcorr_array, gene_list, 0.1, 0.3, timestep=srate)
            true_lags = edge_lags[edge_lags['Edge'].isin(true_edges)]
            lags += true_lags['Lag'].values.tolist()

    lags = np.array(lags)
    lags = [round_to(lag, srate) for lag in lags]
    z = np.nan_to_num(np.array(lags))
    data = pd.DataFrame(np.asarray(list(Counter(z).items())))
    data.sort_values(0, inplace=True)
    # plt.plot(data[0].values, data[1]/len(lags), 'o-', lw=5, ms=10, markerfacecolor='w', mew=3, mec='b')
    plt.figure(figsize=(10, 8))
    bar_width = 0.95
    colors = ['grey'] + ['k']*(len(data)-1)
    plt.bar(range(len(data[0].values)), data[1]/len(lags), width=bar_width, color=colors)
    plt.xlabel(r'Estimated t of true interaction(min)', fontsize=16, weight='bold')
    plt.ylabel('% of true interactions', fontsize=16, weight='bold')
    # labels = [ii if ii==0 else str(data[0].values.astype(int)[ii-1])+"-"+str(val) for ii, val in enumerate(data[0].values.astype(int))]
    labels = [val for val in data[0].values.astype(int)]
    labels[0] = 'N/A'
    plt.xticks(np.arange(len(data[0]))+bar_width/2, labels, rotation=45)
    plt.yticks(np.arange(0, .9, 0.1))
    plt.xlim([0, len(labels)-1])
    plt.ylim([0, .8])
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()

    plt.savefig('lag_distribution_{}.pdf'.format(int(srate)), fmt='pdf')