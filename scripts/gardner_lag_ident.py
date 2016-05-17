import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from nxpd import draw
from collections import Counter
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag, round_to, cc_experiment
from Swing.util.Evaluator import Evaluator

data_folder = "../data/invitro/"
nets = 1
thresh = 0.3
insilico_dict = {ii: {} for ii in range(1, nets + 1)}
lags = []
for net in insilico_dict:
    data_file = data_folder + "gardner_timeseries.tsv"
    gold_file = data_folder + "gardner_goldstandard.tsv"
    # perturb_file = data_folder + "Ecoli-%i_timeseries_perturbations.tsv" % (net)

    # Get the true edges
    evaluator = Evaluator(gold_file, '\t')
    true_edges = evaluator.gs_flat.tolist()
    dg = nx.DiGraph()
    dg.add_edges_from(true_edges)

    # Calculate the xcorr for each gene pair
    df = pd.read_csv(data_file, sep="\t")
    gene_list = df.columns.values[1:].tolist()
    experiment_list = get_experiment_list(data_file, 14, 1)
    ccf = cc_experiment(experiment_list[0].values.T)
    nodes = df.columns.values[1:]
    f = plt.figure()
    for ii, parent in enumerate(ccf):
        for jj, child in enumerate(parent):
            ax = f.add_subplot(ccf.shape[0], ccf.shape[1], ii*8+jj+1)
            if (nodes[ii], nodes[jj]) in true_edges:
                color = 'b'
            elif (nodes[jj], nodes[ii]) in true_edges:
                color = 'r'
            else:
                color = 'k'
            ax.plot(child, '.-', c=color)
            if ii == 0:
                ax.set_title(nodes[jj])
            if jj == 0:
                ax.set_ylabel(nodes[ii])

    xcorr_array = xcorr_experiments(experiment_list)
    edge_lags = calc_edge_lag(xcorr_array, gene_list, 1, 0.3, timestep=1)
    true_lags = edge_lags[edge_lags['Edge'].isin(true_edges)]
    print(true_lags)
    lags += true_lags['Lag'].values.tolist()


# plt.plot(df.Time, df.values[:, 1:], '.-')
plt.show()
sys.exit()
lags = [round_to(lag, 50) for lag in lags]
z = np.nan_to_num(np.array(lags))
data = pd.DataFrame(np.asarray(list(Counter(z).items())))
data.sort_values(0, inplace=True)
# plt.plot(data[0].values, data[1]/len(lags), 'o-', lw=5, ms=10, markerfacecolor='w', mew=3, mec='b')
plt.figure(figsize=(6, 5))
bar_width = 0.95
plt.bar(range(len(data[0].values)), data[1]/len(lags), width=bar_width, color='k')
# plt.hist(lags, bins=71, cumulative=True, normed=1)
plt.xlabel('Apparent Lag (min)', fontsize=20)
plt.ylabel('% of Edges', fontsize=20)
# labels = [ii if ii==0 else str(data[0].values.astype(int)[ii-1])+"-"+str(val) for ii, val in enumerate(data[0].values.astype(int))]
labels = [val for val in data[0].values.astype(int)]
plt.xticks(np.arange(len(data[0]))+bar_width/2, labels)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.show()
# plt.savefig('../manuscript/Figures/true_lag_ecoli10_0.1_0.3.pdf', fmt='pdf')