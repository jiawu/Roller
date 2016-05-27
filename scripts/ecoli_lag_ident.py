import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Swing.util.lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag, round_to
from Swing.util.Evaluator import Evaluator

models = ['Ecoli', 'Yeast']
lags = []
thresh = 0.3
for m in models:
    nets = 20
    data_folder = "../data/gnw_insilico/network_data/%s/" %(m)
    insilico_dict = {ii: {} for ii in range(1, nets + 1)}
    for net in insilico_dict:
        data_file = data_folder + "%s-%i_timeseries.tsv" % (m, net)
        gold_file = data_folder + "%s-%i_goldstandard.tsv" % (m, net)
        perturb_file = data_folder + "%s-%i_timeseries_perturbations.tsv" % (m, net)

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
lags = np.array(lags)
plt.hist(lags[lags==0], color='0.5')
plt.hist(lags[lags != 0], bins=20, color='k')
plt.show()
sys.exit()
lags = [round_to(lag, 50) for lag in lags]
z = np.nan_to_num(np.array(lags))
data = pd.DataFrame(np.asarray(list(Counter(z).items())))
data.sort_values(0, inplace=True)
# plt.plot(data[0].values, data[1]/len(lags), 'o-', lw=5, ms=10, markerfacecolor='w', mew=3, mec='b')
plt.figure(figsize=(5, 4))
bar_width = 0.95
plt.bar(range(len(data[0].values)), data[1]/len(lags), width=bar_width, color='k')
# plt.hist(lags, bins=71, cumulative=True, normed=1)
plt.xlabel(r'Estimated t of true interaction(min)', fontsize=16, weight='bold')
plt.ylabel('% of true interactions', fontsize=16, weight='bold')
# labels = [ii if ii==0 else str(data[0].values.astype(int)[ii-1])+"-"+str(val) for ii, val in enumerate(data[0].values.astype(int))]
labels = [val for val in data[0].values.astype(int)]
plt.xticks(np.arange(len(data[0]))+bar_width/2, labels)
plt.xlim([0, len(labels)-1])
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()

plt.savefig('../manuscript/Figures/true_lag_both10_0.1_0.3.pdf', fmt='pdf')