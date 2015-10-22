import sys
sys.path.append('../pipelines')
import pdb
import Pipelines as pl
import run_pipeline_RF_window_scan_10 as myRF
import run_pipeline_lasso_window_scan_10 as myLasso

#principal component analysis of results


#get the edge list of several different kinds of results

#RF, tdRF, RF

import td_wrapper as tdw

RF_ranked_edge_list_full = myRF.main(21, 1, 'test')
RF_ranked_edge_list_crag = myRF.main(15, 1, 'test')

Lasso_ranked_edge_list_full = myLasso.main(21)

td_edge_lists_lag = []

for lag in range(0, 4):
    td_edge_list = []

    for network_index in range(1,6):
        td_edge_lists_lag.append([])
        target_dataset = "/projects/p20519/Roller/data/dream4/insilico_size10_" + str(network_index) + "_timeseries.tsv"
        roc, pr, RF_td_edge_list = tdw.get_td_stats(target_dataset, min_lag = lag)
        td_edge_list.append(RF_td_edge_list)

    td_edge_lists_lag.append(td_edge_list)
pdb.set_trace()
