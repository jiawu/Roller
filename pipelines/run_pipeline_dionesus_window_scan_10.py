import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Swing_old
import uuid
import pickle
import pdb
import pandas as pd
from Swing_old.util.Analyzer import Analyzer
from datetime import datetime

"""
This pipeline scans a range of window sizes for a given inference method and generates roller objects for post analysis.

This is the 10 node variation of it.
"""
INPUT_PATH = "/home/jjw036/Swing/data/dream4/insilico_size10_"
INF_METHOD = "Dionesus"
OUTPUT_PATH = "/projects/p20519/roller_output/optimizing_window_size/" + INF_METHOD + "/insilico_size10_"
UNIQUE_NAME  = INF_METHOD + "insilico_size10_"
N_BOOT = 200
N_PERM = 5
RANDOM_WINDOWS = False

def main(window_size, n_trials, my_iterating_param):
    all_edge_lists = []
    overall_df = pd.DataFrame()
    for trials in range(0,n_trials):
        for network_index in range(1,6):
            file_path = INPUT_PATH + str(network_index) + "_timeseries.tsv"
            gene_start_column = 1
            time_label = "Time"
            separator = "\t"
            gene_end = None
            gold_standard = INPUT_PATH + str(network_index) + "_goldstandard.tsv"

            roller = Swing_old.Swing(file_path, gene_start_column, gene_end, time_label,separator,window_type=INF_METHOD)
            print("Overall Width: " + str(roller.overall_width))
            roller.zscore_all_data()

            roller.set_window(width=window_size)
            roller.create_windows(random_time = RANDOM_WINDOWS)
            roller.optimize_params()
            if window_size == roller.overall_width:
                crag = False
            else:
                crag = True
            roller.fit_windows(crag=crag)
            roller.rank_edges(permutation_n = N_PERM)

            analyzer = Analyzer(roller)
            result_df = analyzer.get_result_df()
            aggregated_best_windows = analyzer.aggregate_best_windows(result_df)
            
            print(aggregated_best_windows)
            #all_edge_lists.append(roller.window_list[0].results_table)
            all_edge_lists.append(analyzer.aggregated_edge_list)

            current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

            default_params = {'data_folder':OUTPUT_PATH, 'file_path':file_path,  'td_window':window_size,'min_lag':0,'max_lag':0,'n_trees':500,'permutation_n':N_PERM, 'lag_method':'no_lag','n_trials':n_trials, 'run_time':current_time,'iterating_param':my_iterating_param}

            default_params['auroc'] = result_df['auroc'][0]
            default_params['aupr'] = result_df['aupr'][0]
            run_result = pd.Series(default_params)
            overall_df = overall_df.append(run_result, ignore_index=True)
            """
            unique_filename = OUTPUT_PATH + str(network_index) + "/" + str(uuid.uuid4())
            with open(unique_filename, 'wb') as output:
                pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL)
            """
    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    overall_df.to_csv(OUTPUT_PATH+current_time+'.tsv',index=False,sep='\t')
    
    return(all_edge_lists)

    
if __name__ == "__main__":
    window_size = int(sys.argv[1])
    plt.ioff()
    n_trials = 2
    my_iterating_param = sys.argv[2]
    main(window_size, n_trials, my_iterating_param)

