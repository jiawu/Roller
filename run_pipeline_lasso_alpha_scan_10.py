import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Roller
import uuid
import pickle
"""
This pipeline scans a range of alphas for a given window size (LASSO) and identifies the optimal alpha.

This is the 10 node variation of it.

"""
INPUT_PATH = "data/dream4/insilico_size10_"
INF_METHOD = "lasso"
OUTPUT_PATH = "/projects/p20519/roller_output/optimizing_window_size/" + INF_METHOD + "/size10_"
UNIQUE_NAME = INF_METHOD + "size10_"
N_BOOT = 200
N_PERM = 200

if __name__ == "__main__":
    window_size = int(sys.argv[1])
    plt.ioff()

    for network_index in range(1,6):
      file_path = INPUT_PATH + str(network_index) + "_timeseries.tsv"

      gene_start_column = 1
      time_label = "Time"
      separator = "\t"
      gene_end = None
      gold_standard = INPUT_PATH + str(network_index) + "_goldstandard.tsv"

      roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
      print("Overall Width: " + str(roller.overall_width))
      roller.zscore_all_data()

      #### My goal here is to test the whole range of alphas for the full window ####

      roller.set_window(width=window_size)
      roller.create_windows(random_time=False)
      roller.optimize_params()
      for alpha in roller.window_list[0].cv_table['alpha']:
          print("current alpha: " + str(alpha))
          roller.fit_windows(alpha=alpha)
          roller.rank_edges(n_bootstraps=N_BOOT, permutation_n = N_PERM)
          roller.average_rank(rank_by='stability', ascending = False)
          #score some edge lists
          #first score the sorted average edge list
          averaged_score_dict = roller.score(roller.averaged_ranks, gold_standard)
          #next score each individual edge list

          unique_filename = OUTPUT_PATH+ str(network_index) + "/" + str(uuid.uuid4())
          with open(unique_filename, 'wb') as output:
            pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL)
