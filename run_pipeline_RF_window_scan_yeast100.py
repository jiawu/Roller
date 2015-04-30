import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Roller
import uuid
import pickle
"""
This pipeline scans a range of window sizes for a given inference method and generates roller objects for post analysis.

This is the 100 node yeast variation of it.
"""
INPUT_PATH = "data/dream4/yeast_size100_"
INF_METHOD = "RandomForest"
OUTPUT_PATH = "/projects/p20519/roller_output/optimizing_window_size/" + INF_METHOD + "/yeast_size100_"
UNIQUE_NAME  = INF_METHOD + "yeast_size100_"
N_BOOT = 200
N_PERM = 200
RANDOM_WINDOWS = False


if __name__ == "__main__":
    window_size = int(sys.argv[1])
    plt.ioff()

    for network_index in range(1,2):
      file_path = INPUT_PATH + str(network_index) + "_timeseries.tsv"
      gene_start_column = 1
      time_label = "Time"
      separator = "\t"
      gene_end = None
      gold_standard = INPUT_PATH + str(network_index) + "_goldstandard.tsv"

      roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label,separator,window_type=INF_METHOD)
      print("Overall Width: " + str(roller.overall_width))
      roller.zscore_all_data()

      roller.set_window(width=window_size)
      roller.create_windows(random_time = RANDOM_WINDOWS)
      roller.optimize_params()
      roller.fit_windows()
      roller.rank_edges(permutation_n = N_PERM)

      unique_filename = OUTPUT_PATH + str(network_index) + "/" + str(uuid.uuid4())
      with open(unique_filename, 'wb') as output:
        pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL)

