import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Swing_old
import uuid
import pickle
import pdb
"""
This pipeline scans a range of window sizes for a given inference method and generates roller objects for post analysis.

This is a pipeline scanning the Janes data.
"""
INPUT_PATH = "/projects/p20519/Swing/data/invitro/whitfield_muk"
INF_METHOD = "RandomForest"
OUTPUT_PATH = "/projects/p20519/roller_output/optimizing_window_size/" + INF_METHOD + "/whitfield_muk"
UNIQUE_NAME  = INF_METHOD + "whitfield_muk"
N_BOOT = 5
N_PERM = 5
RANDOM_WINDOWS = False


if __name__ == "__main__":
    
    window_size = int(sys.argv[1])
    print("Scanning with a window size of " + str(window_size))
    plt.ioff()

    file_path = INPUT_PATH + "_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    gold_standard = INPUT_PATH + "_goldstandard.tsv"

    roller = Swing_old.Swing(file_path, gene_start_column, gene_end, time_label,separator,window_type=INF_METHOD)
    print("Overall Width: " + str(roller.overall_width))
    roller.zscore_all_data()

    roller.set_window(width=window_size)
    roller.create_windows(random_time = RANDOM_WINDOWS)
    roller.optimize_params()
    if window_size == roller.overall_width:
        roller.fit_windows(crag=False)
        roller.rank_edges(permutation_n = N_PERM, crag=False)
    else:
        roller.fit_windows()
        roller.rank_edges(permutation_n = N_PERM)

    unique_filename = OUTPUT_PATH  + "/" + "_window_size_"+ str(window_size) +"_"+ str(uuid.uuid4())
    with open(unique_filename, 'wb') as output:
    #with gzip.GzipFile(unique_filename,'w') as output:
      pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL)
    #joblib.dump(roller, unique_filename)

