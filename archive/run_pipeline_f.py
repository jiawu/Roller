import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Swing_old
import uuid
import pickle
"""
This pipeline scans a range of
"""




if __name__ == "__main__":
    window_size = int(sys.argv[1])
    plt.ioff()

    for file_path_index in range(1,6):
      file_path = "data/dream4/insilico_size100_" + str(file_path_index) + "_timeseries.tsv"
      #file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
      gene_start_column = 1
      time_label = "Time"
      separator = "\t"
      gene_end = None
      gold_standard = "data/dream4/insilico_size100_" + str(file_path_index) + "_goldstandard.tsv"
      #gold_standard = "data/dream4/insilico_size10_1_goldstandard.tsv"
      image_file_path = "insilico_size100_" + str(file_path_index) + "_alphas"
      #image_file_path = "insilico_size10_1_alphas"
      roller = Swing_old.Swing(file_path, gene_start_column, gene_end, time_label, separator)
      print("Overall Width: " + str(roller.overall_width))
      roller.zscore_all_data()

      #### My goal here is to test the whole range of alphas for the full window ####
      alpha_list = []
      aupr_list = []

      roller.set_window(width=window_size)
      roller.create_windows()
      roller.optimize_params()
      roller.fit_windows()
      roller.rank_edges(n_bootstraps=500, permutation_n = 500)
      roller.average_rank(rank_by='stability', ascending = False)
      #score some edge lists
      #first score the sorted average edge list
      averaged_score_dict = roller.score(roller.averaged_ranks, gold_standard)
      #next score each individual edge list
      score_list = []
      for window in roller.window_list:
          score_dict = roller.score(window.results_table,gold_standard)
          score_list.append(score_dict)
      aupr_list.append(max(score_dict['aupr']))
      unique_filename = "/projects/p20519/Roller_outputs_L_yeast100/"+ str(uuid.uuid4())
      with open(unique_filename, 'wb') as output:
        pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL)
