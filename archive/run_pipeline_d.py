import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Swing
import pdb
import uuid
import pickle
if __name__ == "__main__":
    file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    gold_standard = "data/dream4/insilico_size10_1_goldstandard.tsv"
    max_window_size = 21
    image_file_path = "insilico_size10_1_windows"
    roller = Swing.Swing(file_path, gene_start_column, gene_end, time_label, separator)
    print("Overall Width: " + str(roller.overall_width))
    roller.zscore_all_data()

    #### My goal here is to test the whole range of window sizes and compare aupr ####
    window_list = []
    aupr_list = []
    best_aupr_list = []
    for window_size in range(2,max_window_size+1):
      roller.set_window(width=window_size)
      roller.create_windows()

      roller.optimize_params()
      best_alpha = roller.window_list[0].alpha
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
      window_list.append(window_size)
      aupr_list.append(max(score_dict['aupr']))
      best_score = max([max(x['aupr']) for x in score_list])
      best_aupr_list.append(best_score)
      unique_filename = "/projects/p20519/Roller_outputs/"+ str(uuid.uuid4())
      with open(unique_filename, 'wb') as output:
        pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL) 

    #plot aupr and alpha
    fig1 = plt.figure(0)
    plt.plot(window_list, aupr_list)
    title_string = file_path + " window size: " + str(window_size)
    plt.title(title_string)
    plt.xlabel('window size')
    plt.ylabel('aupr')
    image_save = image_file_path + "w" + str(window_size) + ".png"
    plt.savefig(image_save)

    fig2 = plt.figure(1)
    plt.plot(window_list, best_aupr_list)
    title_string = file_path + " best window size: " + str(window_size)
    plt.title(title_string)
    plt.xlabel('window size')
    plt.ylabel('best aupr')
    image_save2 = image_file_path + "best_w" + str(window_size) + ".png"
    plt.savefig(image_save2)
    
