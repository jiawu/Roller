import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import Roller
import uuid
import pickle
if __name__ == "__main__":
    window_size = int(sys.argv[1])
    plt.ioff()
    file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    gold_standard = "data/dream4/insilico_size10_1_goldstandard.tsv"
    image_file_path = "insilico_size10_1_alphas"
    roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    print("Overall Width: " + str(roller.overall_width))
    roller.zscore_all_data()

    #### My goal here is to test the whole range of alphas for the full window ####
    alpha_list = []
    aupr_list = []

    roller.set_window(width=window_size)
    roller.create_windows()
    roller.optimize_params()
    for alpha in roller.window_list[0].cv_table['alpha']:
        print("current alpha: " + str(alpha))
        roller.fit_windows(alpha=alpha)
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
        alpha_list.append(alpha)
        aupr_list.append(max(score_dict['aupr']))
        unique_filename = "/projects/p20519/Roller_outputs/"+ str(uuid.uuid4())
        with open(unique_filename, 'wb') as output:
          pickle.dump(roller,output, pickle.HIGHEST_PROTOCOL) 
    #plot aupr and alpha
    fig = plt.figure()
    plt.plot(alpha_list, aupr_list)
    title_string = file_path + " window size: " + str(window_size)
    plt.title(title_string)
    plt.xlabel('alpha')
    plt.ylabel('aupr')
    image_save = image_file_path + "w" + str(window_size) + ".png"
    fig.savefig(image_save)

