import Roller
from pylab import *

if __name__ == "__main__":
    file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    gold_standard = "data/dream4/insilico_size10_1_goldstandard.tsv"
    window_size = 21
    image_file_path = "insilico_size10_1_alphas"
    roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    print("Overall Width: " + str(roller.overall_width))
    roller.zscore_all_data()

    #### My goal here is to test the whole range of alphas for the full window ####
    alpha_list = []
    aupr_list = []

    roller.create_windows(width=window_size)
    roller.optimize_params()
    for alpha in roller.window_list[0].cv_table['alpha']:
        print("current alpha: " + str(alpha))
        roller.fit_windows(alpha=alpha)
        roller.rank_edges(n_bootstraps=200, permutation_n = 200)
        roller.average_rank(rank_by='stability', ascending = True)
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
    #plot aupr and alpha
    plot(alpha_list, aupr_list)
    title_string = file_path + " window size: " + str(window_size)
    title(title_string)
    xlabel('alpha')
    ylabel('aupr')
    image_save = image_file_path + "w" + str(window_size) + ".png"
    savefig(image_save)
