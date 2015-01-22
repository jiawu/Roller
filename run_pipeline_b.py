import Roller
import pdb

if __name__ == "__main__":
    file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    gold_standard = "data/dream4/insilico_size10_1_goldstandard.tsv"

    roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    print("Overall Width: " + str(roller.overall_width))
    roller.create_windows(width=21)
    roller.optimize_params()
    roller.fit_windows()
    roller.rank_edges(n_bootstraps=200)
    roller.average_rank(rank_by='stability')
    #score some edge lists
    #first score the sorted average edge list
    averaged_score_dict = roller.score(roller.averaged_ranks, gold_standard)
    #next score each individual edge list
    score_list = []
    for window in roller.window_list:
        score_dict = roller.score(window.results_table,gold_standard)
        score_list.append(score_dict)
    pdb.set_trace()

