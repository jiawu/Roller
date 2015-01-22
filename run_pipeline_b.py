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
    roller.rank_edges(n_bootstraps=100)
    roller.average_rank(rank_by='p-value-perm')
    #score some edge lists
    #first score the sorted average edge list
    averaged_aupr = roller.score(roller.averaged_ranks, gold_standard)
    #next score each individual edge list
    aupr_list = []
    for window in roller.window_list:
        aupr = roller.score(window.results_table,gold_standard)
        aupr_list.append(aupr)
    pdb.set_trace()

