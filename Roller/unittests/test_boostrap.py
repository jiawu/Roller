__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
import numpy as np
import sys
from Roller.util.linear_wrapper import LassoWrapper
from Roller.util import Ranker
from Roller.util import Grapher

if __name__ == '__main__':
    file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    # Initialize Model
    roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    boot = Ranker.Bootstrapper(roll_me)
    alphas, freq_matrix = boot.run_bootstrap(5, 15, 100)
    sums = np.sum(freq_matrix, axis=3)

    #sys.exit()
    # Graph
    Grapher.plot_stability(freq_matrix, alphas, 8)