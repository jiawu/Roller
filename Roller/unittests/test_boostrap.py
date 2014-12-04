__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
import numpy as np
import sys
from util.linear_wrapper import LassoWrapper
from util import Ranker

if __name__ == '__main__':
    file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    window_size = roll_me.overall_width
    alpha=0.003
    coefs = roll_me.fit(window_size, alpha=alpha)
    np.set_printoptions(suppress=True)
    print "Alpha: ", alpha
    print "Number of nonzero Edges:", np.count_nonzero(coefs)