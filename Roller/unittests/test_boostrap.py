__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
import numpy as np
import sys
from Roller.util.linear_wrapper import LassoWrapper
from Roller.util import Ranker

if __name__ == '__main__':
    file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    # Initial Model
    roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    window_size = 5
    alpha = 0.003
    roll_me.set_window(window_size)
    current_window = roll_me.get_window()
    a = roll_me.resample_window(current_window.values)
    b = roll_me.add_noise_to_window(a)
    print a[0]
    print b[0]


    sys.exit()
    o_coefs = roll_me.fit(window_size, alpha=alpha, resample=True)
    np.set_printoptions(suppress=True)
    print "Alpha: ", alpha
    print "Number of nonzero Edges:", np.count_nonzero(o_coefs)