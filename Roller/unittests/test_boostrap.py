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
    window_size = roll_me.overall_width
    alpha = 0.003
    roll_me.set_window(window_size)
    current_window = roll_me.get_window()
    coef = roll_me.fit(window_size, alpha=alpha, resamples=100)
    print coef.shape