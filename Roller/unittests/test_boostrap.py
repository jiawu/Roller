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
    roll_me.set_window(roll_me.overall_width)
    window = roll_me.get_window()
    lasso = LassoWrapper(window.values)
    max_alpha = lasso.get_max_alpha()
    np.set_printoptions(suppress=True)
    coefs = lasso.get_coeffs(max_alpha-1e-9)

    boots = 10000
    max_random = 1e-7
    o = np.ones((5,5))
    o = np.multiply(o, np.random.uniform(low=1-max_random, high=1+max_random, size=o.shape))
    print o
    exists = np.empty((boots))
    for ii in range(boots):
        n, p = window.values.shape
        resample_values = window.values.copy()
        # For each column randomly choose samples
        #resample_values = np.array([np.random.choice(window.values[:,ii], size=n) for ii in range(p)]).T
        noise = np.random.uniform(low=1-max_random, high=1+max_random, size=window.values.shape)
        noisy_values = np.multiply(resample_values, noise)
        lasso = LassoWrapper(noisy_values)
        coefs = lasso.get_coeffs(max_alpha)
        exists[ii] = coefs[0,4]
    print np.count_nonzero(exists)/float(boots)





    sys.exit()
    alphas, freq_matrix = boot.run_bootstrap(roll_me.overall_width, 1, 100, noise=0)
    edge_nonzeros = np.sum(boot.bootstrap_matrix[4,0,0] !=0, axis=0)
    print edge_nonzeros
    sums = np.sum(freq_matrix, axis=3)


    sys.exit()
    # Graph
    Grapher.plot_stability(freq_matrix, alphas, 0)
    sys.exit()
    for ii in range(17):
        Grapher.plot_stability(freq_matrix, alphas, ii)