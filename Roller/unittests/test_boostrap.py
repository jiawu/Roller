__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
import numpy as np
import sys
from Roller.util.linear_wrapper import LassoWrapper
from Roller.util import Ranker
from Roller.util import Grapher
from scipy import stats
import matplotlib.pyplot as plt

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
    norm_values = stats.zscore(window.values.copy(), axis=1, ddof=1)
    lasso = LassoWrapper(norm_values)
    max_alpha = lasso.get_max_alpha()

    print max_alpha

    boots = 100
    max_random = 0.1
    alphas = 10
    #alpha_range = [max_alpha/10**aa for aa in range(alphas)]
    alpha_range = np.linspace(0, max_alpha, 100)
    print alpha_range
    freqs = np.empty((len(alpha_range),4))
    for kk, alpha in enumerate(alpha_range):
        print kk
        exists = np.empty((boots, 4))
        for ii in range(boots):
            n, p = norm_values.shape
            # For each column randomly choose samples
            resample_values = np.array([np.random.choice(norm_values[:, jj], size=n) for jj in range(p)]).T
            noise = np.random.uniform(low=1-max_random, high=1+max_random, size=resample_values.shape)
            noisy_values = np.multiply(norm_values, noise)
            noisy_resample_values = np.multiply(resample_values, noise)

            #Vanilla values
            v_coefs = LassoWrapper(norm_values).get_coeffs(alpha)

            # Noisy values
            n_coefs = LassoWrapper(noisy_values).get_coeffs(alpha)

            # Resample values
            r_coefs = LassoWrapper(resample_values).get_coeffs(alpha)

            # Noisy Resample
            n_r_coefs = LassoWrapper(noisy_resample_values).get_coeffs(alpha)

            exists[ii] = [v_coefs[0,4]!=0, n_coefs[0,4]!=0, r_coefs[0,4]!=0, n_r_coefs[0,4]!=0]
        freqs[kk] = np.sum(exists, axis=0)/float(boots)

    print freqs
    plt.plot(alpha_range, freqs, 'o-')
    plt.legend(['Reg', 'N', 'R', 'N+R'])
    plt.show()




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