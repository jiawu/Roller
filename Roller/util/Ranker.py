__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

"""
NOTE TO JIA:

I never made this into a module, so some of the variables aren't passed to the function.
Therefore, not all of the functions are immediately usable. Still, you should get the idea.


"""

import numpy as np
from scipy import stats
import sys
from sklearn.utils import shuffle
import time
import dionesus as dio

def save_permutes(x, num_permutes, path):
    for ii in range(num_permutes):
        x_permuted = shuffle(x)
        np.save(path+"/permuted_%i"%ii, x_permuted)

def save_ptest_output(savepath, p, o_values, avg, stdev, zscores):
    permutations = p
    np.save(savepath+'/betas_%i' % permutations, o_values)
    np.save(savepath+'/avg_%i' % permutations, avg)
    np.save(savepath+'/stdev_%i' % permutations, stdev)
    np.save(savepath+'/zscores_%i' % permutations, zscores)

def low_mem_ptest(p, o_values, permutepath, savepath):
    permutations = p
    averages = np.zeros(o_values.shape)
    stdevs = np.zeros(o_values.shape)
    zscores = np.zeros(o_values.shape)
    for ii in range(num_variables):
        tic = time.clock()
        print "Doing %i permutations for gene %i" %(permutations,ii)
        cur_y = yMatrix[:, ii]
        cur_betas = np.array(np.matrix(o_values[:, ii]).T)
        for pp in range(permutations):
            cur_path = permutepath+"/permuted_%i.npy" %pp
            cur_x = np.load(cur_path)
            d.run_dionesus(cur_x, cur_y, x_match_y=False, show_output=False)
            perm_betas = np.array(d.beta_matrix.copy())
            cur_betas = np.hstack((cur_betas, perm_betas))
        averages[:, ii] = np.mean(cur_betas[:, 1:], 1)
        stdevs[:, ii] = np.std(cur_betas[:, 1:], 1, ddof=1)
        print time.clock()-tic


def boostrap():
    # Method specific
    # Generate n models
    # Keep beta value ranges for each model
    # Count how many times a beta value appears