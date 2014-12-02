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

    # Calculate zscores
    zscores = (o_values-averages)/stdevs

    save_ptest_output(savepath, p, o_values, averages, stdevs, zscores)

def ptest(p, o_values, savepath):
    permutations = p
    print 'Starting %i permutations' %permutations
    for p in range(permutations):
        print 'Permutation %i' %p
        x_permuted = shuffle(xMatrix)
        d.run_dionesus(x_permuted, yMatrix, show_output=False)
        o_values = np.dstack((o_values, np.array(d.beta_matrix.copy())))
    print o_values.shape

    # Calculate average
    o_values = o_values[:, :, 0]
    avg = np.mean(o_values[:, :, 1:], 2)

    # Calculate standard deviation
    stdev = np.std(o_values[:, :, 1:], 2, ddof=1)

    # Calculate zscore
    zscores = (o_values-avg)/stdev

    save_ptest_output(savepath, p, o_values, avg, stdev, zscores)