__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

"""
NOTE TO JIA:

I never made this into a module, so some of the variables aren't passed to the function.
Therefore, not all of the functions are immediately usable. Still, you should get the idea.



Brainstorming:
    1. get a Roller
    2. For each window
    3. Permute the data
    4. calculate coefficients for the permuted data
    5. create a X by X matrix. within the X by X matrix, each unit is a dict with a mean and variance.
    6. append to the X by X matrix by doing the inline calculations
    7. repeat steps 3-6 N number of times.
    8. plot distribution of mean and standard deviation for each edge.
    9. compare cutoffs

"""

import numpy as np
from scipy import stats
import sys
from sklearn.utils import shuffle
import time
import dionesus as dio
from Roller.util.linear_wrapper import LassoWrapper

class Permuter():
    def __init__(self):
        self.permutation_means = None
        self.permutation_sd = None

    def set_data(self):
        pass

    def run_permutation_test(self, roller_object, alpha = 10, permutation_n=1000):
        total_window_number = roller_object.get_n_windows()
        n_genes = roller_object.get_n_genes()
        #initialize permutation results array
        self.permutation_means = np.empty((n_genes,n_genes, total_window_number))
        self.permutation_sd = np.empty((n_genes,n_genes,total_window_number))
        roller_object.reset()
        for nth_window in range(0, total_window_number):
            current_window = roller_object.get_window()
            permuted_window = current_window.copy()
            for nth_perm in range(0, permutation_n):
                self.permute_data(permuted_window.values)
                current_lasso = LassoWrapper(permuted_window.values)
                permuted_coeffs = current_lasso.get_coeffs(alpha)
            roller_object.next()

#        coeff_matrix_3d = np.empty((n_genes,n_genes,permutation_n, total_window_number))

    def calculate_inline(nth_window,):
        pass

    def permute_data(self, array):
        """Warning: Modifies data in place. also remember the """
        #new_array   = array
        _ = [np.random.shuffle(i) for i in array]
        return array

    def update_variance(self, prev_result, new_samples):
        """incremental calculation of means: accepts new_samples, which is a list of samples. then calculates a new mean. this is a useful function for calculating the means of large arrays"""
        n = float(prev_result["n"])
        mean = float(prev_result["mean"])
        sum_squares = float(prev_result["ss"])

        for x in new_samples:
            n = n + 1
            #delta = float(x) - mean
            old_mean = mean
            mean = old_mean + (float(x)-old_mean)/n
            sum_squares = sum_squares + (float(x)-mean)*(float(x)-old_mean)

        if (n < 2):
            return 0

        variance = sum_squares/(n-1)
        result = {  "mean": mean,
                    "ss": sum_squares,
                    "variance": variance,
                    "n": n}
        return result


    def get_true_edges(self, fp_cutoff, edges):
        pass

    def get_edge_list(self):
        pass

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
