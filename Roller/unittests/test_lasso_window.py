__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
import numpy as np
import Roller
import pandas as pd
import numpy as np
from random import randint
import numpy.testing as npt
import random

class TestLassoWindow(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
        self.test_lassoWindow = Roller.LassoWindow(self.roller.current_window)

    def test_initialize_params_default(self):
        """ Test parameter initialization with default arguments """
        expected_alpha = self.test_lassoWindow.cv_select_alpha()
        self.test_lassoWindow.initialize_params()
        self.assertTrue(expected_alpha, self.test_lassoWindow.alpha)

    def test_initialize_params_manual_alpha(self):
        """ Test parameter initialization if passed an alpha value """
        expected_alpha = 5
        self.test_lassoWindow.initialize_params(expected_alpha)
        self.assertTrue(expected_alpha, self.test_lassoWindow.alpha)

    def test_sum_of_squares(self):
        data = np.reshape(np.arange(6), (3, 2))
        expected_ss = np.array([8, 8])
        calc_ss = self.test_lassoWindow.sum_of_squares(data)
        npt.assert_array_equal(calc_ss, expected_ss)

    def test_get_null_alpha(self):
        alpha_precision = 1e-9
        max_alpha = self.test_lassoWindow.get_null_alpha()
        num_coef_at_max_alpha = np.count_nonzero(self.test_lassoWindow.get_coeffs(max_alpha))
        num_coef_less_max_alpha = np.count_nonzero(self.test_lassoWindow.get_coeffs(max_alpha-alpha_precision))
        self.assertTrue(num_coef_at_max_alpha == 0)
        self.assertTrue(num_coef_less_max_alpha > 0)

    def test_cross_validate_alpha(self):
        test_alpha = 0.9 * self.test_lassoWindow.get_null_alpha()

        q_squared, model_q = self.test_lassoWindow.cross_validate_alpha(test_alpha, 3)

        # At least make sure there is a q_squared value for each gene
        self.assertTrue(len(q_squared) == len(self.test_lassoWindow.genes))

    def test_cv_select_alpha(self):
        calc_alpha, calc_cv_table = self.test_lassoWindow.cv_select_alpha()

        # Make sure alpha and cv_table are correct types
        self.assertTrue(type(calc_alpha) is np.float64 and type(calc_cv_table) is pd.DataFrame)

    def test_get_coeffs(self):
        # With alpha at 0 everything should be nonzero except the diagonal values
        expected_non_zero = len(self.test_lassoWindow.genes)**2-len(self.test_lassoWindow.genes)
        calc_coeffs = self.test_lassoWindow.get_coeffs(0)
        calc_non_zero = np.count_nonzero(calc_coeffs)
        self.assertTrue(expected_non_zero == calc_non_zero)

    def test_bootstrap_alpha(self):
        # The model must first be initialized
        self.test_lassoWindow.initialize_params()
        self.test_lassoWindow.fit_window()
        n_genes = len(self.test_lassoWindow.genes)
        n_boots = 100
        test_boot = self.test_lassoWindow.bootstrap_alpha(0.02, n_boots, 0.2)
        self.assertTrue(test_boot.shape == (n_genes, n_genes, n_boots))

    def test_run_bootstrap(self):
        # The model must first be initialized
        self.test_lassoWindow.initialize_params()
        self.test_lassoWindow.fit_window()
        n_boots = 13
        n_alphas = 20
        n_genes = len(self.test_lassoWindow.genes)
        self.test_lassoWindow.run_bootstrap(n_boots, n_alphas)
        self.assertTrue(self.test_lassoWindow.bootstrap_matrix.shape == (n_genes, n_genes, n_boots, n_alphas))

    def test_mean_calculation(self):
        # initialize result dict, which is a dict with mean, n, ss

        result = {'n': 0, 'mean': 0, 'ss': 0}
        new_samples = [randint(0, 100) for r in xrange(5)]
        new_result = self.test_lassoWindow.update_variance_1D(result, new_samples)

        correct_mean = np.mean(new_samples)
        self.assertEqual(new_result["mean"], correct_mean)

    def test_large_mean_calculation(self):
        # initialize result dict, which is a dict with mean, n, ss

        result = {'n': 0, 'mean': 0, 'ss': 0}
        new_samples = [randint(0, 100) for r in xrange(10000)]
        new_result = self.test_lassoWindow.update_variance_1D(result, new_samples)

        correct_mean = np.mean(new_samples)
        #get within 6 decimals?
        self.assertEqual("%.6f" % new_result["mean"], "%.6f" % correct_mean)

    def test_iterative_mean_calculation(self):
        new_result = {'n': 0, 'mean': 0, 'ss': 0}
        gold_samples = []

        for i in xrange(1000):
            new_samples = [random.uniform(0.001, 0.009) for r in xrange(10)]
            gold_samples = gold_samples + new_samples
            new_result = self.test_lassoWindow.update_variance_1D(new_result, new_samples)
        correct_mean = np.mean(gold_samples)
        print 'Correct Mean: ', (correct_mean)
        print 'Calculated Mean: ', (new_result['mean'])
        self.assertEqual("%.12f" % new_result["mean"], "%.12f" % correct_mean)

    def test_iterative_variance_calculation(self):
        new_result = {'n': 0, 'mean': 0, 'ss': 0}
        gold_samples = []

        for i in xrange(1000):
            new_samples = [random.uniform(0.001, 0.009) for r in xrange(10)]
            gold_samples = gold_samples + new_samples
            new_result = self.test_lassoWindow.update_variance_1D(new_result, new_samples)
        correct_variance = np.var(gold_samples, ddof=1)
        print 'Correct Var: ', (correct_variance)
        print 'Calculated Var: ', (new_result['variance'])
        self.assertEqual("%.12f" % new_result["variance"], "%.12f" % correct_variance)


    def test_variance_calculation(self):
        # initialize result dict, which is a dict with mean, n, ss

        result = {'n': 0, 'mean': 0, 'ss': 0}
        new_samples = [randint(0, 100) for r in xrange(100)]
        new_result = self.test_lassoWindow.update_variance_1D(result, new_samples)

        correct_var = np.var(new_samples, ddof=1)
        self.assertEqual("%.6f" % new_result["variance"], "%.6f" % correct_var)

    def test_incremental_variance_calculation(self):
        # initialize result dict, which is a dict with mean, n, ss

        result = {'n': 0, 'mean': 0, 'ss': 0}
        new_samples_A = [randint(0, 100) for r in xrange(2)]
        new_result = self.test_lassoWindow.update_variance_1D(result, new_samples_A)

        new_samples_B = [randint(0, 100) for r in xrange(2)]
        new_result = self.test_lassoWindow.update_variance_1D(new_result, new_samples_B)
        combined_samples = new_samples_A + new_samples_B

        correct_var = np.var(combined_samples, ddof=1)
        self.assertEqual("%.6f" % new_result["variance"], "%.6f" % correct_var)

    def test_incremental_mean_calculation(self):
        # initialize result dict, which is a dict with mean, n, ss

        result = {'n': 0, 'mean': 0, 'ss': 0}
        new_samples_A = [randint(0, 10000) for r in xrange(2)]
        new_result = self.test_lassoWindow.update_variance_1D(result, new_samples_A)

        new_samples_B = [randint(0, 10000) for r in xrange(2)]
        new_result = self.test_lassoWindow.update_variance_1D(new_result, new_samples_B)
        combined_samples = new_samples_A + new_samples_B

        correct_var = np.mean(combined_samples)
        self.assertEqual("%.6f" % new_result["mean"], "%.6f" % correct_var)

    def test_incremental_2D_array_mean(self):
        zeros = np.zeros((10, 10))
        result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        new_samples_A = np.random.random((10, 10))
        A_list = []
        A_list.append(new_samples_A)
        new_result = self.test_lassoWindow.update_variance_2D(result, A_list)

        new_samples_B = np.random.random((10, 10))
        B_list = []
        B_list.append(new_samples_B)
        new_result = self.test_lassoWindow.update_variance_2D(new_result, B_list)
        combined_samples = np.dstack((new_samples_A, new_samples_B))
        # dstack coordinates are formatted as follows:
        # (y coord, x coord, z coord)
        # or alternatively, (row index, col index, sample index)

        correct_var = np.mean(combined_samples, axis=2)

        npt.assert_array_almost_equal(new_result["mean"], correct_var, decimal=12)

    def test_iterative_2D_array_mean(self):
        zeros = np.zeros((10, 10))
        new_result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        A_list = []

        combined_samples = np.empty([10, 10])
        for i in xrange(1000):
            new_samples_A = np.random.random((10, 10))
            A_list.append(new_samples_A)
            temp_list = []
            temp_list.append(new_samples_A)
            new_result = self.test_lassoWindow.update_variance_2D(new_result, temp_list)

        combined_samples = np.dstack(A_list)
        # dstack coordinates are formatted as follows:
        # (y coord, x coord, z coord)
        # or alternatively, (row index, col index, sample index)

        correct_mean = np.mean(combined_samples, axis=2)

        npt.assert_array_almost_equal(new_result["mean"], correct_mean, decimal=12)

    def test_iterative_2D_array_variance(self):
        zeros = np.zeros((10, 10))
        new_result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        A_list = []

        combined_samples = np.empty([10, 10])
        for i in xrange(1000):
            new_samples_A = np.random.random((10, 10))
            A_list.append(new_samples_A)
            temp_list = []
            temp_list.append(new_samples_A)
            new_result = self.test_lassoWindow.update_variance_2D(new_result, temp_list)

        combined_samples = np.dstack(A_list)
        # dstack coordinates are formatted as follows:
        # (y coord, x coord, z coord)
        # or alternatively, (row index, col index, sample index)

        correct_var = np.var(combined_samples, axis=2, ddof=1)

        npt.assert_array_almost_equal(new_result["variance"], correct_var, decimal=12)

    def test_incremental_2D_array_variance(self):
        zeros = np.zeros((10, 10))
        result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        new_samples_A = np.random.random((10, 10))
        A_list = []
        A_list.append(new_samples_A)
        new_result = self.test_lassoWindow.update_variance_2D(result, A_list)

        new_samples_B = np.random.random((10, 10))
        B_list = []
        B_list.append(new_samples_B)

        new_result = self.test_lassoWindow.update_variance_2D(new_result, B_list)
        combined_samples = np.dstack((new_samples_A, new_samples_B))
        # dstack coordinates are formatted as follows:
        # (y coord, x coord, z coord)
        # or alternatively, (row index, col index, sample index)

        correct_var = np.var(combined_samples, axis=2, ddof=1)

        npt.assert_array_almost_equal(new_result["variance"], correct_var, decimal=12)

if __name__ == '__main__':
    unittest.main()