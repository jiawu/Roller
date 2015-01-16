__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
import numpy as np
import Roller
import pandas as pd
import numpy.testing as npt

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

if __name__ == '__main__':
    unittest.main()