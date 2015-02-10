__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
import Roller
import pandas as pd
import numpy as np
import numpy.testing as npt


class TestRFRWindow(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
        index = 0
        window_info = {"time_label": self.roller.time_label, "gene_start": self.roller.gene_start,
                       "gene_end": self.roller.gene_end, 'nth_window': index}
        self.test_rfr = Roller.RandomForestRegressionWindow(self.roller.get_window_raw(index),window_info)

    def test_get_coeffs(self):
        # With alpha at 0 everything should be nonzero except the diagonal values
        expected_non_zero = len(self.test_rfr.genes)**2-len(self.test_rfr.genes)
        calc_coeffs = self.test_rfr.get_coeffs(10)
        calc_non_zero = np.count_nonzero(calc_coeffs)
        self.assertTrue(expected_non_zero == calc_non_zero)

    def test_run_permutation_test(self):
        # The model must first be initialized
        self.test_rfr.initialize_params(10)
        self.test_rfr.fit_window()
        self.test_rfr.run_permutation_test(10)
        n_genes = len(self.test_rfr.genes)
        self.assertTrue(self.test_rfr.permutation_means.shape == (n_genes, n_genes))
        self.assertTrue(self.test_rfr.permutation_sd.shape == (n_genes, n_genes))

    def test_make_edge_table(self):
        self.test_rfr.initialize_params(n_trees=10)
        self.test_rfr.fit_window()
        self.test_rfr.run_permutation_test(10)
        n_boots = 13
        n_alphas = 20
        self.test_rfr.run_bootstrap(n_boots, n_alphas)
        self.test_rfr.make_edge_table()
        original_edge_order = self.test_rfr.edge_list
        new_edge_order = self.test_rfr.rank_edges()
        #todo: need a way to assert that these lists are not equal

if __name__ == '__main__':
    unittest.main()