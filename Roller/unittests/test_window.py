__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import unittest
import numpy as np
import Roller
import pandas as pd
import numpy.testing as npt

class TestWindow(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
        self.test_window = Roller.Window(self.roller.current_window)

    def test_resample_window(self):

        resampled_values = self.test_window.resample_window()

        # Confirm shapes are true
        self.assertTrue(self.test_window.window_values.shape == resampled_values.shape)
        num_rows, num_columns = self.test_window.window_values.shape

        # Verify that the resampled matrix doesn't equal the original matrix
        self.assertFalse(np.array_equal(resampled_values, self.test_window.window_values))

        # Verify that values in each column of resampled matrix are values in the same column of the original window
        truth_table = np.array(
            [[value in self.test_window.window_values[:, column] for value in resampled_values[:, column]] for column in
             range(num_columns)]).T
        self.assertTrue(np.all(truth_table))

    def test_add_noise_to_window(self):
        # Generate test data frame
        original_values = self.test_window.window_values
        max_random = 0.3

        # Get noisy values
        noise_values = self.test_window.add_noise_to_values(original_values, max_random=max_random)

        # Make sure noise is within set range
        noise_magnitude = np.abs((noise_values-original_values)/original_values)
        self.assertTrue(np.all(noise_magnitude <= max_random))


class TestLasso_Window(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
        self.test_lasso_window = Roller.Lasso_Window(self.roller.current_window)

    def test_sum_of_squares(self):
        data = np.reshape(np.arange(6), (3, 2))
        expected_ss = np.array([8, 8])
        calc_ss = self.test_lasso_window.sum_of_squares(data)
        npt.assert_array_equal(calc_ss, expected_ss)

    def test_get_null_alpha(self):
        alpha_precision = 1e-9
        max_alpha = self.test_lasso_window.get_null_alpha()
        num_coef_at_max_alpha = np.count_nonzero(self.test_lasso_window.fit_window(max_alpha))
        num_coef_less_max_alpha = np.count_nonzero(self.test_lasso_window.fit_window(max_alpha-alpha_precision))
        self.assertTrue(num_coef_at_max_alpha == 0)
        self.assertTrue(num_coef_less_max_alpha > 0)

    def test_cross_validate_alpha(self):
        test_alpha = 0.9 * self.test_lasso_window.get_null_alpha()

        q_squared, model_q = self.test_lasso_window.cross_validate_alpha(test_alpha, 3)

        # At least make sure there is a q_squared value for each gene
        self.assertTrue(len(q_squared) == len(self.test_lasso_window.genes))

    def test_cv_select_alpha(self):
        alpha_range = np.linspace(0, self.test_lasso_window.null_alpha)
        folds = 3
        cv_results = np.array(
            [self.test_lasso_window.cross_validate_alpha(alpha, folds, True) for alpha in alpha_range])
        column_labels = np.append(self.test_lasso_window.genes+"_Q^2", "Model_Q^2")
        df = pd.DataFrame(cv_results, columns=column_labels)
        df.insert(0, 'alpha', alpha_range)
        print df.head()
        #self.test_lasso_window.cv_select_alpha(alpha_range)

if __name__ == '__main__':
    unittest.main()
