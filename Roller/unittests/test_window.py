__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
import numpy as np
from Roller import Window
import pandas as pd
import numpy.testing as npt

class TestWindow(unittest.TestCase):
    def setUp(self):
        df = pd.DataFrame(np.random.random([10,10]))
        self.window = Window(df)

    def test_resample_window(self):
        # Generate test data matrix
        window = pd.DataFrame(np.random.random([5, 3]), columns=['a', 'b', 'c'], index=['A', 'B', 'C', 'D', 'E'])
        resample_window = self.roller.resample_window(window)

        # Make sure columns and index values remain the same
        npt.assert_array_equal(window.columns.values, resample_window.columns.values)
        npt.assert_array_equal(window.index.values, resample_window.index.values)

        window_values = window.values
        resample_window_values = resample_window.values

        # Confirm shapes are true
        self.assertTrue(window_values.shape == resample_window_values.shape)
        num_rows, num_columns = window_values.shape

        # Verify that the resampled matrix doesn't equal the original matrix
        self.assertFalse(np.array_equal(resample_window_values, window_values))

        # Verify that values in each column of resampled matrix are values in the same column of the original window
        truth_table = np.array([[value in window_values[:, column] for value in resample_window_values[:, column]]
                                for column in range(num_columns)]).T
        self.assertTrue(np.all(truth_table))

    def test_add_noise_to_window(self):
        # Generate test data frame
        window = pd.DataFrame(np.random.random([10, 5]))
        max_random = 0.3

        # Get noisy values
        noise_values = self.roller.add_noise_to_window(window, max_random=max_random)

        # Make sure noise is within set range
        noise_magnitude = np.abs((noise_values-window)/window)
        self.assertTrue(np.all(noise_magnitude <= max_random))


if __name__ == '__main__':
    unittest.main()
