__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
import numpy as np
from Roller import Window
import pandas as pd
import numpy.testing as npt

class TestWindow(unittest.TestCase):
    def setUp(self):
        pass

    def test_resample_window(self):
        # Generate test data matrix
        window = Window(pd.DataFrame(np.random.random([5, 3]), columns=['a', 'b', 'c'], index=['A', 'B', 'C', 'D', 'E']))
        window.resample_window()

        # Confirm shapes are true
        self.assertTrue(window.window_values.shape == window.resample_values.shape)
        num_rows, num_columns = window.window_values.shape

        # Verify that the resampled matrix doesn't equal the original matrix
        self.assertFalse(np.array_equal(window.resample_values, window.window_values))

        # Verify that values in each column of resampled matrix are values in the same column of the original window
        truth_table = np.array(
            [[value in window.window_values[:, column] for value in window.resample_values[:, column]] for column in
             range(num_columns)]).T
        self.assertTrue(np.all(truth_table))

    def test_add_noise_to_window(self):
        # Generate test data frame
        window = Window(pd.DataFrame(np.random.random([10, 5])))
        max_random = 0.3

        # Get noisy values
        noise_values = self.roller.add_noise_to_window(window, max_random=max_random)

        # Make sure noise is within set range
        noise_magnitude = np.abs((noise_values-window)/window)
        self.assertTrue(np.all(noise_magnitude <= max_random))


if __name__ == '__main__':
    unittest.main()
