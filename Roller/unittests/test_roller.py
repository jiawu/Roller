import unittest
import Roller
import numpy as np
import pandas as pd
import numpy.testing as npt

class TestRoller(unittest.TestCase):
    def setUp(self):
        self.roller = Roller.Roller('../../data/emt/compressed_katrina_data.txt', 5, None, "time", " ")

    def test_next(self):
        self.roller.next()
        window = self.roller.get_window_raw()
        time_slice = window['time'].unique()
        correct_window = [1,2,3]
        self.assertTrue(np.array_equal(correct_window, time_slice))

    def test_get_only_genes(self):
        only_genes = self.roller.get_window()
        header = only_genes.columns.values
        correct_header = ['AP1','AP2','AP3','AP4', 'AR', 'Bcat', 'Brachyury', 'cmyc', 'CRE', 'E2F','ELK1', 'ER', 'ETS1', 'FOXA', 'FOXO3A', 'GATA1', 'GATA2', 'GATA3','GLI', 'GR', 'HIF1', 'HNF1A', 'HOXA1', 'HSE', 'KLF4', 'LHX8','MEF2', 'MNX', 'MNX1', 'MYB', 'NANOG', 'NFAT', 'NFKb', 'NOBOX','Notch1', 'OCT', 'P53', 'PAX1', 'PEA3', 'PR', 'PTTG', 'RAR','RUNX1', 'RUNX2', 'SMAD1', 'SMAD3', 'SOX', 'SP1', 'SRF', 'STAT1','STAT3', 'STAT4', 'STAT5', 'VDR', 'WT1', 'YY1']
        self.assertTrue(np.array_equal(correct_header, header))

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
        self.assertTrue(window_values.shape==resample_window_values.shape)
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

    def test_cross_validate_window_alpha(self):
        window_values = pd.DataFrame(np.random.random([10, 5]))
        self.roller.cross_validate_window_alpha(window_values)

if __name__ == '__main__':
    unittest.main()
