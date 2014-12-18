import unittest
import Roller
import numpy as np
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
        window_values = np.random.random([10, 5])
        resample_window_values = self.roller.resample_window(window_values)

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
        # Generate test data matrix
        window_values = np.random.random([10, 5])
        max_random = 0.3

        # Get noisy values
        noise_values = self.roller.add_noise_to_window(window_values, max_random=max_random)

        # Make sure noise is within set range
        noise_magnitude = np.abs((noise_values-window_values)/window_values)
        self.assertTrue(np.all(noise_magnitude <= max_random))

if __name__ == '__main__':
    unittest.main()
