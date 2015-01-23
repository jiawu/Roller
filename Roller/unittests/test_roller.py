import unittest
import Roller
import numpy as np
import pandas as pd
import numpy.testing as npt
import matplotlib.pyplot as plt

class TestRoller(unittest.TestCase):
    def setUp(self):
        # Setup a different roller
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.dream_roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)

        # Only make 2 windows, so that that testing doesn't take forever
        self.dream_roller.set_window(self.dream_roller.overall_width-1)

    def test_get_window(self):
        print self.dream_roller.window_width
        print self.dream_roller.overall_width
        print self.dream_roller.get_n_windows()
        print self.dream_roller.get_window_raw(self.dream_roller.current_step)
        print self.dream_roller.get_window(self.dream_roller.current_step).head()


    def test_get_only_genes(self):
        only_genes = self.dream_roller.get_window(self.dream_roller.current_step)
        header = only_genes.columns.values
        correct_header = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10']
        self.assertTrue(np.array_equal(correct_header, header))

    def test_create_windows(self):
        self.dream_roller.create_windows_no_next()
        correct_n_windows = self.dream_roller.get_n_windows()
        n_windows = len(self.dream_roller.window_list)
        self.assertTrue(correct_n_windows == n_windows)

    def test_initialize_windows(self):
        self.dream_roller.create_windows_no_next()
        self.dream_roller.initialize_windows()
        for window in self.dream_roller.window_list:
            self.assertTrue(window.alpha is not None)
            self.assertTrue(window.beta_coefficients is not None)

    def test_rank_windows(self):
        self.dream_roller.create_windows_no_next()
        self.dream_roller.initialize_windows()
        self.dream_roller.rank_windows(10, 10)
        for window in self.dream_roller.window_list:
            self.assertTrue(window.edge_table is not None)

if __name__ == '__main__':
    unittest.main()
