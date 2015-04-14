import unittest
import numpy as np
import Roller
import random
from random import randint
import numpy.testing as npt
import pdb

class TestWindow(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator, window_type = "Lasso")
        self.roller.set_window(width=20)
        self.roller.create_windows()
        self.test_window = self.roller.window_list[0]

    def test_model_is_saved(self):
        model_list = self.test_window.model
        n_genes = self.test_window.n_genes
        self.assertTrue(len(model_list),n_genes)


    def test_prediction(self):

if __name__ == '__main__':
    unittest.main()

