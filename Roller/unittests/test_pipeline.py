import unittest
import numpy as np
import Roller

class TestWindow(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)

    def test_crap(self):
        self.roller.create_windows()
        self.roller.optimize_params()
        self.roller.fit_windows()
        self.roller.rank_edges()

if __name__ == '__main__':
    unittest.main()
