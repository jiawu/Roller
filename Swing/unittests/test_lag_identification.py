import unittest
import Swing.util.lag_identification as lag_id
from Swing.util.Evaluator import Evaluator
import numpy as np
from random import randint
import numpy.testing as npt
import random
import pdb
import itertools
import pandas as pd

"""this test generally checks to see if aupr calculations are robust and if trivial cases evaluate as expected"""

class TestLagIdentification(unittest.TestCase):

    def setUp(self):
        self.experiments = lag_id.get_experiment_list("../../data/invitro/cantone_switchon_interpolated_timeseries.tsv", timepoints=18, perturbs=5)

        gold_standard_file = "../../data/invitro/cantone_switchon_interpolated_goldstandard.tsv"

        self.evaluator = Evaluator(gold_standard_file, sep = '\t')
        self.genes = list(self.experiments[0].columns.values)

    def test_calculate_edge_lag(self):
        exp_xcor = lag_id.xcorr_experiments(self.experiments, 1)
        edges = itertools.product(self.genes,self.genes)
        signed_edge_list = []
        for edge in edges:
            signed_edge_list.append(edge)
        signed_edge_df = pd.DataFrame({'regulator-target':signed_edge_list})
        signs = ['+' if bool(x<13) else '-' for x in range(25)]
        signed_edge_df['signs'] = signs


        filtered_lags=lag_id.calc_edge_lag(exp_xcor, self.genes, 0.1, 0.8, timestep=1, signed_edge_list = signed_edge_df, flat=False)
        print(filtered_lags)

if __name__ == '__main__':
    unittest.main()
