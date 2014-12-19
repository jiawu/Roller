__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
from Roller.util.linear_wrapper import LassoWrapper
import numpy.testing as npt
import numpy as np
import pandas as pd

class TestLassoWrapper(unittest.TestCase):
    def setUp(self):
        # Load data
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        df = pd.DataFrame.from_csv(file_path, sep='\t')
        times = df.index.values[~np.isnan(df.index.values)]
        times_set = set(times)
        genes = df.columns.values
        replicates = len(times)/float(len(times_set))
        data = df.values

        # Remove NaNs from TSV
        data = data[~np.isnan(data).all(axis=1)]
        self.lassowrapper = LassoWrapper(data)

    def test_get_coef(self):
        n_samples, n_features = 5, 20
        X = np.random.randn(n_samples, n_features)
        coef = 3 * np.random.randn(n_features)
        inds = np.arange(n_features)
        y = np.dot(X, coef)

    def test_get_max_alpha(self):
        alpha_precision = 1e-9
        max_alpha = self.lassowrapper.get_max_alpha()
        num_coef_at_max_alpha = np.count_nonzero(self.lassowrapper.get_coeffs(max_alpha))
        num_coef_less_max_alpha = np.count_nonzero(self.lassowrapper.get_coeffs(max_alpha-alpha_precision))
        self.assertTrue(num_coef_at_max_alpha == 0)
        self.assertTrue(num_coef_less_max_alpha > 0)


if __name__ == '__main__':
    unittest.main()

"""
Old test code

if __name__ == '__main__':
    # Load data
    file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
    df = pd.DataFrame.from_csv(file_path, sep='\t')
    times = df.index.values[~np.isnan(df.index.values)]
    times_set = set(times)
    genes = df.columns.values
    replicates = len(times)/float(len(times_set))
    data = df.values

    # Remove NaNs from TSV
    data = data[~np.isnan(data).all(axis=1)]

    # Initialize lassowrapper
    lasso_wrapper = LassoWrapper(data)
    alpha = 0.0
    coef = lasso_wrapper.get_coeffs(alpha)
    m = lasso_wrapper.get_max_alpha()
    print m

"""