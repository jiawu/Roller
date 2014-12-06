__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import unittest
from Roller.util.linear_wrapper import LassoWrapper
import numpy as np
import pandas as pd

""""
########################################################################################################################
########################################################################################################################

NOTE: This will be made into a complete unittest when Justin learns how to do them. For now its just going to
be a script used to make sure the LassoWrapper.get_maximum_alpha function works properly

########################################################################################################################
########################################################################################################################
class TestLassoWrapper(unittest.TestCase):
    def setUp(self):
        self.lassowrapper = LassoWrapper()



if __name__ == '__main__':
    unittest.main()

"""

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