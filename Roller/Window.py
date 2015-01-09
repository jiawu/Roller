__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import numpy as np
import pandas as pd


class Window(object):
    def __init__(self, dataframe, alpha=None):
        self.df = dataframe
        self.alpha = alpha
        self.window_values = dataframe.values

    def resample_window(self):
        """
        Resample window values, along a specific axis
        :param window_values: array

        :return: array
        """
        n, p = self.window_values.shape

        # For each column randomly choose samples
        resample_values = np.array([np.random.choice(self.window_values[:, ii], size=n) for ii in range(p)]).T

        resample_window = pd.DataFrame(resample_values, columns=self.df.columns.values.copy(),
                                       index=self.df.index.values.copy())

        return resample_window
