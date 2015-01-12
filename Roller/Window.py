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

    def add_noise_to_window(self, window, max_random=0.2):
        """
        Add uniform noise to each value
        :param window: dataframe

        :param max_random: float
            Amount of noise to add to each value, plus or minus
        :return: array

        """
        noise = np.random.uniform(low=1-max_random, high=1+max_random, size=window.shape)
        noisy_values = np.multiply(window, noise)
        return noisy_values