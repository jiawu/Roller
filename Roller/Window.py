__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import numpy as np
import pandas as pd
from util.linear_wrapper import LassoWrapper


class Window(object):
    def __init__(self, dataframe):
        self.df = dataframe
        self.window_values = dataframe.values
        self.samples = dataframe.index.values
        self.genes = dataframe.index.values
        self.edge_table = pd.DataFrame()

    def fit_window(self):
        """
        Fit the window with the specified window

        """
        pass

    def resample_window(self):
        """
        Resample window values, along a specific axis
        :param window_values: array

        :return: array
        """
        n, p = self.window_values.shape

        # For each column randomly choose samples
        resample_values = np.array([np.random.choice(self.window_values[:, ii], size=n) for ii in range(p)]).T

        #resample_window = pd.DataFrame(resample_values, columns=self.df.columns.values.copy(),
        #                               index=self.df.index.values.copy())
        return resample_values

    def add_noise_to_values(self, window_values, max_random=0.2):
        """
        Add uniform noise to each value
        :param window: dataframe

        :param max_random: float
            Amount of noise to add to each value, plus or minus
        :return: array

        """

        noise = np.random.uniform(low=1-max_random, high=1+max_random, size=window_values.shape)
        noisy_values = np.multiply(window_values, noise)
        return noisy_values

    def permutation_test(self):
        pass

    def bootstrap(self):
        pass

    def rank_edges(self, method):
        pass

    def initialize_params(self):
        pass

class Lasso_Window(Window):
    def __init__(self, dataframe, alpha=None):
        super(Lasso_Window, self).__init__(dataframe)
        self.alpha = alpha

    def fit_window(self, window_values, alpha):
        """
        Given a window get the lasso coefficients
        :param window_values: array-like
            The values to use for the fit
        :param alpha: float
            Value to use for lasso regression
        :return: array
            Array of lasso beta regression coefficients

        """
        lasso = LassoWrapper(window_values)
        beta_coef = lasso.get_coeffs(alpha)
        return beta_coef

    def get_null_alpha(self):
        pass

    def cross_validate_alpha(self):
        pass

    def fit_window(self):
        print "This method needs to override the super class"

    def permutation_test(self):
        print "This method needs to override the super class"

    def bootstrap(self):
        print "This method needs to override the super class"

    def initialize_params(self):
        print "This method needs to override the super class"
