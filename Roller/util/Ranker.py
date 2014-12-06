__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
from linear_wrapper import LassoWrapper
import numpy as np

class Bootstrapper(object):
    # Method specific
    # Generate n models
    # Keep beta value ranges for each model
    # Count how many times a beta value appears
    def __init__(self, roller_object):
        self.max_alpha = None
        self.bootstrap_matrix = None
        self.roller_object = roller_object

        # Set maximum alpha
        self.set_max_alpha()

    def run_bootstrap(self, window_size, n_bootstraps, n_alphas):
        alpha_range = np.linspace(0, self.max_alpha, n_alphas)
        self.roller_object.set_window(window_size)
        for ii, alpha in enumerate(alpha_range):
            current_coef = self.roller_object.fit(window_size, alpha=alpha, resamples=n_bootstraps)
            if ii is 0:
                empty_shape = list(current_coef.shape) + [len(alpha_range)]
                self.bootstrap_matrix = np.empty(tuple(empty_shape))
            self.bootstrap_matrix[:, :, :, :, ii] = current_coef

    def set_max_alpha(self):
        # Get maximum possible alpha for the whole data set
        self.roller_object.set_window(self.roller_object.overall_width)
        current_window = self.roller_object.get_window()
        lasso = LassoWrapper(current_window.values)
        self.max_alpha = lasso.get_max_alpha()
