__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import numpy as np
import pandas as pd
import itertools


class Window(object):
    """
    A window object is a snapshot created from the full data. Windows operate independently from each other. Windows
    should be sub-classed to have network inference specific features and methods.
    """

    def __init__(self, raw_dataframe, window_info):
        """
        Initialize a window object. Extract information from the passed data-frame. Generate edge list.

        :param raw_dataframe: pandas data-frame
        :param window_info: dict
            Dictionary that provides data that can be used to uniquely identify a window
        :return:
        """
        self.time_label = window_info['time_label']
        self.gene_start = window_info['gene_start']
        self.gene_end = window_info['gene_end']
        self.nth_window = window_info['nth_window']
        self.raw_data = raw_dataframe
        dataframe = raw_dataframe.iloc[:, self.gene_start:self.gene_end]
        self.df = dataframe
        self.window_values = dataframe.values
        self.samples = dataframe.index.values
        self.n_samples = len(self.samples)
        self.genes = dataframe.columns.values
        self.n_genes = len(self.genes)
        self.results_table = pd.DataFrame()
        self.edge_labels = [x for x in itertools.product(self.genes, repeat=2)]

        self.edge_list = self.possible_edge_list(self.genes, self.genes)
        # Add edge list to edge table
        self.results_table['regulator-target'] = self.edge_list

    def create_linked_list(self, numpy_array_2D, value_label):
        """labels and array should be in row-major order"""
        linked_list = pd.DataFrame({'regulator-target': self.edge_labels, value_label: numpy_array_2D.flatten()})
        return linked_list

    def aggregate(self):
        pass

    def get_window_stats(self):
        """for each window, get a dict:
            N : the number of datapoints in this window,
            time_labels: the names of the timepoints in a roller model
            step_size: the step-size of the current model
            window_size: the size of the window of the current model
            total_windows: the number of windows total
            window_index: the index of the window. counts start at 0. ie if the window index is 0 it is the 1st window. if the window index is 12, it is the 12th window in the series."""

        window_stats = {'n_samples': self.n_samples,
                        'n_genes': self.n_genes}
        return window_stats


    def possible_edge_list(self, parents, children, self_edges=True):
        """
        Create a list of all the possible edges between parents and children

        :param parents: array
            labels for parents
        :param children: array
            labels for children
        :param self_edges:
        :return: array, length = parents * children
            array of parent, child combinations for all possible edges
        """
        parent_index = range(len(parents))
        child_index = range(len(children))
        a, b = np.meshgrid(parent_index, child_index)
        parent_list = parents[a.flatten()]
        child_list = children[b.flatten()]
        possible_edge_list = None
        if self_edges:
            possible_edge_list = zip(parent_list, child_list)

        elif not self_edges:
            possible_edge_list = zip(parent_list[parent_list != child_list], child_list[parent_list != child_list])

        return possible_edge_list

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

        # resample_window = pd.DataFrame(resample_values, columns=self.df.columns.values.copy(),
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

        noise = np.random.uniform(low=1 - max_random, high=1 + max_random, size=window_values.shape)
        noisy_values = np.multiply(window_values, noise)
        return noisy_values

    def run_permutation_test(self):
        """
        Run the permutation test and update the edge_table with p values. It is expected that this method will vary
        depending on the type of method used by the window
        :return:
        """
        pass

    def run_bootstrap(self, n_bootstraps):
        """
        Run bootstrapping and update the edge_table with stability values. It is expected that this method will vary
        depending on the type of method used by the window
        :return:
        """
        pass

    def rank_edges(self, method):
        """
        Rank the edges in the edge table. This may eventually be window type specific.

        :param method: The method to use for ranking the edges in edge_table
        :return: list of tuples
            Sorted list [(regulator1, target1), (regulator2, target2)...] that can be scored with aupr or auroc
        """
        pass

    def initialize_params(self):
        """
        Initialize a model for the window and calculate the necessary parameters
        :return:
        """
        pass

    def get_coeffs(self, *args):
        """
        Get the beta coefficients
        :param args:
        :return:
        """
        pass

    def permute_data(self, array):
        """Warning: Modifies data in place. also remember the """
        new_array = array.copy()
        _ = [np.random.shuffle(i) for i in new_array]
        return new_array

    def update_variance_1D(self, prev_result, new_samples):
        """
        incremental calculation of means: accepts new_samples, which is a list of samples. then calculates a new mean.
        this is a useful function for calculating the means of large arrays
        """

        n = float(prev_result["n"])
        mean = float(prev_result["mean"])
        sum_squares = float(prev_result["ss"])

        for x in new_samples:
            n = n + 1
            # delta = float(x) - mean
            old_mean = mean
            mean = old_mean + (float(x) - old_mean) / n
            sum_squares = sum_squares + (float(x) - mean) * (float(x) - old_mean)

        if (n < 2):
            return 0

        variance = sum_squares / (n - 1)
        result = {"mean": mean,
                  "ss": sum_squares,
                  "variance": variance,
                  "n": n}
        return result

    def update_variance_2D(self, prev_result, new_samples):
        """incremental calculation of means: accepts new_samples, which is a list of samples.
        then calculates a new mean. this is a useful function for calculating the means of large arrays"""
        n = prev_result["n"]  # 2D numpy array with all zeros or watev
        mean = prev_result["mean"]  # 2D numpy array
        sum_squares = prev_result["ss"]  # 2D numpy array

        # new_samples is a list of arrays
        #x is a 2D array
        for x in new_samples:
            n += 1
            #delta = float(x) - mean
            old_mean = mean.copy()
            mean = old_mean + np.divide((x - old_mean), n)
            sum_squares = sum_squares + np.multiply((x - mean), (x - old_mean))

        if n[0, 0] < 2:
            result = {"mean": mean,
                      "ss": sum_squares,
                      "n": n}
            return result

        variance = np.divide(sum_squares, (n - 1))
        result = {"mean": mean,
                  "ss": sum_squares,
                  "variance": variance,
                  "n": n}
        return result
