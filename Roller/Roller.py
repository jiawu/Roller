import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
import numpy as np
from util.linear_wrapper import LassoWrapper

import pdb
class Roller:
    """
    A thing that grabs different timepoints of data, can set window and step size.

    To do list:
        -i need to make two modules: a file processing module and a rolling module.
        -accept different table formats

        -add permute_window()
        -add bootstrape_window()
    """
    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator = "\t"):
        """
        Initialize the roller object. Read the file and put it into a pandas dataframe
        :param file_path: file-like object or string
                        The file to read
        :param gene_start: int
        :param gene_end: int
        """
        # Read the raw data into a pandas dataframe object
        self.raw_data = pd.read_csv(file_path, sep=separator)
        self.raw_data = self.raw_data.dropna(axis=0, how='all')

        # Set roller defaults
        self.current_step = 0
        self.window_width = 3
        self.step_size = 1
        self.time_label = time_label
        #pdb.set_trace()

        # Get overall width of the time-course
        self.time_vec = self.raw_data[self.time_label].unique()
        self.overall_width = len(self.time_vec)

        if gene_end is not None:
            self.gene_end = gene_end
        else:
            self.gene_end = len(self.raw_data.columns)
        if gene_start is not None:
            self.gene_start = gene_start
        else:
            self.gene_start = 0

        self.gene_list = self.raw_data.columns.values[self.gene_start:self.gene_end]

        self.current_window = self.get_window()

    def get_n_windows(self):
        total_windows = (self.overall_width - self.window_width+1)/(self.step_size)
        return total_windows

    def get_window(self):
        raw_window = self.get_window_raw()
        only_genes = raw_window.iloc[:, self.gene_start:self.gene_end]
        return only_genes

    def get_window_raw(self):
        start_index = self.current_step
        end_index = start_index + self.window_width
        time_window = self.time_vec[start_index:end_index]
        data = self.raw_data[self.raw_data[self.time_label].isin(time_window)]
        return data

    def next(self):
        end_index = self.current_step + self.window_width
        if end_index <= self.overall_width:
            self.current_step += self.step_size
            self.current_window = self.get_window()
            return self.current_window
        else:
            return "end"

    def set_window(self, width):
        self.window_width = width

    def set_step(self, step):
        self.step_size = step

    def reset(self):
        self.current_step = 0

    def remove_blank_rows(self):
        """calculates sum of rows. if sum is NAN, then remove row"""
        coln = len(self.raw_data.columns)
        sums = [self.raw_data.iloc[:,x].sum() for x in range(0,coln)]
        ind = np.where(np.isnan(sums))[0]
        self.raw_data.iloc[:,ind]=0

    def get_n_genes(self):
        return(len(self.raw_data.columns) -1)

    def fit(self, window_size, method='lasso', alpha=0.2, resamples=0, noise=0.2):
        """
        Fit the rolling model

        :return: 3D matrix (n_windows, n_parents, n_children)
            Matrix of coefficients. Each slice is the beta coefficient matrix for a given window
        """

        # Set the window size
        self.set_window(window_size)

        # Calculate total number of windows, and regressors, then initialize coefficient matrix
        total_window_number = self.get_n_windows()
        n_genes = len(self.gene_list)

        if not resamples:
            coeff_matrix_3d = np.empty((n_genes, n_genes, total_window_number))

            for nth_window in range(total_window_number):
                #loop gets the window, gets the coefficients for that window, then increments the window
                current_window = self.get_window()

                filled_matrix = current_window.values
                #filled_matrix = imputer.fit_transform(current_window)
                current_lasso = LassoWrapper(filled_matrix)
                coeff_mat = current_lasso.get_coeffs(alpha)
                coeff_matrix_3d[:, :, nth_window] = coeff_mat
                #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
                self.next()

            return coeff_matrix_3d

        else:
            coeff_matrix_4d = np.empty((n_genes, n_genes, total_window_number, resamples))

            for nth_window in range(total_window_number):
                #loop gets the window, gets the coefficients for that window, then increments the window
                current_window = self.get_window()
                filled_matrix = current_window.values
                for sample in range(resamples):
                    sample_matrix = self.resample_window(filled_matrix)
                    noisy_sample = self.add_noise_to_window(sample_matrix, noise)
                    current_lasso = LassoWrapper(noisy_sample)
                    coeff_mat = current_lasso.get_coeffs(alpha)
                    coeff_matrix_4d[:, :, nth_window, sample] = coeff_mat
                self.next()

            return coeff_matrix_4d

    def resample_window(self, window_values):
        """
        Resample window values, along a specific axis
        :param window_values: array

        :return: array
        """
        n, p = window_values.shape

        # For each column randomly choose samples
        resample_values = np.array([np.random.choice(window_values[:,ii], size=n) for ii in range(p)]).T

        return resample_values

    def add_noise_to_window(self, window_values, max_random=0.2):
        """
        Add uniform noise to each value
        :param window_values: array

        :param max_random: float
            Amount of noise to add to each value, plus or minus
        :return: array

        """
        noise = np.random.uniform(low=1-max_random, high=1+max_random, size=window_values.shape)
        noisy_values = np.multiply(window_values, noise)
        return noisy_values
