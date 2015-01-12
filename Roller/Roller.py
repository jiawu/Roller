import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
import numpy as np
from util.linear_wrapper import LassoWrapper


class Roller(object):
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

    def fit_window(self, window, alpha):
        """
        Given a window get the lasso coefficients
        :param window: data-frame
            A roller window
        :param alpha: float
            Value to use for lasso regression
        :return: array
            Array of lasso beta regression coefficients

        """
        window_values = window.values
        lasso = LassoWrapper(window_values)
        beta_coef = lasso.get_coeffs(alpha)
        return beta_coef

    def fit_model(self, window_size, method='lasso', alpha=0.2, resamples=0, noise=0.2):
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
                coeff_mat = self.fit_window(current_window, alpha)
                coeff_matrix_3d[:, :, nth_window] = coeff_mat
                #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
                self.next()

            return coeff_matrix_3d

        else:
            coeff_matrix_4d = np.empty((n_genes, n_genes, total_window_number, resamples))
            for nth_window in range(total_window_number):
                #loop gets the window, gets the coefficients for that window, then increments the window
                current_window = self.get_window()

                for sample in range(resamples):
                    sample_window = self.resample_window(current_window)
                    noisy_window = self.add_noise_to_window(sample_window, noise)
                    coeff_mat = self.fit_window(noisy_window, alpha)
                    coeff_matrix_4d[:, :, nth_window, sample] = coeff_mat
                self.next()

            return coeff_matrix_4d

    def resample_window(self, window):
        """
        Resample window values, along a specific axis
        :param window_values: array

        :return: array
        """
        window_values = window.values
        n, p = window_values.shape

        # For each column randomly choose samples
        resample_values = np.array([np.random.choice(window_values[:, ii], size=n) for ii in range(p)]).T

        resample_window = pd.DataFrame(resample_values, columns=window.columns.values.copy(),
                                       index=window.index.values.copy())

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

    def zscore_all_data(self):
        #zscores all the data
        dataframe =  self.raw_data

        #for each column, calculate the zscore
        #zscore is calculated as X - meanX / std(ddof = 1)
        for item in dataframe.columns:
            if item != self.time_label:
                dataframe[item] = (dataframe[item] - dataframe[item].mean())/dataframe[item].std(ddof=1)
        self.raw_data = dataframe

    def get_window_stats(self):
        """for each window, get a dict:
            N : the number of datapoints in this window,
            time_labels: the names of the timepoints in a roller model
            step_size: the step-size of the current model
            window_size: the size of the window of the current model
            total_windows: the number of windows total
            window_index: the index of the window. counts start at 0. ie if the window index is 0 it is the 1st window. if the window index is 12, it is the 12th window in the series."""
        current_window = self.get_window_raw()

        """calculate the window index. todo: move into own function later"""
        min_time = np.amin(current_window[self.time_label])
        window_index = np.where(self.time_vec == min_time)/self.step_size
        # to calculate the nth window, time vector
        # index of the time-vector, step size of 2? window 4, step size 2
        #
        #total windows = total width (10) - window_width (2) +1 / step size
        # 10 time points 0 1 2 3 4 5 6 7 8 9
        #width is 2: 0 and 1
        # step size is 2
        # 01, 12, 23, 34, 45, 56, 67, 78, 89

        #todo: so the issue is that total windows (get n windows) is the true number of windows, and window index is the nth -1 window... it would be great to consolidate these concepts but no big deal if they can't be.


        window_stats = {'N': len(current_window.index),
                        'time_labels': current_window[self.time_label].unique(),
                        'step_size': self.step_size,
                        'window_size': self.window_width,
                        'total_windows': self.get_n_windows(),
                        'window_index': window_index}
        return window_stats

    def cross_validate_window(self, window, n_alphas=50, n_folds=3):
        lasso = LassoWrapper(window.values)
        alpha_range = np.linspace(0, lasso.get_max_alpha(), n_alphas)
        q_squared_array = np.array([lasso.cross_validate_alpha(alpha, n_folds) for alpha in alpha_range])
        return alpha_range, q_squared_array

    def select_alpha(self, window, method='max', alpha_list=None):
        """
        Select the alpha values for each window to use for fitting the initial model

        :param window:
        :param method: str
            currently supports methods 'max' and 'manual'. 'max' selects alpha values that maximize the cross validation
            score for each window. 'manual' is just a dummy method that will pass through the list provided for it after
            verifying that is the right size.
        :param alpha_list:
        :return:
        """
        total_window_number = self.get_n_windows()
        #todo: determine where alpha overrides should be
        if method is 'max':
            alpha_range, q_list = self.cross_validate_window(window)
            alpha_table = zip(alpha_range, q_list)
            alpha_list.append(alpha_table)
            if pd['override-alpha'] == "false":
                (best_alpha,Qs) = max(alpha_table, key = lambda t: t[1])
            elif pd['override-alpha'] == "true":
                best_alpha = alpha
            coeff_mat = current_lasso.get_coeffs(best_alpha)


