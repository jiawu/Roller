__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import warnings
import numpy as np
import pandas as pd
from util.linear_wrapper import LassoWrapper
from sklearn import linear_model
from sklearn.cross_validation import KFold


class Window(object):
    def __init__(self, dataframe):
        self.df = dataframe
        self.window_values = dataframe.values
        self.samples = dataframe.index.values
        self.genes = dataframe.columns.values
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
        """
        Run the permutation test and update the edge_table with p values. It is expected that this method will vary
        depending on the type of method used by the window
        :return:
        """
        pass

    def bootstrap(self):
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

class Lasso_Window(Window):
    def __init__(self, dataframe):
        super(Lasso_Window, self).__init__(dataframe)
        self.alpha = None
        self.null_alpha = None
        self.beta_coefficients = None
        self.cv_alpha_range = None
        self.cv_scores = None

    def permutation_test(self):
        print "This method needs to override the super class"

    def bootstrap(self):
        print "This method needs to override the super class"

    def initialize_params(self, alpha=None):
        """
        Choose the value of alpha to use for fitting
        :param alpha: float, optional
            The alpha value to use for the window. If none is entered the alpha will be chosen by cross validation
        :return:
        """
        if alpha is None:
            self.alpha = self.cross_validate_alpha()

    def fit_window(self, alpha=0.2):
        """returns a 2D array with target as rows and regulators as columns"""
        clf = linear_model.Lasso(alpha)
        #loop that iterates through the target genes
        all_data = self.window_values
        coeff_matrix = np.array([],dtype=np.float_).reshape(0, all_data.shape[1])

        for col_index,column in enumerate(all_data.T):
            #delete the column that is currently being tested
            X_matrix = np.delete(all_data, col_index, axis=1)
            #take out the column so that the gene does not regress on itself
            target_TF = all_data[:,col_index]
            clf.fit(X_matrix, target_TF)
            coeffs = clf.coef_
            #artificially add a 0 to where the col_index is
            #to prevent self-edges
            coeffs = np.insert(coeffs,col_index,0)
            coeff_matrix=np.vstack((coeff_matrix,coeffs))
        return coeff_matrix

    def get_null_alpha(self, max_expected_alpha=1e4, min_step_size=1e-9):
        """
        Get the smallest value of alpha that returns a lasso coefficient matrix of all zeros

        :param max_expected_alpha: float, optional

            Largest expected value of alpha that will return all zeros. This is a guess and is dependent on the data.
            The function will step from the minimum to this value with a step size one power of 10 less than this value

        :param min_step_size: float, optional

            The smallest step size that will be taken. This also defines the precision of the alpha value returned

        :return: float

            The smallest alpha value that will return all zero beta values, subject to the precision of min_step_size

        """
        warnings.simplefilter("ignore")
        # Get maximum edges, assuming all explanors are also response variables and no self edges
        [n, p] = self.window_values.shape
        max_edges = p * (p-1)

        # Raise exception if Lasso doesn't converge with alpha == 0
        if np.count_nonzero(self.fit_window(0)) != max_edges:
            raise ValueError('Lasso does not converge with alpha = 0')

        # Raise exception if max_expected_alpha does not return all zero betas
        if np.count_nonzero(self.fit_window(max_expected_alpha)) != 0:
            raise ValueError('max_expected_alpha not high enough, coefficients still exist. Guess higher')

        # Set ranges of step sizes, assumed to be powers of 10
        powers = int(np.log10(max_expected_alpha/min_step_size))
        step_sizes = [max_expected_alpha/(10**ii) for ii in range(powers+1)]

        # Intialize loop values
        cur_min = 0
        alpha_max = step_sizes[0]

        # Start stepping with forward like selection
        for ii, cur_max in enumerate(step_sizes[:-1]):

            # Set the maximum for the range to scan
            if alpha_max > cur_max:
                cur_max = alpha_max

            # Set the current step size and new range to look through
            cur_step = step_sizes[ii+1]
            cur_range = np.linspace(cur_min, cur_max, (cur_max-cur_min)/cur_step+1)

            # In the current range, check when coefficients start popping up
            for cur_alpha in cur_range:
                num_coef = np.count_nonzero(self.fit_window(cur_alpha))
                if num_coef > 0:
                    cur_min = cur_alpha
                elif num_coef == 0:
                    # Found a new maximum that eliminates all betas, no need to keep stepping
                    alpha_max = cur_alpha
                    break
        return alpha_max

    def cv_select_alpha(self, alpha_range, n_folds=3):
        self.cv_alpha_range = alpha_range
        self.cv_scores = [self.cross_validate_alpha(alpha, n_folds) for alpha in alpha_range]
        #print self.cv_scores

    def cross_validate_alpha(self, alpha, n_folds):
        '''
        Get a Q^2 value for each explanatory value (column) at the given alpha value
        :param alpha: float
        :param n_folds: int
            when number of folds is the same as number of samples this is equivalent to leave-one-out
        :return:
        '''
        data = self.window_values.copy()
        n_elements = len(data)
        kf = KFold(n_elements, n_folds)

        press = 0.0
        ss = 0.0

        for train_index, test_index in kf:
            x_train = data[train_index]
            x_test = data[test_index]
            y_test = x_test.copy()

            # Run Lasso
            current_coef = self.get_coeffs(x_train, alpha)

            y_predicted = np.dot(x_test, current_coef)

            # Calculate PRESS and SS
            current_press = np.sum(np.power(y_predicted-y_test, 2), axis=0)
            current_ss = self.sum_of_squares(y_test)

            press += current_press
            ss += current_ss
        q_squared = 1-press/ss
        return q_squared

    def sum_of_squares(self, X):
        """
        Calculate the sum of the squares for each column
        :param X: array-like
            The data matrix for which the sum of squares is taken
        :return: float or array-like
            The sum of squares, columnwise or total
        """
        column_mean = np.mean(X, axis=0)
        ss = np.sum(np.power(X-column_mean, 2), axis=0)
        return ss

    def get_coeffs(self, data, alpha):
        """returns a 2D array with target as rows and regulators as columns"""
        clf = linear_model.Lasso(alpha)
        #loop that iterates through the target genes
        all_data = data.copy()
        coeff_matrix = np.array([],dtype=np.float_).reshape(0, all_data.shape[1])

        for col_index,column in enumerate(all_data.T):
            #delete the column that is currently being tested
            X_matrix = np.delete(all_data, col_index, axis=1)
            #take out the column so that the gene does not regress on itself
            target_TF = all_data[:,col_index]
            clf.fit(X_matrix, target_TF)
            coeffs = clf.coef_
            #artificially add a 0 to where the col_index is
            #to prevent self-edges
            coeffs = np.insert(coeffs,col_index,0)
            coeff_matrix=np.vstack((coeff_matrix,coeffs))

        return coeff_matrix