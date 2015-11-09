import warnings
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.cross_validation import KFold
from scipy import integrate
from scipy import stats
from sklearn.metrics import mean_squared_error
from util.utility_module import sum_of_squares
import scipy
import sys

from Window import Window


class LassoWindow(Window):
    """
    A window that runs Lasso regression as the network inference method
    """

    def __init__(self, dataframe, window_info, roller_data, td_window, explanatory_dict, response_dict):
        #todo: unclear if roller_data is necessary
        """
        Initialize a LassoWindow with necessary data
        :param dataframe: pandas data-frame
        :param window_info: dict
            Dictionary contain the parameters of the window
        :param roller_data:
        :return:
        """
        super(LassoWindow, self).__init__(dataframe, window_info, roller_data, td_window, explanatory_dict,
                                          response_dict)

        # Initialize window attributes
        self.alpha = None
        self.null_alpha = None
        self.beta_coefficients = None
        self.cv_table = None
        self.bootstrap_matrix = None
        self.freq_matrix = None
        self.edge_stability_auc = None

    def sort_edges(self, method="p_value"):
        """
        Sort the edge table based on the selected edge ranking method
        :param method: str
            Which method to use for ranking the edges
        :return:

        Called by:
            Swing.average_rank()
        """
        if self.edge_table is None:
            raise ValueError("The edge table must be created before getting edges")
        temp_edge_table = self.edge_table.copy()
        if method == "p_value":
            temp_edge_table.sort(columns=['p_value', 'stability'], ascending=[True, False], inplace=True)
        elif method == "stability":
            temp_edge_table.sort(columns=['stability', 'p_value'], ascending=[False, True], inplace=True)

        return temp_edge_table['regulator-target'].values

    def rank_results(self, rank_by, ascending=False):
        rank_column_name = rank_by + "-rank"
        ##rank edges with an actual beta value first until further notice ##
        valid_indices = self.results_table['B'] != 0
        valid_window = self.results_table[valid_indices]
        valid_window[rank_column_name] = valid_window[rank_by].rank(method="dense", ascending=ascending)
        edge_n = len(valid_window.index)

        invalid_indices = self.results_table['B'] == 0
        invalid_window = self.results_table[invalid_indices]
        invalid_window[rank_column_name] = invalid_window[rank_by].rank(method="dense", ascending=ascending)
        invalid_window[rank_column_name] = invalid_window[rank_column_name] + edge_n
        self.results_table = valid_window.append(invalid_window)
        self.results_table = self.results_table.sort(columns=rank_column_name, axis=0)

        return self.results_table

    def _permute_coeffs(self, zeros, crag=False, n_permutations=10):
        result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}
        # inner loop: permute the window N number of times
        for nth_perm in range(0, n_permutations):
            # if (nth_perm % 200 == 0):
            # print 'Perm Run: ' +str(nth_perm)

            # permute data
            permuted_data = self.permute_data(self.explanatory_data)

            # fit the data and get coefficients
            permuted_coeffs, _ = self.get_coeffs(self.alpha, x_data=permuted_data, crag=crag)
            dummy_list = [permuted_coeffs]
            result = self.update_variance_2D(result, dummy_list)

        self.permutation_means = result['mean'].copy()
        self.permutation_sd = np.sqrt(result['variance'].copy())
        self.permutation_p_values = self.calc_p_value()

    def run_permutation_test(self, n_permutations=10, crag=False):
        # initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        zeros = np.zeros(self.beta_coefficients.shape)

        # initialize running calculation
        self._permute_coeffs(zeros=zeros, n_permutations=n_permutations)

    def calc_p_value(self, value=None, mean=None, sd=None):
        if value is None:
            value = self.beta_coefficients.copy()
        if mean is None:
            mean = self.permutation_means.copy()
        if sd is None:
            sd = self.permutation_sd.copy()

        z_scores = (value - mean) / sd
        cdf = stats.norm.cdf((-1 * abs(z_scores)))
        p_values = 2 * cdf
        return p_values

    def run_bootstrap(self, n_bootstraps=10, n_alphas=20, noise=0.2):
        if self.null_alpha is None:
            # Try set the null alpha value using default parameters.
            try:
                self.null_alpha = self.get_null_alpha()
            except ValueError:
                warnings.warn("Could not set null_alpha with default parameters. Set manually")
        alpha_range = np.linspace(0, self.null_alpha, n_alphas)

        n_columns = self.explanatory_data.shape[1]

        self.bootstrap_matrix = np.empty((self.n_genes, n_columns, n_bootstraps, n_alphas))

        for ii, alpha in enumerate(alpha_range):
            self.bootstrap_matrix[:, :, :, ii] = self.bootstrap_alpha(alpha, n_bootstraps, noise)

        self.freq_matrix = self.calc_edge_freq() / float(n_bootstraps)

        self.edge_stability_auc = self.auc(self.freq_matrix, alpha_range)

    def bootstrap_alpha(self, alpha, resamples, noise):
        """
        Bootstrap a particular alpha, return a stack of beta matrices

        :param alpha:
        :param resamples:
        :param noise:
        :return:
        """
        boot_matrix = None
        for ii, sample in enumerate(range(resamples)):
            sample_window = self.resample_window()
            noisy_window = self.add_noise_to_values(sample_window, noise)
            if ii == 0:
                boot_matrix, _ = self.get_coeffs(alpha, x_data=noisy_window)
            else:
                current_coeff, _ = self.get_coeffs(alpha, x_data=noisy_window)
                boot_matrix = np.dstack((boot_matrix, current_coeff))

        return boot_matrix

    def calc_edge_freq(self):
        "This is agnostic to the edge sign, only whether it exists"
        edge_exists = self.bootstrap_matrix != 0
        freq_matrix = np.sum(edge_exists, axis=2)
        return freq_matrix

    def auc(self, matrix, xvalues, axis=-1):
        """
        Calculate area under the curve
        :return:
        """
        edge_auc = integrate.trapz(matrix, xvalues, axis=axis)
        return edge_auc

    def get_nth_window_auc(self, nth):
        auc = self.edge_stability_auc[:, :, nth]
        return auc

    def initialize_params(self, alpha=None):
        """
        Choose the value of alpha to use for fitting
        :param alpha: float, optional
            The alpha value to use for the window. If none is entered the alpha will be chosen by cross validation
        :return:
        """
        if alpha is None:
            "Select alpha with default parameters"
            self.alpha, self.cv_table = self.cv_select_alpha()
        elif alpha >= 0:
            self.alpha = alpha
        else:
            raise ValueError("alpha must be float (>=0) or None")
        return

    def fit_window(self, crag=False, calc_mse=False):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        # Make sure an alpha value has been selected to use for fitting the window
        if self.alpha is None:
            raise ValueError("window alpha value must be set before the window can be fit")

        self.beta_coefficients, self.edge_mse_diff = self.get_coeffs(self.alpha, crag=crag, calc_mse=calc_mse)

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

        child_nodes = self.response_data.shape[1]
        parent_nodes = self.explanatory_data.shape[1]
        if self.nth_window in self.explanatory_window:
            # Can't allow self edges
            max_edges = parent_nodes*child_nodes-child_nodes
        else:
            max_edges = parent_nodes*child_nodes

        # Raise exception if Lasso doesn't converge with alpha == 0
        zero_alpha_coeffs, _ = self.get_coeffs(0)
        if np.count_nonzero(zero_alpha_coeffs) != max_edges:
            raise ValueError('Lasso does not converge with alpha = 0')

        # Raise exception if max_expected_alpha does not return all zero betas
        high_alpha_coeffs, _ = self.get_coeffs(max_expected_alpha)
        if np.count_nonzero(high_alpha_coeffs) != 0:
            raise ValueError('max_expected_alpha not high enough, coefficients still exist. Guess higher')

        # Set ranges of step sizes, assumed to be powers of 10
        powers = int(np.log10(max_expected_alpha / min_step_size))
        step_sizes = [max_expected_alpha / (10 ** ii) for ii in range(powers + 1)]

        # Intialize loop values
        cur_min = 0
        alpha_max = step_sizes[0]

        # Start stepping with forward like selection
        for ii, cur_max in enumerate(step_sizes[:-1]):

            # Set the maximum for the range to scan
            if alpha_max > cur_max:
                cur_max = alpha_max

            # Set the current step size and new range to look through
            cur_step = step_sizes[ii + 1]
            cur_range = np.linspace(cur_min, cur_max, (cur_max - cur_min) / cur_step + 1)

            # In the current range, check when coefficients start popping up
            for cur_alpha in cur_range:
                cur_coeffs, _ = self.get_coeffs(cur_alpha)
                num_coef = np.count_nonzero(cur_coeffs)
                if num_coef > 0:
                    cur_min = cur_alpha
                elif num_coef == 0:
                    # Found a new maximum that eliminates all betas, no need to keep stepping
                    alpha_max = cur_alpha
                    break
        return alpha_max

    def cv_select_alpha(self, alpha_range=None, method='modelQ2', n_folds=5):
        """
        Calculate the cross-validation metrics for a range of alpha values and select the best one based on the chosen
        method
        :param alpha_range: list, optional
            A range of alpha values to test. Default is a linspace of 50 alpha values from 0.0 to null_alpha
        :param method: str, optional
            The method to use to set the alpha value. Current methods are 'modelQ2' and 'max_posQ2'. 'modelQ2' uses the
            alpha that has the highest Q-squared value for the whole model. 'max_posQ2' uses the alpha that has the
            maximum number of genes with positive Q-squared values. Default is modelQ2
        :param n_folds: int, optional
            The number of folds to use during cross validation. Default is 3. If n_folds equals the number of samples
            this is equivalent to leave-one-out-validation
        :return: tuple
            (alpha, cv_table)
        """
        if self.null_alpha is None:
            # Try set the null alpha value using default parameters.
            self.null_alpha = self.get_null_alpha()

        if alpha_range is None:
            alpha_range = np.linspace(0, self.null_alpha, num=25)

        # Calculate the cv_table values
        cv_results = np.array([self.cross_validate_alpha(alpha, n_folds, True) for alpha in alpha_range])
        column_labels = np.append(self.genes + "_Q^2", ["Model_Q^2", "positive_q2"])

        # Put the results into a dataframe
        df = pd.DataFrame(cv_results, columns=column_labels)

        # Insert the alpha range as the first column
        df.insert(0, 'alpha', alpha_range)

        if method == 'modelQ2':
            # Set the window alpha to the alpha that produced the highest Q^2 value
            cv_selected_alpha = alpha_range[df["Model_Q^2"].idxmax(1)]

        elif method == 'max_posQ2':
            # todo: an secondary criteria needs to be employed if this is used because there will likely be multiple alphas with the same number of positive Q2
            cv_selected_alpha = alpha_range[df["positive_q2"].idxmax(1)]

        else:
            raise ValueError("User entered method %s is not valid" % method)

        return cv_selected_alpha, df.copy()

    def cross_validate_alpha(self, alpha, n_folds, condensed=False):
        """
        Get a Q^2 value for each explanatory value (column) at the given alpha value
        :param alpha: float
        :param n_folds: int
            when number of folds is the same as number of samples this is equivalent to leave-one-out
        :return:
        """
        x_data = self.explanatory_data.copy()
        y_data = self.response_data.copy()
        n_elements = len(x_data)
        kf = KFold(n_elements, n_folds)

        press = 0.0
        ss = 0.0

        for train_index, test_index in kf:
            x_train = x_data[train_index]
            y_train = y_data[train_index]
            x_test = x_data[test_index]
            y_test = y_data[test_index]

            # Run Lasso
            current_coef, _ = self.get_coeffs(alpha, x_data=x_train, y_data=y_train)
            
            y_predicted = np.dot(x_test, current_coef.T)

            # Calculate PRESS and SS
            current_press = np.sum(np.power(y_predicted - y_test, 2), axis=0)
            current_ss = sum_of_squares(y_test)

            press += current_press
            ss += current_ss

        # Calculate the Q^2 value of each gene
        q_squared = 1 - press / ss

        # Calculate the Q^2 of the whole model. This is different than averaging the individual q_squared values
        model_q_squared = 1 - np.sum(press) / np.sum(ss)

        if condensed:
            return np.append(q_squared, [model_q_squared, np.sum(q_squared > 0)])
        else:
            return q_squared, model_q_squared

    def _fitstack_coeffs(self, alpha, coeff_matrix, model_list, x_matrix, target_y, col_index, crag=False):
        """
                                      
        example call:
        coeff_matrix, model_list = self._fitstack_coeffs(coeff_matrix, model_list, all_data, col_index, n_trees, n_jobs,crag)
        """
        # Initialize the model
        clf = linear_model.Lasso(alpha)

        # Fit the model
        clf.fit(x_matrix, target_y)

        model_params = {'col_index': col_index,
                        'response': target_y,
                        'predictor': x_matrix,
                        'model': clf}
        model_list.append(model_params)
        coeffs = clf.coef_

        # artificially add a 0 to where the col_index is
        # to prevent self-edges
        if coeff_matrix.shape[1] - len(coeffs) == 1:
            coeffs = np.insert(coeffs, col_index, 0)

        # Stack coefficients
        coeff_matrix = np.vstack((coeff_matrix, coeffs))

        if crag:
            training_scores, test_scores = self.crag_window(model_params)
            self.training_scores.append(training_scores)
            self.test_scores.append(test_scores)

        return coeff_matrix, model_list

    def get_coeffs(self, alpha, crag=False, x_data=None, y_data=None, calc_mse=False):
        """
        :param x_data:
        :param n_trees:
        :return: array-like
            An array in which the rows are children and the columns are the parents
        """
        # initialize items
        if y_data is None:
            y_data = self.response_data.copy()
        if x_data is None:
            x_data = self.explanatory_data.copy()

        y_labels = self.response_labels.copy()
        x_windows = self.explanatory_window.copy()
        x_labels = self.explanatory_labels.copy()
        coeff_matrix, model_list = self._initialize_coeffs(x_data)
        mse_matrix = None

        # Calculate a model for each target column
        for col_index, column in enumerate(y_data.T):
            target_y = column.copy()
            x_matrix = x_data.copy()
            insert_index = col_index

            # If the current window values are in the x_data, remove them
            if self.nth_window in x_windows:
                keep_columns = ~((x_windows == self.nth_window) & (x_labels == y_labels[col_index]))
                insert_index = list(keep_columns).index(False)
                x_matrix = x_matrix[:, keep_columns].copy()

            coeff_matrix, model_list = self._fitstack_coeffs(alpha, coeff_matrix, model_list, x_matrix, target_y,
                                                             insert_index, crag=crag)

            base_mse = mean_squared_error(model_list[col_index]['model'].predict(x_matrix), target_y)

            if calc_mse:
                f_coeff_matrix, f_model_list = self._initialize_coeffs(data=x_matrix)
                mse_list = []
                for idx in range(x_matrix.shape[1]):
                    adj_x_matrix = np.delete(x_matrix, idx, axis=1)
                    f_coeff_matrix, f_model_list = self._fitstack_coeffs(alpha, f_coeff_matrix, f_model_list,
                                                                         adj_x_matrix, target_y, idx, crag)
                    mse_diff = base_mse - mean_squared_error(f_model_list[idx]['model'].predict(adj_x_matrix), target_y)
                    mse_list.append(mse_diff)
                if mse_matrix is None:
                    mse_matrix = np.array(mse_list)
                else:
                    mse_matrix = np.vstack((mse_matrix, np.array(mse_list)))

        importance_dataframe = pd.DataFrame(coeff_matrix, index=y_labels, columns=x_labels)
        importance_dataframe.index.name = 'Child'
        importance_dataframe.columns.name = 'Parent'
        return importance_dataframe, mse_matrix
    
    def make_edge_table(self, calc_mse=False):
        """

        :return:

        Called by:
            Swing.rank_edges()
        """
        # Build indexing method for all possible edges. Length = number of parents * number of children
        parent_index = range(self.beta_coefficients.shape[1])
        child_index = range(self.beta_coefficients.shape[0])
        a, b = np.meshgrid(parent_index, child_index)

        # Flatten arrays to be used in link list creation
        df = pd.DataFrame()
        df['Parent'] = self.beta_coefficients.columns.values[a.flatten()]
        df['Child'] = self.beta_coefficients.index.values[b.flatten()]
        df['Importance'] = self.beta_coefficients.values.flatten()
        df['Stability'] = self.edge_stability_auc.flatten()
        df['P_window'] = self.explanatory_window[a.flatten()]

        # Calculate the window of the child node, which is equivalent to the current window index
        child_values = np.array([self.nth_window] * self.beta_coefficients.shape[0])
        df['C_window'] = child_values[b.flatten()]

        if self.permutation_p_values is not None:
            df["p_value"] = self.permutation_p_values.flatten()

        # Remove any self edges
        df = df[~((df['Parent'] == df['Child']) & (df['P_window'] == df['C_window']))]

        if calc_mse:
            df['MSE_diff'] = self.edge_mse_diff.flatten()

        return df


