import warnings
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.cross_validation import KFold
from scipy import integrate
import scipy
from Window import Window

class LassoWindow(Window):
    def __init__(self, dataframe, window_info):
        super(LassoWindow, self).__init__(dataframe, window_info)
        self.alpha = None
        self.beta_coefficients = None
        self.cv_table = None
        self.bootstrap_matrix = None
        self.freq_matrix = None
        self.edge_stability_auc = None
        self.permutation_means = None
        self.permutation_sd = None
        self.permutation_pvalues = None

        # Try set the null alpha value using default parameters.
        try:
            self.null_alpha = self.get_null_alpha()
        except ValueError:
            warnings.warn("Could not set null_alpha with default parameters. Set manually")
            self.null_alpha = None

    def generate_results_table(self):

        #generate edges for initial model
        initial_edges = self.create_linked_list(self.beta_coefficients, 'B')
        #permutation edges
        permutation_mean_edges =self.create_linked_list(self.permutation_means, 'p-means')
        permutation_sd_edges = self.create_linked_list(self.permutation_sd, 'p-sd')
        stability_edges = self.create_linked_list(self.edge_stability_auc, 'stability')

        aggregated_edges = initial_edges.merge(permutation_mean_edges, on='regulator-target').merge(permutation_sd_edges, on='regulator-target').merge(stability_edges, on='regulator-target')

        #sorry, it is a little messy to do the p-value calculations for permutation tests here...
        valid_indices = aggregated_edges['B'] != 0
        valid_window = aggregated_edges[valid_indices]
        initial_B = valid_window['B']
        sd = valid_window['p-sd']
        mean = valid_window['p-means']
        valid_window['final-z-scores-perm']=(initial_B - mean)/sd
        valid_window['cdf-perm'] = (-1*abs(valid_window['final-z-scores-perm'])).apply(scipy.stats.norm.cdf)
        #calculate t-tailed pvalue
        valid_window['p-value-perm'] = (2*valid_window['cdf-perm'])
        self.results_table = valid_window
        return(self.results_table)

    def rank_results(self, rank_by):
        rank_column_name = rank_by + "-rank"
        self.results_table[rank_column_name] = self.results_table[rank_by].rank(method="dense",ascending = True)
        self.results_table = self.results_table.sort(columns=rank_column_name, axis = 0)

        return(self.results_table)

    def permutation_test(self, permutation_n=1000):
        #initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        nth_window = self.nth_window
        zeros = np.zeros((self.n_genes, self.n_genes))
        #initialize running calculation
        result = {'n':zeros.copy(), 'mean':zeros.copy(), 'ss':zeros.copy()}
        #inner loop: permute the window N number of times
        permuted_window = self.df.copy()

        for nth_perm in range(0, permutation_n):
            if (nth_perm % 200 == 0):
                print('Perm Window: '+ str(nth_window) + ' Perm Run: ' +str(nth_perm))
            #permute data
            self.permute_data(permuted_window.values)
            #fit the data and get coefficients
            permuted_coeffs = self.get_coeffs(alpha=self.alpha, data=permuted_window.values)
            dummy_list = []
            dummy_list.append(permuted_coeffs)
            result = self.update_variance_2D(result,dummy_list)

        self.permutation_means = result['mean'].copy()
        self.permutation_sd= result['variance'].copy()

    def update_variance_2D(self, prev_result, new_samples):
        """incremental calculation of means: accepts new_samples, which is a list of samples. then calculates a new mean. this is a useful function for calculating the means of large arrays"""
        n = prev_result["n"] #2D numpy array with all zeros or watev
        mean = prev_result["mean"] #2D numpy array
        sum_squares = prev_result["ss"] #2D numpy array

        #new_samples is a list of arrays
        #x is a 2D array
        for x in new_samples:
            n = n + 1
            #delta = float(x) - mean
            old_mean = mean.copy()
            mean = old_mean + np.divide( (x-old_mean) , n)
            sum_squares = sum_squares + np.multiply((x-mean),(x-old_mean))

        if (n[0,0] < 2):
            result = {  "mean": mean,
                        "ss": sum_squares,
                        "n": n}
            return result

        variance = np.divide(sum_squares,(n-1))
        result = {  "mean": mean,
                    "ss": sum_squares,
                    "variance": variance,
                    "n": n}
        return result
    def run_bootstrap(self, n_bootstraps=1000, n_alphas=20, noise=0.2):
        alpha_range = np.linspace(0, self.null_alpha, n_alphas)
        self.bootstrap_matrix = np.empty((self.n_genes, self.n_genes, n_bootstraps, n_alphas))
        for ii, alpha in enumerate(alpha_range):
            self.bootstrap_matrix[:, :, :, ii] = self.bootstrap_alpha(alpha, n_bootstraps, noise)

        self.freq_matrix = self.calc_edge_freq()/float(n_bootstraps)

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
                boot_matrix = self.get_coeffs(alpha, noisy_window)
            else:
                boot_matrix = np.dstack((boot_matrix, self.get_coeffs(alpha, noisy_window)))

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
        auc = self.edge_stability_auc[:,:, nth]
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

    def fit_window(self):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        # Make sure an alpha value has been selected to use for fitting the window
        if self.alpha is None:
            raise ValueError("window alpha value must be set before the window can be fit")

        self.beta_coefficients = self.get_coeffs(self.alpha)

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
        if np.count_nonzero(self.get_coeffs(0)) != max_edges:
            raise ValueError('Lasso does not converge with alpha = 0')

        # Raise exception if max_expected_alpha does not return all zero betas
        if np.count_nonzero(self.get_coeffs(max_expected_alpha)) != 0:
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
                num_coef = np.count_nonzero(self.get_coeffs(cur_alpha))
                if num_coef > 0:
                    cur_min = cur_alpha
                elif num_coef == 0:
                    # Found a new maximum that eliminates all betas, no need to keep stepping
                    alpha_max = cur_alpha
                    break
        return alpha_max

    def cv_select_alpha(self, alpha_range=None, method='modelQ2', n_folds=3):
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

        if alpha_range is None:
            alpha_range = np.linspace(0, self.null_alpha)

        # Calculate the cv_table values
        cv_results = np.array([self.cross_validate_alpha(alpha, n_folds, True) for alpha in alpha_range])
        column_labels = np.append(self.genes+"_Q^2", ["Model_Q^2", "positive_q2"])

        # Put the results into a dataframe
        df = pd.DataFrame(cv_results, columns=column_labels)

        # Insert the alpha range as the first column
        df.insert(0, 'alpha', alpha_range)

        if method == 'modelQ2':
            # Set the window alpha to the alpha that produced the highest Q^2 value
            cv_selected_alpha = alpha_range[df["Model_Q^2"].idxmax(1)]

        elif method == 'max_posQ2':
            #todo: an secondary criteria needs to be employed if this is used because there will likely be multiple alphas with the same number of positive Q2
            cv_selected_alpha = alpha_range[df["positive_q2"].idxmax(1)]

        else:
            raise ValueError("User entered method %s is not valid" %method)

        return cv_selected_alpha, df.copy()

    def cross_validate_alpha(self, alpha, n_folds, condensed=False):
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
            current_coef = self.get_coeffs(alpha, x_train)

            y_predicted = np.dot(x_test, current_coef)

            # Calculate PRESS and SS
            current_press = np.sum(np.power(y_predicted-y_test, 2), axis=0)
            current_ss = self.sum_of_squares(y_test)

            press += current_press
            ss += current_ss

        # Calculate the Q^2 value of each gene
        q_squared = 1-press/ss

        # Calculate the Q^2 of the whole model. This is different than averaging the individual q_squared values
        model_q_squared = 1 - np.sum(press)/np.sum(ss)

        if condensed:
            return np.append(q_squared, [model_q_squared, np.sum(q_squared > 0)])
        else:
            return q_squared, model_q_squared

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

    def get_coeffs(self, alpha, data=None):
        """
        returns a 2D array with target as rows and regulators as columns
        :param alpha: float
            value to use for lasso fitting
        :param data: array-like, optional
            Data to fit. If none, will use the window values. Default is None
        :return:
        """
        clf = linear_model.Lasso(alpha)
        #loop that iterates through the target genes
        if data is None:
            all_data = self.window_values.copy()
        else:
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
