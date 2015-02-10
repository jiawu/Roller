import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import KFold
from scipy import integrate
from scipy import stats
import scipy

from Window import Window

class RandomForestRegressionWindow(Window):
    def __init__(self, dataframe, window_info):
        super(RandomForestRegressionWindow, self).__init__(dataframe, window_info)
        self.beta_coefficients = None
        self.n_trees = None
        self.bootstrap_matrix = None
        self.freq_matrix = None
        self.edge_stability_auc = None
        self.permutation_means = None
        self.permutation_sd = None
        self.permutation_p_values = None
        self.permutation_pvalues = None

    def make_edge_table(self):
        if self.permutation_p_values is None:
            raise ValueError("p values must be set before making the edge table. Use method run_permutation test")

        if self.edge_stability_auc is None:
            raise ValueError("edge stability values must be set before making the edge table. "
                             "Use method run_permutation test")
        self.edge_table["P_Value"] = self.permutation_p_values.flatten()
        self.edge_table["Stability"] = self.edge_stability_auc.flatten()

    def rank_edges(self, method="p_value"):
        if self.edge_table is None:
            raise ValueError("The edge table must be created before getting edges")
        temp_edge_table = self.edge_table.copy()
        if method == "p_value":
            temp_edge_table.sort(columns=['P_Value', 'Stability'], ascending=[True, False], inplace=True)
        elif method == "stability":
            temp_edge_table.sort(columns=['Stability', 'P_Value'], ascending=[False, True], inplace=True)

        return temp_edge_table['Edge'].values

    def generate_results_table(self):

        #generate edges for initial model
        initial_edges = self.create_linked_list(self.beta_coefficients, 'B')
        #permutation edges
        permutation_mean_edges =self.create_linked_list(self.permutation_means, 'p-means')
        permutation_sd_edges = self.create_linked_list(self.permutation_sd, 'p-sd')
        stability_edges = self.create_linked_list(self.edge_stability_auc, 'stability')

        aggregated_edges = initial_edges.merge(permutation_mean_edges, on='regulator-target').merge(permutation_sd_edges, on='regulator-target').merge(stability_edges, on='regulator-target')

        #sorry, it is a little messy to do the p-value calculations for permutation tests here...
        valid_indices = aggregated_edges['p-sd'] != 0
        #valid_indices = aggregated_edges['B'] != 0
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

    def rank_results(self, rank_by, ascending=False):
        rank_column_name = rank_by + "-rank"
        ##rank edges with an actual beta value first until further notice ##
        valid_indices = self.results_table['B'] != 0
        valid_window = self.results_table[valid_indices]
        valid_window[rank_column_name] = valid_window[rank_by].rank(method="dense",ascending = ascending)
        edge_n=len(valid_window.index)

        invalid_indices = self.results_table['B'] == 0
        invalid_window = self.results_table[invalid_indices]
        invalid_window[rank_column_name] = invalid_window[rank_by].rank(method="dense",ascending = ascending)
        invalid_window[rank_column_name] = invalid_window[rank_column_name] + edge_n
        self.results_table = valid_window.append(invalid_window)
        self.results_table = self.results_table.sort(columns=rank_column_name, axis = 0)

        return(self.results_table)

    def run_permutation_test(self, n_permutations=1000):
        #initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        nth_window = self.nth_window
        zeros = np.zeros((self.n_genes, self.n_genes))
        #initialize running calculation
        result = {'n':zeros.copy(), 'mean':zeros.copy(), 'ss':zeros.copy()}
        #inner loop: permute the window N number of times
        for nth_perm in range(0, n_permutations):
            #if (nth_perm % 200 == 0):
                #print 'Perm Run: ' +str(nth_perm)

            #permute data
            permuted_data = self.permute_data(self.window_values)

            #fit the data and get coefficients

            permuted_coeffs = self.get_coeffs(self.n_trees, permuted_data)
            dummy_list = []
            dummy_list.append(permuted_coeffs)
            result = self.update_variance_2D(result, dummy_list)

        self.permutation_means = result['mean'].copy()
        self.permutation_sd = np.sqrt(result['variance'].copy())
        self.permutation_p_values = self.calc_p_value()

    def calc_p_value(self, value=None, mean=None, sd=None):
        if value is None:
            value = self.beta_coefficients.copy()
        if mean is None:
            mean = self.permutation_means.copy()
        if sd is None:
            sd = self.permutation_sd.copy()

        z_scores = (value - mean)/sd
        p_values = (1-stats.norm.cdf(z_scores))*2
        return p_values

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

    def initialize_params(self, n_trees=None):
        """
        Choose the value of alpha to use for fitting
        :param n_trees: float, optional
            The alpha value to use for the window. If none is entered the alpha will be chosen by cross validation
        :return:
        """
        if n_trees is None:
            # Select number of trees with default parameters
            self.n_trees = 100
        elif n_trees >= 0 and type(n_trees) == int:
            self.n_trees = n_trees
        else:
            raise ValueError("Number of trees must be int (>=0) or None")
        return

    def fit_window(self):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        # Make sure an alpha value has been selected to use for fitting the window
        if self.n_trees is None:
            raise ValueError("window alpha value must be set before the window can be fit")

        self.beta_coefficients = self.get_coeffs(self.n_trees)

    def get_coeffs(self, n_trees, data=None):
        """

        :param data:
        :param n_trees:
        :return: array-like
            An array in which the rows are childred and the columns are the parents
        """
        if data is None:
            all_data = self.window_values.copy()
        else:
            all_data = data.copy()

        coeff_matrix = np.array([], dtype=np.float_).reshape(0, all_data.shape[1])

        for col_index, column in enumerate(all_data.T):
            #print "Inferring parents for gene %i of %i" % (col_index, self.n_labels)
            #delete the column that is currently being tested
            X_matrix = np.delete(all_data, col_index, axis=1)

            #take out the column so that the gene does not regress on itself
            target_TF = all_data[:, col_index]

            rfr = RandomForestRegressor(n_estimators=n_trees, n_jobs=-1)
            rfr.fit(X_matrix, target_TF)
            coeffs = rfr.feature_importances_
            #artificially add a 0 to where the col_index is
            #to prevent self-edges
            coeffs = np.insert(coeffs, col_index, 0)
            coeff_matrix = np.vstack((coeff_matrix, coeffs))

        return coeff_matrix
