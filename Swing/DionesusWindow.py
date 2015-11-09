__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from util.pls_nipals import vipp
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
import numpy as np
from Window import Window
from sklearn.metrics import mean_squared_error
import pandas as pd
import sys


class DionesusWindow(Window):
    """
    A window that runs Dionesus as the network inference algorithm. The PLSR function is from sci-kit learn for
    implementation consistency between window types

    For more information about Dionesus see:

    Ciaccio, Mark F., et al. "The DIONESUS algorithm provides scalable and accurate reconstruction of dynamic
    phosphoproteomic networks to reveal new drug targets." Integrative Biology (2015).

    """

    def __init__(self, dataframe, window_info, roller_data, td_window, explanatory_dict, response_dict):
        super(DionesusWindow, self).__init__(dataframe, window_info, roller_data,  td_window, explanatory_dict,
                                             response_dict)
        self.beta_coefficients = None
        self.vip = None
        self.cv_table = None
        self.bootstrap_matrix = None
        self.freq_matrix = None
        self.edge_stability_auc = None
        self.permutation_means = None
        self.permutation_sd = None
        self.permutation_p_values = None
        self.permutation_pvalues = None

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
        df['Importance'] = self.vip.values.flatten()
        df['Beta'] = self.beta_coefficients.values.flatten()
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
        
    def sort_edges(self, method="importance"):
        if self.results_table is None:
            raise ValueError("The edge table must be created before getting edges")
        if method == "p_value":
            self.results_table.sort(columns=['p_value', 'importance'], ascending=[True, False], inplace=True)
        elif method == "importance":
            self.results_table.sort(columns=['importance', 'p_value'], ascending=[False, True], inplace=True)

        return self.results_table['regulator-target'].values

    def generate_results_table(self):

        # generate edges for initial model
        initial_edges = self.create_linked_list(self.beta_coefficients, 'B')

        # permutation edges
        permutation_mean_edges = self.create_linked_list(self.permutation_means, 'p-means')
        permutation_sd_edges = self.create_linked_list(self.permutation_sd, 'p-sd')
        stability_edges = self.create_linked_list(self.edge_stability_auc, 'stability')

        aggregated_edges = initial_edges.merge(permutation_mean_edges, on='regulator-target').merge(
            permutation_sd_edges, on='regulator-target').merge(stability_edges, on='regulator-target')

        # sorry, it is a little messy to do the p-value calculations for permutation tests here...
        # valid_indices = aggregated_edges['p-sd'] != 0
        # valid_indices = aggregated_edges['B'] != 0
        valid_window = aggregated_edges
        initial_B = valid_window['B']
        sd = valid_window['p-sd']
        mean = valid_window['p-means']
        valid_window['final-z-scores-perm'] = (initial_B - mean) / sd
        valid_window['cdf-perm'] = (-1 * abs(valid_window['final-z-scores-perm'])).apply(stats.norm.cdf)
        # calculate t-tailed pvalue
        valid_window['p-value-perm'] = (2 * valid_window['cdf-perm'])
        self.results_table = valid_window
        return (self.results_table)

    def rank_results(self, rank_by, ascending=False):
        rank_column_name = rank_by + "-rank"
        # rank edges with an actual beta value first until further notice ##
        valid_indices = self.results_table['B'] != 0
        valid_window = self.results_table[valid_indices]
        valid_window[rank_column_name] = valid_window[rank_by].rank(method="dense", ascending=ascending)
        edge_n = len(valid_window.index)

        invalid_indices = self.results_table['B'] == 0
        invalid_window = self.results_table[invalid_indices]
        invalid_window[rank_column_name] = invalid_window[rank_by].rank(method="dense", ascending=ascending)
        invalid_window[rank_column_name] += edge_n
        self.results_table = valid_window.append(invalid_window)
        self.results_table = self.results_table.sort(columns=rank_column_name, axis=0)

        return (self.results_table)

    def run_permutation_test(self, n_permutations=1000, crag=False):

        # initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        zeros = np.zeros(self.beta_coefficients.shape)

        # initialize running calculation
        result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        # inner loop: permute the window N number of times
        for nth_perm in range(0, n_permutations):
            # if (nth_perm % 200 == 0):
            # print 'Perm Run: ' +str(nth_perm)

            # permute data
            permuted_data = self.permute_data(self.explanatory_data)

            # fit the data and get coefficients

            permuted_coeffs, permuted_vip, _ = self.get_coeffs(x_data=permuted_data)
            dummy_list = [permuted_coeffs]
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

        z_scores = (value - mean) / sd
        cdf = stats.norm.cdf((-1 * abs(z_scores)))
        p_values = 2 * cdf
        return p_values

    def initialize_params(self):
        """
        Nothing to initialize for Dionesus
        :return:
        """
        pass

    def fit_window(self, pcs=3, crag=False, calc_mse=False):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values

        :return:
        """

        self.beta_coefficients, self.vip, self.edge_mse_diff = self.get_coeffs(pcs, crag=crag, calc_mse=calc_mse)
    
    def _fitstack_coeffs(self, n_pcs, coeff_matrix, vip_matrix, model_list, x_matrix, target_y, col_index, crag=False):
        """
        :param n_pcs:
        :param coeff_matrix:
        :param vip_matrix:
        :param model_list:
        :param x_matrix:
        :param target_y:
        :param col_index:
        :param crag:
        :return:
        """
        pls = PLSRegression(n_pcs, False)

        # Fit the model
        pls.fit(x_matrix, target_y)

        model_params = {'col_index': col_index,
                          'response': target_y,
                          'predictor': x_matrix,
                          'model': pls}

        model_list.append(model_params)

        # artificially add a 0 to where the col_index is to prevent self-edges
        coeffs = pls.coefs
        coeffs = np.reshape(coeffs, (len(coeffs),))
        vips = vipp(x_matrix, target_y, pls.x_scores_, pls.x_weights_)
        vips = np.reshape(vips, (len(vips),))
        if coeff_matrix.shape[1] - len(coeffs) == 1:
            coeffs = np.insert(coeffs, col_index, 0)
            vips = np.insert(vips, col_index, 0)

        coeff_matrix = np.vstack((coeff_matrix, coeffs))
        vip_matrix = np.vstack((vip_matrix, vips))

        # scoping issues
        if crag:
            training_scores, test_scores = self.crag_window(model_params)
            self.training_scores.append(training_scores)
            self.test_scores.append(test_scores)
        return coeff_matrix, vip_matrix, model_list
    
    def get_coeffs(self, num_pcs=3, x_data=None, y_data=None, crag=False, calc_mse=False):
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
        vip_matrix = coeff_matrix.copy()
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

            coeff_matrix, vip_matrix, model_list = self._fitstack_coeffs(num_pcs, coeff_matrix, vip_matrix, model_list,
                                                                         x_matrix, target_y, insert_index, crag=crag)

            base_mse = mean_squared_error(model_list[col_index]['model'].predict(x_matrix), target_y)

            if calc_mse:
                f_coeff_matrix, f_model_list = self._initialize_coeffs(data=x_matrix)
                f_vip_matrix = f_coeff_matrix.copy()
                mse_list = []
                for idx in range(x_matrix.shape[1]):
                    adj_x_matrix = np.delete(x_matrix, idx, axis=1)
                    f_coeff_matrix, f_vip_matrix, f_model_list = self._fitstack_coeffs(num_pcs, f_coeff_matrix,
                                                                                       f_vip_matrix, f_model_list,
                                                                                       adj_x_matrix, target_y,
                                                                                       idx, crag)
                    mse_diff = base_mse - mean_squared_error(f_model_list[idx]['model'].predict(adj_x_matrix), target_y)
                    mse_list.append(mse_diff)
                if mse_matrix is None:
                    mse_matrix = np.array(mse_list)
                else:
                    mse_matrix = np.vstack((mse_matrix, np.array(mse_list)))

        coeff_dataframe = pd.DataFrame(coeff_matrix, index=y_labels, columns=x_labels)
        coeff_dataframe.index.name = 'Child'
        coeff_dataframe.columns.name = 'Parent'

        importance_dataframe = pd.DataFrame(vip_matrix, index=y_labels, columns=x_labels)
        importance_dataframe.index.name = 'Child'
        importance_dataframe.columns.name = 'Parent'
        return coeff_dataframe, importance_dataframe, mse_matrix