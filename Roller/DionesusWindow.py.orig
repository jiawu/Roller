__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from util.pls_nipals import vipp
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
import numpy as np
from Window import Window
import pandas as pd


class DionesusWindow(Window):
    """
    A window that runs Dionesus as the network inference algorithm. The PLSR function is from sci-kit learn for
    implementation consistency between window types

    For more information about Dionesus see:

    Ciaccio, Mark F., et al. "The DIONESUS algorithm provides scalable and accurate reconstruction of dynamic
    phosphoproteomic networks to reveal new drug targets." Integrative Biology (2015).

    """

    def __init__(self, dataframe, window_info, roller_data):
        super(DionesusWindow, self).__init__(dataframe, window_info, roller_data)
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

    def make_edge_table(self, calc_mse = False):
        """
        if self.permutation_p_values is None:
            raise ValueError("p values must be set before making the edge table. Use method run_permutation test")

        if self.edge_stability_auc is None:
            raise ValueError("edge stability values must be set before making the edge table. "
                             "Use method run_permutation test")
        self.edge_table["p_value"] = self.permutation_p_values.flatten()
        self.edge_table["stability"] = self.edge_stability_auc.flatten()
        """
        # For now, the edge table will only consist of VIP scores
        self.results_table['p_value'] = self.permutation_p_values.flatten()
        self.results_table['importance'] = np.asarray(self.vip).flatten()
        
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

    def run_permutation_test(self, n_permutations=1000, crag = False):

        # initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        nth_window = self.nth_window
        zeros = np.zeros((self.n_genes, self.n_genes))

        # initialize running calculation
        result = {'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}

        # inner loop: permute the window N number of times
        for nth_perm in range(0, n_permutations):
            # if (nth_perm % 200 == 0):
            # print 'Perm Run: ' +str(nth_perm)

            # permute data
            permuted_data = self.permute_data(self.window_values)

            # fit the data and get coefficients

            permuted_coeffs, permuted_vip = self.get_coeffs(data=permuted_data)
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

        z_scores = (value - mean) / sd
        cdf = stats.norm.cdf((-1 * abs(z_scores)))
        p_values = 2 * cdf
        return p_values

    def update_variance_2D(self, prev_result, new_samples):
        """incremental calculation of means: accepts new_samples, which is a list of samples. then calculates a new mean. this is a useful function for calculating the means of large arrays"""
        n = prev_result["n"]  # 2D numpy array with all zeros or watev
        mean = prev_result["mean"]  # 2D numpy array
        sum_squares = prev_result["ss"]  # 2D numpy array

        # new_samples is a list of arrays
        # x is a 2D array
        for x in new_samples:
            n = n + 1
            # delta = float(x) - mean
            old_mean = mean.copy()
            mean = old_mean + np.divide((x - old_mean), n)
            sum_squares = sum_squares + np.multiply((x - mean), (x - old_mean))

        if (n[0, 0] < 2):
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

        self.beta_coefficients, self.vip = self.get_coeffs(pcs, crag=crag)

    def sum_of_squares(self, X):
        """
        Calculate the sum of the squares for each column
        :param X: array-like
            The data matrix for which the sum of squares is taken
        :return: float or array-like
            The sum of squares, columnwise or total
        """
        column_mean = np.mean(X, axis=0)
        ss = np.sum(np.power(X - column_mean, 2), axis=0)
        return ss

    def _initialize_coeffs(self, data):
        """
        example call:
        all_data, coeff_matrix, model_list, max_nodes =self._initialize_coeffs(data=data)
        """
        if data is None:
            all_data = self.window_values.copy()
        else:
            all_data = data.copy()
        max_nodes = self.window_values.shape[1]

        coeff_matrix = np.array([], dtype=np.float_).reshape(0, all_data.shape[1])
        model_list = []
        return((all_data, coeff_matrix, model_list, max_nodes))
    
    def _fitstack_coeffs(self, num_pcs, coeff_matrix, vip_matrix, model_list, all_data, col_index, crag=False):
      """
      example call:
      coeff_matrix, model_list = self._fitstack_coeffs(coeff_matrix, model_list, all_data, col_index, n_trees, n_jobs,crag)
      """
      pls = PLSRegression(num_pcs, False)

      # delete the column that is currently being tested
      X_matrix = np.delete(all_data, col_index, axis=1)

      # take out the column so that the gene does not regress on itself
      target_TF = all_data[:, col_index]
      pls.fit(X_matrix, target_TF)
      model_params = {'col_index': col_index,
                      'response': target_TF,
                      'predictor': X_matrix,
                      'model': pls}

      model_list.append(model_params)

      # artificially add a 0 to where the col_index is to prevent self-edges
      coeffs = pls.coefs
      coeffs = np.insert(coeffs, col_index, 0)
      coeff_matrix = np.vstack((coeff_matrix, coeffs))

      # Calculate and store VIP scores
      vips = vipp(X_matrix, target_TF, pls.x_scores_, pls.x_weights_)
      vips = np.insert(vips, col_index, 0)
      vip_matrix = np.vstack((vip_matrix, vips))

      # scoping issues
      if crag:
          training_scores, test_scores = self.crag_window(model_params)
          self.training_scores.append(training_scores)
          self.test_scores.append(test_scores)
      return((coeff_matrix, vip_matrix, model_list))
    
    def get_coeffs(self, num_pcs=3, data=None, crag=False):
        """
        returns a 2D array with target as rows and regulators as columns
        :param alpha: float
            value to use for lasso fitting
        :param data: array-like, optional
            Data to fit. If none, will use the window values. Default is None
        :return:
        """
        all_data, coeff_matrix, model_list, max_nodes = self._initialize_coeffs(data=data)
        vip_matrix = coeff_matrix.copy()

        for col_index, column in enumerate(all_data.T):
            # Instantiate a new PLSR object
            coeff_matrix, vip_matrix, model_list = self._fitstack_coeffs(num_pcs = num_pcs, coeff_matrix=coeff_matrix, vip_matrix=vip_matrix, model_list=model_list, all_data=all_data, col_index=col_index, crag=crag)

        return coeff_matrix, vip_matrix
    
    def _permute_coeffs(self, zeros, crag=False, n_permutations=10):
        result={'n': zeros.copy(), 'mean': zeros.copy(), 'ss': zeros.copy()}
        for nth_perm in range(0, n_permutations):
            if self.x_data is not None:
                permuted_data = self.permute_data(self.x_data)
            else:
                permuted_data = self.permute_data(self.window_values)
            importance_df, permuted_coeffs = self.get_coeffs(data=permuted_data)
            dummy_list = []
            dummy_list.append(permuted_coeffs)
            result = self.update_variance_2D(result, dummy_list)
        self.permutation_means = result['mean'].copy()
        self.permutation_sd = np.sqrt(result['variance'].copy())
        self.permutation_p_values = self.calc_p_value()

class tdDionesusWindow(DionesusWindow):
    def __init__(self, dataframe, window_info, roller_data):
        super(tdDionesusWindow, self).__init__(dataframe, window_info, roller_data)
        self.x_data = None
        self.x_labels = None 
        self.x_times = None
        self.edge_table = None
        self.include_window = True
        self.earlier_window_idx = None
        
    def create_linked_list(self, numpy_array_2D, value_label):
        """labels and array should be in row-major order"""
        linked_list = pd.DataFrame({'regulator-target': self.edge_labels, value_label: numpy_array_2D.flatten()})
        return linked_list
    
    def resample_window(self):
        """
        Resample window values, along a specific axis
        :param window_values: array

        :return: array
        """
        n, p = self.x_data.shape

        # For each column randomly choose samples
        resample_values = np.array([np.random.choice(self.x_data[:, ii], size=n) for ii in range(p)]).T

        # resample_window = pd.DataFrame(resample_values, columns=self.df.columns.values.copy(),
        #                               index=self.df.index.values.copy())
        return resample_values 

    def fit_window(self, pcs=3, crag=False, calc_mse=False):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        self.vip, self.beta_coefficients = self.get_coeffs(pcs, crag=crag, data=self.x_data)
        self.edge_importance = self.vip.copy()
        
    def get_coeffs(self, num_pcs=3,crag=False, data=None):
        """
        :param data:
        :param n_trees:
        :return: array-like
            An array in which the rows are children and the columns are the parents
        """
        if data is None:
            data = self.x_data
        ## initialize items
        all_data, coeff_matrix, model_list, max_nodes = self._initialize_coeffs(data=data)
        vip_matrix = coeff_matrix.copy()
      

        for col_index, column in enumerate(all_data[:,:max_nodes].T):
            # Once we get through all the nodes at this timepoint we can stop
            if col_index == max_nodes:
                break
            coeff_matrix, vip_matrix, model_list = self._fitstack_coeffs(num_pcs=num_pcs,coeff_matrix=coeff_matrix, vip_matrix=vip_matrix, model_list=model_list, all_data=all_data, col_index=col_index,crag=False) 
            
        """
        if self.x_labels == None:
            label = self.raw_data.columns[1:]
            importance_dataframe = pd.DataFrame(coeff_matrix, index=label[:max_nodes], columns=label)
        else:
        """
        importance_dataframe = pd.DataFrame(vip_matrix, index=self.x_labels[:max_nodes], columns=self.x_labels)
        importance_dataframe.index.name = 'Child'
        importance_dataframe.columns.name = 'Parent'
        return(importance_dataframe,coeff_matrix)

    def make_edge_table(self, calc_mse = False):
        """
        Make the edge table
        :return:
        """

        if not self.include_window:
            return

        # Build indexing method for all possible edges. Length = number of parents * number of children
        parent_index = range(self.edge_importance.shape[1])
        child_index = range(self.edge_importance.shape[0])
        a, b = np.meshgrid(parent_index, child_index)

        # Flatten arrays to be used in link list creation
        df = pd.DataFrame()
        df['Parent'] = self.edge_importance.columns.values[a.flatten()]
        df['Child'] = self.edge_importance.index.values[b.flatten()]
        df['Importance'] = self.edge_importance.values.flatten()
        df['P_window'] = self.x_times[a.flatten()]
        df['C_window'] = self.x_times[b.flatten()]
        if self.permutation_p_values is not None:
            df["p_value"] = self.permutation_p_values.flatten()

        return df

    def run_permutation_test(self, n_permutations=10, crag=False):
        #if not self.include_window:
        #    return
        #initialize permutation results array
        self.permutation_means = np.empty(self.edge_importance.shape)
        self.permutation_sd = np.empty(self.edge_importance.shape)

        zeros = np.zeros(self.edge_importance.shape)
        self._permute_coeffs(zeros, crag=False, n_permutations=n_permutations)
