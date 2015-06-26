__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from util.pls_nipals import vipp
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
import numpy as np
from Window import Window

class DionesusWindow(Window):
    """
    A window that runs Dionesus as the network inference algorithm

    For more information about Dionesus see:

    Ciaccio, Mark F., et al. "The DIONESUS algorithm provides scalable and accurate reconstruction of dynamic
    phosphoproteomic networks to reveal new drug targets." Integrative Biology (2015).

    """

    def __init__(self, dataframe, window_info, roller_data):
        super(DionesusWindow, self).__init__(dataframe, window_info, roller_data)
        self.alpha = None
        self.beta_coefficients = None
        self.cv_table = None
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
        self.edge_table["p_value"] = self.permutation_p_values.flatten()
        self.edge_table["stability"] = self.edge_stability_auc.flatten()

    def rank_edges(self, method="p_value"):
        if self.edge_table is None:
            raise ValueError("The edge table must be created before getting edges")
        temp_edge_table = self.edge_table.copy()
        if method == "p_value":
            temp_edge_table.sort(columns=['p_value', 'stability'], ascending=[True, False], inplace=True)
        elif method == "stability":
            temp_edge_table.sort(columns=['stability', 'p_value'], ascending=[False, True], inplace=True)

        return temp_edge_table['regulator-target'].values

    def generate_results_table(self):

        #generate edges for initial model
        initial_edges = self.create_linked_list(self.beta_coefficients, 'B')
        #permutation edges
        permutation_mean_edges =self.create_linked_list(self.permutation_means, 'p-means')
        permutation_sd_edges = self.create_linked_list(self.permutation_sd, 'p-sd')
        stability_edges = self.create_linked_list(self.edge_stability_auc, 'stability')

        aggregated_edges = initial_edges.merge(permutation_mean_edges, on='regulator-target').merge(permutation_sd_edges, on='regulator-target').merge(stability_edges, on='regulator-target')

        #sorry, it is a little messy to do the p-value calculations for permutation tests here...
        #valid_indices = aggregated_edges['p-sd'] != 0
        #valid_indices = aggregated_edges['B'] != 0
        valid_window = aggregated_edges
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

            permuted_coeffs = self.get_coeffs(self.alpha, permuted_data)
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
        cdf = stats.norm.cdf((-1*abs(z_scores)))
        p_values = 2*cdf
        return p_values

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

    def get_coeffs(self, num_pcs=3, data=None):
        """
        returns a 2D array with target as rows and regulators as columns
        :param alpha: float
            value to use for lasso fitting
        :param data: array-like, optional
            Data to fit. If none, will use the window values. Default is None
        :return:
        """

        #loop that iterates through the target genes
        if data is None:
            all_data = self.window_values.copy()
        else:
            all_data = data.copy()

        coeff_matrix = np.array([],dtype=np.float_).reshape(0, all_data.shape[1])

        model_list = []

        for col_index,  column in enumerate(all_data.T):
            # Instantiate a new PLSR object
            pls = PLSRegression(num_pcs, False)

            #delete the column that is currently being tested
            X_matrix = np.delete(all_data, col_index, axis=1)

            #take out the column so that the gene does not regress on itself
            target_TF = all_data[:,col_index]
            pls.fit(X_matrix, target_TF)
            model_params = {'col_index':col_index,
                            'response':target_TF,
                            'predictor':X_matrix,
                            'model':pls}

            model_list.append(model_params)
            coeffs = pls.coef_
            #artificially add a 0 to where the col_index is
            #to prevent self-edges
            coeffs = np.insert(coeffs, col_index, 0)
            coeff_matrix=np.vstack((coeff_matrix, coeffs))

            #scoping issues
            training_scores, test_scores = self.crag_window(model_params)
            self.training_scores.append(training_scores)
            self.test_scores.append(test_scores)

        return coeff_matrix
