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
        self.edge_importance = None
        self.n_trees = None
        self.bootstrap_matrix = None
        self.freq_matrix = None
        self.permutation_means = None
        self.permutation_sd = None
        self.permutation_p_values = None
        self.permutation_pvalues = None

    def make_edge_table(self):
        if self.permutation_p_values is None:
            raise ValueError("p values must be set before making the edge table. Use method run_permutation test")

        if self.edge_importance is None:
            raise ValueError("edge importance values must be set before making the edge table. "
                             "Use method run_permutation test")
        self.edge_table["P_Value"] = self.permutation_p_values.flatten()
        self.edge_table["Importance"] = self.edge_importance.flatten()


    def rank_edges(self, method="p_value"):
        if self.edge_table is None:
            raise ValueError("The edge table must be created before getting edges")
        if method == "p_value":
            self.edge_table.sort(columns=['P_Value', 'Importance'], ascending=[True, False], inplace=True)
        elif method == "importance":
            self.edge_table.sort(columns=['Stability', 'Importance'], ascending=[False, True], inplace=True)
        return

    def run_permutation_test(self, n_permutations=1000, n_jobs=-1):
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

            permuted_coeffs = self.get_coeffs(self.n_trees, permuted_data, n_jobs=n_jobs)
            dummy_list = []
            dummy_list.append(permuted_coeffs)
            result = self.update_variance_2D(result, dummy_list)

        self.permutation_means = result['mean'].copy()
        self.permutation_sd = np.sqrt(result['variance'].copy())
        self.permutation_p_values = self.calc_p_value()

    def calc_p_value(self, value=None, mean=None, sd=None):
        if value is None:
            value = self.edge_importance.copy()
        if mean is None:
            mean = self.permutation_means.copy()
        if sd is None:
            sd = self.permutation_sd.copy()

        z_scores = (value - mean)/sd
        p_values = (1-stats.norm.cdf(z_scores))*2
        return p_values

    def get_nth_window_auc(self, nth):
        auc = self.edge_importance[:,:, nth]
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

    def fit_window(self, n_jobs=-1):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        # Make sure an alpha value has been selected to use for fitting the window
        if self.n_trees is None:
            raise ValueError("window alpha value must be set before the window can be fit")

        self.edge_importance = self.get_coeffs(self.n_trees, n_jobs=n_jobs)

    def get_coeffs(self, n_trees, data=None, n_jobs=-1):
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

            rfr = RandomForestRegressor(n_estimators=n_trees, n_jobs=n_jobs)
            rfr.fit(X_matrix, target_TF)
            coeffs = rfr.feature_importances_
            #artificially add a 0 to where the col_index is
            #to prevent self-edges
            coeffs = np.insert(coeffs, col_index, 0)
            coeff_matrix = np.vstack((coeff_matrix, coeffs))

        return coeff_matrix
