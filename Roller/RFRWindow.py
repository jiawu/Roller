import warnings
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import KFold
from scipy import integrate
from scipy import stats
import scipy
import time

from Window import Window

class RandomForestRegressionWindow(Window):
    def __init__(self, dataframe, window_info, roller_data):
        super(RandomForestRegressionWindow, self).__init__(dataframe, window_info, roller_data)
        self.edge_importance = None
        self.n_trees = None
        self.n_jobs = None
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
        self.results_table["p_value"] = self.permutation_p_values.flatten()
        self.results_table["importance"] = self.edge_importance.flatten()


    def sort_edges(self, rank_by):
        rank_column_name = rank_by + "-rank"
        if self.results_table is None:
            raise ValueError("The edge table must be created before getting edges")
        if rank_by == "p_value":
            self.results_table.sort(columns=['p_value', 'importance'], ascending=[True, False], inplace=True)
        elif rank_by == "importance":
            self.results_table.sort(columns=['importance','p_value'], ascending=[False, True], inplace=True)
        valid_window = self.results_table
        valid_window[rank_column_name] = valid_window[rank_by].rank(method="dense",ascending=False)
        self.results_table = self.results_table.sort(columns=rank_column_name, axis=0)
        return self.results_table

    def _permute_coeffs(self, zeros, crag,n_permutations, n_jobs, td=False):
        """
        Internal method that is shared between window and tdWindow

        """
        #initialize running calculation
        result = {'n':zeros.copy(), 'mean':zeros.copy(), 'ss':zeros.copy()}
        #inner loop: permute the window N number of times
        for nth_perm in range(0, n_permutations):
            #if (nth_perm % 200 == 0):
                #print 'Perm Run: ' +str(nth_perm)

            #permute data
            if td == True:
                permuted_data = self.permute_data(self.x_data)
            else:
                permuted_data = self.permute_data(self.window_values)

            #fit the data and get coefficients

            permuted_coeffs = self.get_coeffs(self.n_trees, crag=crag, data=permuted_data, n_jobs=n_jobs)
            dummy_list = []
            dummy_list.append(permuted_coeffs)
            result = self.update_variance_2D(result, dummy_list)

        self.permutation_means = result['mean'].copy()
        self.permutation_sd = np.sqrt(result['variance'].copy())
        self.permutation_p_values = self.calc_p_value()
    
    def run_permutation_test(self, crag=False, n_permutations=1000, n_jobs=-1):
        #initialize permutation results array
        self.permutation_means = np.empty((self.n_genes, self.n_genes))
        self.permutation_sd = np.empty((self.n_genes, self.n_genes))
        zeros = np.zeros((self.n_genes, self.n_genes))
        self._permute_coeffs(zeros, crag=crag, n_permutations=n_permutations, n_jobs=n_jobs)

    def calc_p_value(self, value=None, mean=None, sd=None):
        if value is None:
            value = self.edge_importance.copy()
        if mean is None:
            mean = self.permutation_means.copy()
        if sd is None:
            sd = self.permutation_sd.copy()

        z_scores = (value - mean)/sd
        cdf = stats.norm.cdf((-1*abs(z_scores)))
        p_values = 2*cdf
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
            self.n_trees = 500
        elif n_trees >= 0 and type(n_trees) == int:
            self.n_trees = n_trees
        else:
            raise ValueError("Number of trees must be int (>=0) or None")
        return

    def fit_window(self, crag=True, n_jobs=-1):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        # Make sure an alpha value has been selected to use for fitting the window
        if self.n_trees is None:
            raise ValueError("window alpha value must be set before the window can be fit")

        self.edge_importance = self.get_coeffs(self.n_trees, crag=crag, n_jobs=n_jobs)

    def _initialize_coeffs(self, data):
        """ Returns a copy of the vector, an empty array with a defined shape, an empty list, and the maximum number of nodes
        """
        if data is None:
            all_data = self.window_values.copy()
        else:
            all_data = data.copy()

        max_nodes = self.window_values.shape[1]

        coeff_matrix = np.array([], dtype=np.float_).reshape(0, all_data.shape[1])

        model_list = []
        return((all_data,coeff_matrix,model_list,max_nodes))
  
    def _fitstack_coeffs(self, coeff_matrix, model_list, x_matrix, target_y, col_index, n_trees, n_jobs, crag):
        #print "Inferring parents for gene %i of %i" % (col_index, self.n_labels)
        #X_matrix = np.delete(all_data, col_index, axis=1)

        #take out the column so that the gene does not regress on itself
        #target_TF = all_data[:, col_index]

        rfr = RandomForestRegressor(n_estimators=n_trees, n_jobs=n_jobs, max_features="sqrt")
        rfr.fit(x_matrix, target_y)
        model_params = {'col_index':col_index,
                        'response':target_y,
                        'predictor':x_matrix,
                        'model':rfr}
        model_list.append(model_params)
        coeffs = rfr.feature_importances_

        #artificially add a 0 to where the col_index is
        #to prevent self-edges
        if coeff_matrix.shape[1]-len(coeffs) == 1:
            coeffs = np.insert(coeffs, col_index, 0)
        coeff_matrix = np.vstack((coeff_matrix, coeffs))
        # there's some scoping issues here. cragging needs the roller's raw data but the window does not know what roller contains (outside scope). have to pass in the roller's raw data and save it somehow :/
        if crag == True:
            training_scores, test_scores = self.crag_window(model_params)
            self.training_scores.append(training_scores)
            self.test_scores.append(test_scores)
        
        return (coeff_matrix, model_list)
        
    def get_coeffs(self, n_trees, crag=True, data=None, n_jobs=-1):
        """

        :param data:
        :param n_trees:
        :return: array-like
            An array in which the rows are childred and the columns are the parents
        """
        ## initialize items
        all_data, coeff_matrix, model_list, max_nodes =self._initialize_coeffs(data=data)

        for col_index, column in enumerate(all_data.T):
            coeff_matrix, model_list = self._fitstack_coeffs(coeff_matrix, model_list, all_data, col_index, n_trees, n_jobs,crag) 
        return coeff_matrix

class tdRFRWindow(RandomForestRegressionWindow):
    def __init__(self, dataframe, window_info, roller_data):
        super(tdRFRWindow, self).__init__(dataframe, window_info, roller_data)
        self.x_data = None
        self.x_labels = None
        self.x_times = None
        self.edge_table = None
        self.include_window = True
        self.earlier_window_idx = None

    def fit_window(self, crag=False ,n_jobs=1):
        """
        Set the attributes of the window using expected pipeline procedure and calculate beta values
        :return:
        """
        if self.include_window:
            print "Regressing window index %i against the following window indices: "%self.nth_window,\
                self.earlier_window_idx
            self.edge_importance = self.get_coeffs(self.n_trees, crag=False, data=self.x_data, n_jobs=n_jobs)

    def get_coeffs(self, n_trees, crag=False, data=None, n_jobs=-1):
        """

        :param data:
        :param n_trees:
        :return: array-like
            An array in which the rows are children and the columns are the parents
        """

        ## initialize items
        all_data, coeff_matrix, model_list, max_nodes = self._initialize_coeffs(data=data)

        for col_index, column in enumerate(self.window_values.T):
            target_y = column.copy()
            x_matrix = all_data.copy()
            if self.nth_window in self.x_times:
                x_matrix = np.delete(x_matrix, col_index, axis=1)
            coeff_matrix, model_list = self._fitstack_coeffs(coeff_matrix, model_list, x_matrix, target_y, col_index, n_trees, n_jobs,crag)
        importance_dataframe = pd.DataFrame(coeff_matrix, index=self.x_labels[:max_nodes], columns=self.x_labels)
        importance_dataframe.index.name = 'Child'
        importance_dataframe.columns.name = 'Parent'
        return importance_dataframe

    def make_edge_table(self):
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

    def run_permutation_test(self, n_permutations=1000, n_jobs=1,crag=False):
        if not self.include_window:
            return
        #initialize permutation results array
        self.permutation_means = np.empty(self.edge_importance.shape)
        self.permutation_sd = np.empty(self.edge_importance.shape)

        zeros = np.zeros(self.edge_importance.shape)
        self._permute_coeffs(zeros, crag=False, n_permutations=n_permutations, n_jobs=n_jobs, td=True)

        
