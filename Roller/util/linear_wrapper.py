from sklearn import linear_model
import numpy as np

#Note: if the optimal alpha value is very small, then another method other than LASSO should be used.

class LassoWrapper:
    #I want this to spit out a matrix of coefficients, where row is the gene target and columns are the regulators


    def __init__(self, data_frame):
        #todo: data_frame is really a numpy array. Make a better name
        self.data = data_frame

    def get_coeffs(self, alpha=0.2):
        """returns a 2D array with target as rows and regulators as columns"""
        clf = linear_model.Lasso(alpha)
        #loop that iterates through the target genes
        all_data = self.data
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

    def get_max_alpha(self):

        # Get maximum edges, assuming all explanors are also response variables and no self edges
        [n, p] = self.data.shape
        max_edges = p * (p-1)

        # Raise exception if Lasso doesn't converge with alpha == 0
        if np.count_nonzero(self.get_coeffs(0)) != max_edges:
            raise ValueError('Lasso does not converge with alpha = 0')

        # Start stepping with forward like selection
        max_step_size = 1e4
        min_step_size = 1e-9
        powers = int(np.log10(max_step_size/min_step_size))
        step_sizes = [max_step_size/(10**ii) for ii in range(powers+1)]
        cur_min = 0
        alpha_max = step_sizes[0]
        for ii, cur_max in enumerate(step_sizes[:-1]):
            if alpha_max > cur_max:
                cur_max = alpha_max
            cur_step = step_sizes[ii+1]
            cur_range = np.linspace(cur_min, cur_max, (cur_max-cur_min)/cur_step+1)
            for cur_alpha in cur_range:
                num_coef = np.count_nonzero(self.get_coeffs(cur_alpha))
                if num_coef > 0:
                    cur_min = cur_alpha
                elif num_coef == 0:
                    alpha_max = cur_alpha
                    break
        return alpha_max





