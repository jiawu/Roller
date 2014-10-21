from sklearn import linear_model
import numpy as np

class LassoWrapper:
    #I want this to spit out a matrix of coefficients, where row is the gene target and columns are the regulators

    #potentially I want to do cross validation too to see if the cross validation score is better rolling vs non-rolling. I need to think of the strategy for this.

    def __init__(self, data_frame):
        self.data = data_frame

    def get_coeffs(self, alpha=0.2):
        """returns a 2D array with target as rows and regulators as columns"""
        clf = linear_model.Lasso(alpha)
        #loop that iterates through the target genes
        all_data = self.data
        coeff_matrix = np.array([],dtype=np.float_).reshape(0,all_data.shape[1])

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
        return(coeff_matrix)
