__author__ = 'jfinkle'


from Swing import Swing
from sklearn.preprocessing import Imputer
from linear_wrapper import LassoWrapper
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import matplotlib as mpl
import sys

file_path = "goldbetter_data.txt"
gene_start_column = 2
roll_me = Swing(file_path, gene_start_column)
window_size = 10000-1
#get only TFs data, window size of 4
roll_me.set_window(window_size)
#impute missing values
imputer = Imputer(missing_values="NaN")
mpl.rcParams.update({'font.size':8})
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

total_window_number = roll_me.get_n_windows()
#todo: fix the below two lines, I don't need this to get the matrix size...
current_window = roll_me.get_window()
filled_matrix = imputer.fit_transform(current_window)
n_genes = filled_matrix.shape[1]
coeff_matrix_3d = np.empty((n_genes, n_genes, total_window_number))
gene_names = list(current_window.columns.values)

current_lasso = LassoWrapper(filled_matrix)
coeff_mat = current_lasso.get_coeffs(alpha = 0.005)

print np.count_nonzero(coeff_mat)/(16.0**2)
