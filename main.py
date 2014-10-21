from Roller import Roller
from sklearn.preprocessing import Imputer
from linear_wrapper import LassoWrapper
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

import pdb

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot_figure(coeff_matrix,nth_window):
    df = pd.DataFrame(coeff_matrix)
    figure1 = plt.figure()
    heatmap = figure1.add_subplot(1,1,1)
    my_axis = heatmap.imshow(df,interpolation='nearest',cmap=cm.OrRd)
    my_axi = my_axis.get_axes()
    clean_axis(my_axi)
    figure1.savefig('figure'+str(nth_window)+'.png')

file_path = "compressed_katrina_data.txt"
gene_start_column = 5
roll_me = Roller(file_path, gene_start_column)

#get only TFs data, window size of 4
roll_me.set_window = 4
#impute missing values
imputer = Imputer(missing_values="NaN")

total_window_number = roll_me.get_n_windows()
#todo: fix the below two lines, I don't need this to get the matrix size...
current_window = roll_me.get_window()
filled_matrix = imputer.fit_transform(current_window)
n_genes = filled_matrix.shape[1]
coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))

for nth_window in range(0, total_window_number):
    #loop gets the window, gets the coefficients for that window, then increments the window
    current_window = roll_me.get_window()

    #check if any values are completely blank
    sums = [current_window.iloc[:,x].sum() for x in range(0,56)]
    ind = np.where(np.isnan(sums))[0]
    current_window.iloc[:,ind]=0

    filled_matrix = imputer.fit_transform(current_window)
    current_lasso = LassoWrapper(filled_matrix)
    coeff_mat = current_lasso.get_coeffs()
    pdb.set_trace()
    coeff_matrix_3d[:,:,nth_window] = coeff_mat
    plot_figure(coeff_mat,nth_window)
    roll_me.next()

gene_names=list(current_window.columns.values)
# hmm maybe make a heatmap for each slice...



