from __future__ import absolute_import
import Roller
from sklearn.preprocessing import Imputer
from Roller.util.linear_wrapper import LassoWrapper
import numpy as np
import matplotlib as mpl
import pdb
import Roller.util.Grapher as gr
from Roller.util.permutation_test import Permuter

file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
#file_path = "/Users/jjw036/Roller/goldbetter_model/goldbetter_data.txt"
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None

roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
window_size = 5
#get only TFs data, window size of 4
roll_me.set_window(window_size)
#impute missing values
imputer = Imputer(missing_values="NaN")
mpl.rcParams.update({'font.size':8})
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
roll_me.remove_blank_rows()
total_window_number = roll_me.get_n_windows()
#todo: fix the below two lines, I don't need this to get the matrix size...
current_window = roll_me.get_window()
filled_matrix = imputer.fit_transform(current_window)
n_genes = filled_matrix.shape[1]
coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))
gene_names=list(current_window.columns.values)

for nth_window in range(0, total_window_number):
    #loop gets the window, gets the coefficients for that window, then increments the window
    current_window = roll_me.get_window()

    #check if any values are completely blank
    current_window = current_window
    filled_matrix = current_window.values
    #filled_matrix = imputer.fit_transform(current_window)
    current_lasso = LassoWrapper(filled_matrix)
    coeff_mat = current_lasso.get_coeffs(0.002)
    coeff_matrix_3d[:,:,nth_window] = coeff_mat
    #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
    roll_me.next()

permuter = Permuter()

#give it a roller object
permuter.run_permutation_test(roll_me, alpha = 0.002)
pdb.set_trace()

# hmm maybe make a heatmap for each slice...
# get the binary coefficients of each index
binary_matrix = np.empty((n_genes,n_genes))
for row_index in range(n_genes):
    for col_index in range(n_genes):
        coeffs_over_time = coeff_matrix_3d[row_index,col_index,:]
        if sum(coeffs_over_time) == 0:
            binary_matrix[row_index,col_index]=0
        else:
            binary_matrix[row_index,col_index]=1

ii,jj = np.where(binary_matrix==1)
present = np.append(ii,jj)
present = np.unique(present)
present_ii = np.unique(ii)
present_jj = np.unique(jj)

present_mat = coeff_matrix_3d[present_ii,:,:]
present_mat = present_mat[:,present_jj,:]
present_genes_ii = [gene_names[item] for item in present_ii]
present_genes_jj = [gene_names[item] for item in present_jj]
for nth_window in range(0, total_window_number):
    coeff_mat = present_mat[:,:,nth_window]
    #gr.plot_figure(coeff_mat,nth_window,present_genes_ii, present_genes_jj, window_size,str("compressed_window"+str(window_size)))

    #move all the non-zero coefficients into a binary matrix
series_list = []

#targets of plot
for gene in present_genes_jj:
    col_index = present_genes_jj.index(gene)
    present_coeffs_over_time = present_mat[:,col_index,:]
    binary_matrix_tf = np.empty((present_coeffs_over_time.shape[0]))
    tf_list = []
    for row_index in range(len(present_genes_ii)):
        coeffs_TF = present_coeffs_over_time[row_index,:]
        if sum(coeffs_TF) == 0:
            binary_matrix_tf[row_index]=0
        else:
            binary_matrix_tf[row_index]=1
            tf_list.append(coeffs_TF)

    present_targets = np.where(binary_matrix_tf==1)[0]
    tfs_of_interest = [ (col_index,y) for y in present_targets]
    regulator_labels = [gene]*len(tfs_of_interest)
    target_labels = [present_genes_ii[item] for item in present_targets]
    #pdb.set_trace()
    #np.savetxt('binary_mat')
    #plot_lines(tf_list, regulator_labels, target_labels, window_size, suffix=str("targets_of_"+gene))

