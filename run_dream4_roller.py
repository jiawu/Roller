from __future__ import absolute_import
import Roller
from sklearn.preprocessing import Imputer
from Roller.util.linear_wrapper import LassoWrapper
import numpy as np
import matplotlib as mpl
from Roller.util.permutation_test import Permuter
import itertools
from Roller.util import utility_module as utility
from Roller.util.Ranker import Bootstrapper
import pdb

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

#zscore all data

roll_me.zscore_all_data()

#regulator-target labels
edge_labels = [x for x in itertools.product(gene_names,repeat=2)]

##initialize pandas panel with edge stats
##todo: delete this part
column_names = ['regulator-target', 'B', 'p-test-B-mean', 'p-test-B-std','p-test-p-values', 'stability']


##build initial model##
print("Creating initial model...")
for nth_window in range(0,total_window_number):
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

##convert coeff model to edge list
##takes a 3D matrix and gene names list, converts it to a 3D edge list. it assumes that the 3d matrix is symmetric and square.

#initial model conversion
initial_model = utility.create_3D_linked_list(edge_labels, coeff_matrix_3d, 'B')

roll_me.reset()
print("Running permutation test...")
#start permutation test
permuter = Permuter()
#give it a roller object
permuter.run_permutation_test(roll_me, alpha = 0.002)
perm_means=permuter.permutation_means
perm_sd=permuter.permutation_sd

permuted_model_means = utility.create_3D_linked_list(edge_labels, perm_means, 'p-means')
permuted_model_sd = utility.create_3D_linked_list(edge_labels, perm_sd, 'p-sd')

#get the permutation test beta matrix, mean and standard deviation.

roll_me.reset()
print("Running bootstrapping test...")
booter = Bootstrapper(roll_me)
boots = 100
max_random = 0.1
n_alphas = 100

booted_alphas = booter.run_bootstrap(window_size, boots, n_alphas, noise=max_random)
sums = np.sum(booter.freq_matrix,axis=3)
auc = []
for nth_window in range(0,total_window_number):
    auc.append(booter.get_nth_window_auc(0))
auc = np.dstack(auc)
#get 3d coeff matrix for stability scores
stability_model = utility.create_3D_linked_list(edge_labels, auc, 'stability')


#merge panels into one large panel with B, p-means, p-sd, stability, and p-value
all_panels = [initial_model,permuted_model_means, permuted_model_sd, stability_model]
for nth_window in range(0,total_window_number):
    #merge panels
    for panel_index in range(1,len(all_panels)):
        all_panels[0].merge(all_panels, panel_index, on='regulator-target')

pdb.set_trace()
#Research Question: Does rolling regression improve the sensitivity of inference methods? Compare lasso_rolling with regular lasso

#combine 3D panels into one dataframe, consisting of the edge and average rank through windows. rank standard deviation.


#final_model.sort_index(axis=0)



#load gold standard model and compare averaged model to gold stardard with AUROC

