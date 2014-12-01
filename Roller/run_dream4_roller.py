from Roller import Roller
from sklearn.preprocessing import Imputer
from linear_wrapper import LassoWrapper
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import matplotlib as mpl
import pdb
#import auroc

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot_lines(series_list, regulator_labels, target_labels, window_size, suffix=""):
    figure2 = plt.figure()
    lineplot = figure2.add_subplot(1,1,1)
    lineplot.set_xlabel('start day')
    lineplot.set_ylabel('Beta')
    lines = []
    time = [x for x in range(0,22-window_size)]
    label_list = []
    for counter,series in enumerate(series_list):
        my_label = str(regulator_labels[counter]+" -> "+target_labels[counter])
        label_list.append(my_label)
        pdb.set_trace()
        line, = lineplot.plot(time, series, label = my_label)
        lines.append(line)
    figure2.legend(lines,label_list)
    figure2.savefig('line_figure'+str(window_size)+suffix+'.png')

def plot_figure(coeff_matrix,nth_window, row_labels, col_labels, window_size,prefix=""):
    df = pd.DataFrame(coeff_matrix)
    figure1 = plt.figure()
    heatmap = figure1.add_subplot(1,1,1)
    my_axis = heatmap.imshow(df,interpolation='nearest',cmap=cm.OrRd, vmin=0, vmax=1.2)
    my_axi = my_axis.get_axes()
    clean_axis(my_axi)
    heatmap.set_yticks(np.arange(df.shape[0]))
    heatmap.yaxis.set_ticks_position('left')
    heatmap.set_yticklabels(row_labels)
    heatmap.set_xticks(np.arange(df.shape[1]))
    heatmap.xaxis.set_ticks_position('top')
    xlabelsL = heatmap.set_xticklabels(col_labels)
    for label in xlabelsL:
        label.set_rotation(90)
    for l in heatmap.get_xticklines() + heatmap.get_yticklines():
        l.set_markersize(0)
    title=heatmap.set_title("Day " +str(nth_window) + " to Day " +str(nth_window+window_size))
    title.set_x(1.2)
    figure1.savefig(prefix+'figure'+str(nth_window)+'.png')


#file_path = "/Users/jjw036/Roller/experiments/dream5_experiment/raw_data/insilico_size10_1_timeseries.tsv"
file_path = "/Users/jjw036/Roller/goldbetter_model/goldbetter_data.txt"
gene_start_column = 1
roll_me = Roller(file_path, gene_start_column)
window_size = 1000
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
coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))
gene_names=list(current_window.columns.values)

for nth_window in range(0, total_window_number):
    #loop gets the window, gets the coefficients for that window, then increments the window
    current_window = roll_me.get_window()

    #check if any values are completely blank
    sums = [current_window.iloc[:,x].sum() for x in range(0,10)]
    ind = np.where(np.isnan(sums))[0]
    current_window.iloc[:,ind]=0
    current_window = current_window *100
    filled_matrix = current_window.values
    #filled_matrix = imputer.fit_transform(current_window)
    current_lasso = LassoWrapper(filled_matrix)
    coeff_mat = current_lasso.get_coeffs(10)
    coeff_matrix_3d[:,:,nth_window] = coeff_mat
    #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
    roll_me.next()

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
    plot_figure(coeff_mat,nth_window,present_genes_ii, present_genes_jj, window_size,str("compressed_window"+str(window_size)))

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
    pdb.set_trace()
    np.savetxt('binary_mat')
    #plot_lines(tf_list, regulator_labels, target_labels, window_size, suffix=str("targets_of_"+gene))

