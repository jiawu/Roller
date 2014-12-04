__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import Roller
from sklearn.preprocessing import Imputer
import matplotlib as mpl
import numpy as np
import pandas as pd

if __name__ == '__main__':
    file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)


    window_size = 6
    #get only TFs data, window size of 4
    roll_me.set_window(window_size)
    #impute missing values
    imputer = Imputer(missing_values="NaN")
    mpl.rcParams.update({'font.size':8})
    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    total_window_number = roll_me.get_n_windows()
    #todo: fix the below two lines, I don't need this to get the matrix size...
    current_window = roll_me.get_window()
    print current_window.shape
    filled_matrix = imputer.fit_transform(current_window)
    n_genes = filled_matrix.shape[1]
    coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))
    gene_names=list(current_window.columns.values)
    print filled_matrix.shape