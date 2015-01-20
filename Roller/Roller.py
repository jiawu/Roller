import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
import numpy as np
from util.linear_wrapper import LassoWrapper
import Window
from LassoWindow import LassoWindow


class Roller(object):
    """
    A thing that grabs different timepoints of data, can set window and step size.

    To do list:
        -i need to make two modules: a file processing module and a rolling module.
        -accept different table formats

        -add permute_window()
        -add bootstrape_window()
    """
    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator = "\t"):
        """
        Initialize the roller object. Read the file and put it into a pandas dataframe
        :param file_path: file-like object or string
                        The file to read
        :param gene_start: int
        :param gene_end: int
        """
        # Read the raw data into a pandas dataframe object
        self.raw_data = pd.read_csv(file_path, sep=separator)
        self.raw_data = self.raw_data.dropna(axis=0, how='all')

        # Set roller defaults
        self.current_step = 0
        self.window_width = 3
        self.step_size = 1
        self.time_label = time_label

        # Get overall width of the time-course
        self.time_vec = self.raw_data[self.time_label].unique()
        self.overall_width = len(self.time_vec)

        if gene_end is not None:
            self.gene_end = gene_end
        else:
            self.gene_end = len(self.raw_data.columns)
        if gene_start is not None:
            self.gene_start = gene_start
        else:
            self.gene_start = 0

        self.gene_list = self.raw_data.columns.values[self.gene_start:self.gene_end]

        self.current_window = self.get_window(self.current_step)
        self.window_list = None

    def get_n_windows(self):
        total_windows = (self.overall_width - self.window_width+1)/(self.step_size)
        return total_windows

    def get_window(self, start_index):
        raw_window = self.get_window_raw(start_index)
        only_genes = raw_window.iloc[:, self.gene_start:self.gene_end]
        return only_genes

    def get_window_raw(self, start_index):
        end_index = start_index + self.window_width
        time_window = self.time_vec[start_index:end_index]
        data = self.raw_data[self.raw_data[self.time_label].isin(time_window)]
        return data

    def set_window(self, width):
        self.window_width = width

    def set_step(self, step):
        self.step_size = step

    def reset(self):
        self.current_step = 0

    def remove_blank_rows(self):
        """calculates sum of rows. if sum is NAN, then remove row"""
        coln = len(self.raw_data.columns)
        sums = [self.raw_data.iloc[:,x].sum() for x in range(0,coln)]
        ind = np.where(np.isnan(sums))[0]
        self.raw_data.iloc[:,ind]=0

    def get_n_genes(self):
        return(len(self.raw_data.columns) -1)

    def create_windows(self):
        window_list = [LassoWindow(self.get_window(index)) if (index + self.window_width <= self.overall_width) else ''
                       for index in range(self.get_n_windows())]
        self.window_list = window_list

    def initialize_windows(self):
        for window in self.window_list:
            window.initialize_params()

    def zscore_all_data(self):
        #zscores all the data
        dataframe =  self.raw_data

        #for each column, calculate the zscore
        #zscore is calculated as X - meanX / std(ddof = 1)
        for item in dataframe.columns:
            if item != self.time_label:
                dataframe[item] = (dataframe[item] - dataframe[item].mean())/dataframe[item].std(ddof=1)
        self.raw_data = dataframe

    def get_window_stats(self):
        """for each window, get a dict:
            N : the number of datapoints in this window,
            time_labels: the names of the timepoints in a roller model
            step_size: the step-size of the current model
            window_size: the size of the window of the current model
            total_windows: the number of windows total
            window_index: the index of the window. counts start at 0. ie if the window index is 0 it is the 1st window. if the window index is 12, it is the 12th window in the series."""
        current_window = self.get_window_raw()

        """calculate the window index. todo: move into own function later"""
        min_time = np.amin(current_window[self.time_label])
        window_index = np.where(self.time_vec == min_time)/self.step_size
        # to calculate the nth window, time vector
        # index of the time-vector, step size of 2? window 4, step size 2
        #
        #total windows = total width (10) - window_width (2) +1 / step size
        # 10 time points 0 1 2 3 4 5 6 7 8 9
        #width is 2: 0 and 1
        # step size is 2
        # 01, 12, 23, 34, 45, 56, 67, 78, 89

        #todo: so the issue is that total windows (get n windows) is the true number of windows, and window index is the nth -1 window... it would be great to consolidate these concepts but no big deal if they can't be.


        window_stats = {'N': len(current_window.index),
                        'time_labels': current_window[self.time_label].unique(),
                        'step_size': self.step_size,
                        'window_size': self.window_width,
                        'total_windows': self.get_n_windows(),
                        'window_index': window_index}
        return window_stats


