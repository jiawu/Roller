__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
import numpy as np
import pandas as pd
from scipy import stats
import sys

class tdRoller(Roller):
    """
    A roller that incorporates time delays
    """

    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator="\t",
                 window_type="RandomForest", threshold=1, q=None):

        super(tdRoller, self).__init__(file_path, gene_start, gene_end, time_label, separator,
                 window_type)

        # Zscore the data
        self.norm_data = self.zscore_all_data()

        # transform the matrix
        self.binned_data = self.bin_data()

    def zscore_all_data(self):
        """
        Z-score the data by column
        :return:
        """
        # Get the raw data values
        raw_dataset = self.raw_data.values.copy()

        # z-score the values with 1 degree of freedom
        zscored_datset = pd.DataFrame(stats.zscore(raw_dataset, axis=0, ddof=1), index=self.raw_data.index,
                                      columns=self.raw_data.columns)

        # replace time values with original time values
        zscored_datset[self.time_label] = self.raw_data[self.time_label]
        return zscored_datset

    def bin_data(self, threshold=1):
        """
        Transform the data into O'' then O' using the method described in Liping et. al, Bioinformatics, 2005

        NOTE: The data here is in the form m x n (timepoints by genes) which is the transpose of how it is presented
        in the original paper

        :param threshold: float
            the normalization threshold
        :return:
        """
        # Just to be safe make the threshold positive
        threshold = np.abs(threshold)

        # The last datapoint cannot be assessed because there is no t+1 for it
        o_data = self.norm_data.values.copy() # Ignore the time data
        o_i_j = o_data[:-1].copy()
        o_i_j_plusone = o_data[1:].copy()

        o_doubleprime = (o_i_j_plusone - o_i_j)/np.abs(o_i_j)

        # If o_i_j is 0 then there are three possibilities that are easily fixes
        #   1. o_i_j_plusone is also zero. This will result in a value of nan from numpy. Replace with 0
        #   2. o_i_j_plusone is > zero. This will result in a value of inf from numpy. Replace with 1
        #   3. o_i_j_plusone is < zero. This will result in a value of -inf from numpy. Replace with -1
        o_doubleprime[np.isnan(o_doubleprime)] = 0
        o_doubleprime[o_doubleprime == np.inf] = 1
        o_doubleprime[o_doubleprime == -np.inf] = -1

        # Threshold values
        o_prime = np.zeros(o_doubleprime.shape)
        o_prime[o_doubleprime >= threshold] = 1
        o_prime[o_doubleprime <= -threshold] = -1

        df = pd.DataFrame(o_prime, index=self.raw_data.index[:-1], columns=self.raw_data.columns)
        df[self.time_label] = self.raw_data[self.time_label]

        return df

    def create_q_clusters(self, q=None):
        if q is None:
            q = int(self.overall_width/2)

        self.set_window(q)
        self.create_windows()
        print self.window_list[0].raw_data

    def create_windows(self, random_time=False, ):
        """
        Create window objects for the roller to use

        Called by:
            pipeline

        :return:
        """
        window_list = [self.get_window_object(self.get_window_raw(index, random_time),
                                              {"time_label": self.time_label,
                                               "gene_start": self.gene_start,
                                               "gene_end": self.gene_end,
                                               "nth_window": index}) if (
        index + self.window_width <= self.overall_width) else '' for index in range(self.get_n_windows())]
        self.window_list = window_list

    def get_window_raw(self, start_index, random_time=False):
        """
        Select a window from the full data set. This is fancy data-frame slicing

        Called by:
            create_windows
            get_window_stats
            get_window

        :param start_index: int
            The start of the window
        :param random_time: bool, optional
        :return: data-frame
        """
        if random_time:
            # select three random timepoints
            time_window = self.time_vec[start_index]
            choices = self.time_vec
            choices = np.delete(choices, start_index)
            for x in range(0, self.window_width - 1):
                chosen_time = random.choice(choices)
                time_window = np.append(time_window, chosen_time)
                chosen_index = np.where(choices == chosen_time)
                choices = np.delete(choices, chosen_index)
        else:
            end_index = start_index + self.window_width
            time_window = self.time_vec[start_index:end_index]
        data = self.raw_data[self.raw_data[self.time_label].isin(time_window)]
        return data

if __name__ == "__main__":
    file_path = "../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.create_q_clusters()