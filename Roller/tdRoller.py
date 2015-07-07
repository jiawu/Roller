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

if __name__ == "__main__":
    file_path = "../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)