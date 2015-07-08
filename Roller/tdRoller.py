__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
from RFRWindow import tdRFRWindow
import numpy as np
import pandas as pd
from scipy import stats
import sys
import random

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

if __name__ == "__main__":
    file_path = "../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.zscore_all_data()