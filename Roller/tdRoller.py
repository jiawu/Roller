__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
import numpy as np
import pandas as pd

class tdRoller(Roller):
    """
    A roller that incorporates time delays
    """

    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator="\t",
                 window_type="RandomForest"):

        super(tdRoller, self).__init__(file_path, gene_start, gene_end, time_label, separator,
                 window_type)

if __name__ == "__main__":
    file_path = "../data/dream4/insilico_size10_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.zscore_all_data()
