__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
from RFRWindow import tdRFRWindow
import numpy as np
import pandas as pd
from scipy import stats
import sys
import itertools
import random
import time
import matplotlib.pyplot as plt
from util.Evaluator import Evaluator
from util.utility_module import elbow_criteria

class tdRoller(Roller):
    """
    A roller that incorporates time delays
    """

    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator="\t",
                 window_type="RandomForest", threshold=1, q=None):

        super(tdRoller, self).__init__(file_path, gene_start, gene_end, time_label, separator,
                 window_type)
        self.full_edge_list = None
        self.edge_dict = None

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

    def create_windows(self, random_time=False):
        """
        Create window objects for the roller to use

        Called by:
            pipeline

        :return:
        """
        window_list = [self.get_window_object(self.get_window_raw(index, random_time),
                                              {"time_label": self.time_label+'_'+str(index),
                                               "gene_start": self.gene_start,
                                               "gene_end": self.gene_end,
                                               "nth_window": index}) if (
        index + self.window_width <= self.overall_width) else '' for index in range(self.get_n_windows())]
        self.window_list = window_list

    def get_window_object(self, dataframe, window_info_dict):
        """
        Return a window object from a data-frame

        Called by:
            create_windows

        :param dataframe: data-frame
        :param window_info_dict: dict
            Dictionary containing information needed for window initialization
        :return:
        """
        window_obj = tdRFRWindow(dataframe, window_info_dict, self.raw_data)

        return window_obj

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
        data = self.norm_data[self.norm_data[self.time_label].isin(time_window)]
        return data

    def augment_windows(self):
        """
        Window data is augmented to include data from previous time points and labeled accordingly
        :return:
        """
        # todo: should only allow for regression against earlier timepoints? (Testing seems to indicate no, just correct afterward)
        # Enumerate all of the windows except the first
        for ww, window in enumerate(self.window_list):
            if ww == 0:
                window.x_data = window.window_values.copy()
                window.x_labels = window.raw_data.columns[1:]
                window.x_times = np.array([window.nth_window]*len(window.genes))

            else:
                window.x_data = np.hstack((window.window_values, self.window_list[ww-1].x_data))
                window.x_labels = np.append(window.raw_data.columns[1:], self.window_list[ww-1].x_labels)
                window.x_times = np.append(np.array([window.nth_window]*len(window.genes)), self.window_list[ww-1].x_times)

    def compile_roller_edges(self, self_edges=False):
        """
        Edges across all windows will be compiled into a single edge list
        :return:
        """
        print "Compiling all model edges...",
        for window in self.window_list:

            if window.nth_window == 0:
                df = window.make_edge_table()
            else:
                df = df.append(window.make_edge_table(), ignore_index=True)

        if not self_edges:
            df = df[df.Parent != df.Child]

        df['Edge'] = zip(df.Parent, df.Child)
        self.full_edge_list = df.copy()
        print "[DONE]"
        return

    def make_static_edge_dict(self, self_edges=False):
        """
        Make a dictionary of edges
        :return:
        """
        df = self.full_edge_list.copy()
        if not self_edges:
            df = df[df.Parent != df.Child]
        edge_set = list(set(df.Edge))
        self.edge_dict = {}

        for edge in edge_set:
            current_df = df[df.Edge==edge]
            self.edge_dict[edge] = {"dataframe":current_df, "mean_importance":np.mean(current_df.Importance)}

        return

    def calc_edge_mean(self, df):
        """
        Calculate the mean for each edge
        :param df: dataframe
        :return: dataframe
        """
        temp_dict = {edge:df[edge]['mean_importance'] for edge in df.keys()}
        mean_df = pd.DataFrame.from_dict(temp_dict, orient='index')
        mean_df.columns = ['mean_importance']
        mean_df.sort('mean_importance', ascending=False, inplace=True)
        mean_df['regulator-target'] = mean_df.index
        return mean_df

    def calc_edge_importance_cutoff(self, df):
        """
        Calculate the importance threshold to filter edges on
        :param df:
        :return: dict
        """
        x, y = elbow_criteria(range(len(df.Importance)), df.Importance.values.astype(np.float64))
        elbow_dict = {'num_edges':x, 'importance_threshold':y}

        return elbow_dict

    def score(self, sorted_edge_list, gold_standard_file=None):
        """
        Scores some stuff, I think...
        Called by:
            pipeline
        :param sorted_edge_list:
        :param gold_standard_file:
        :return:
        """
        if gold_standard_file is None:
            current_gold_standard = self.file_path.replace("timeseries.tsv","goldstandard.tsv")
        else:
            current_gold_standard = gold_standard_file
        evaluator = Evaluator(current_gold_standard, '\t')
        tpr, fpr, auroc = evaluator.calc_roc(sorted_edge_list)
        auroc_dict = {'tpr':tpr, 'fpr':fpr, 'auroc': auroc}
        precision, recall, aupr = evaluator.calc_pr(sorted_edge_list)
        aupr_random = [len(evaluator.gs_flat)/float(len(evaluator.full_list))]*len(recall)
        aupr_dict = {"precision": precision, "recall": recall, "aupr": aupr, "aupr_random":aupr_random}
        return auroc_dict, aupr_dict

if __name__ == "__main__":
    file_path = "../data/dream4/insilico_size10_2_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    #pd.set_option('display.width', 1000)

    #np.random.seed(8)

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.zscore_all_data()
    tdr.set_window(13)
    tdr.create_windows()
    tdr.augment_windows()
    tdr.fit_windows(n_trees=10)
    tdr.compile_roller_edges(self_edges=True)
    tdr.make_static_edge_dict()
    df2 = tdr.calc_edge_mean(tdr.edge_dict)
    roc_dict, pr_dict = tdr.score(df2)
    f, axarr = plt.subplots(2)
    axarr[0].plot(roc_dict['fpr'], roc_dict['tpr'])
    axarr[0].plot(roc_dict['fpr'], roc_dict['fpr'])
    axarr[1].plot(pr_dict['recall'], pr_dict['precision'])
    axarr[1].plot(pr_dict['recall'], pr_dict['aupr_random'])
    plt.show()


    sys.exit()

