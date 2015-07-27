__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
from RFRWindow import tdRFRWindow
import numpy as np
import pandas as pd
from scipy import stats
import sys
import warnings
import random
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
        self.lag_set = None

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

    def augment_windows(self, min_lag=0, max_lag=None):
        """
        Window data is augmented to include data from previous time points and labeled accordingly
        :param min_lag: int
        :param max_lag: int or None
            if None, all earlier windows will included
        :return:
        """
        if min_lag>max_lag and max_lag is not None:
            raise Exception('The minimum lag cannot be greater than the maximum lag')
        # Enumerate all of the windows except the first
        for window in self.window_list:
            window_idx = window.nth_window
            if max_lag is None:
                start_idx = 0
            else:
                start_idx = max(window_idx-max_lag,0)
            end_index = max(window_idx-min_lag+1, 0)
            if max_lag is not None and (end_index-start_idx)<(max_lag-min_lag+1):
                window.include_window = False
                continue

            earlier_windows = self.window_list[start_idx:end_index]
            window.earlier_window_idx = [w.nth_window for w in earlier_windows]
            # Add necessary data from earlier windows
            for ww, win in enumerate(earlier_windows[::-1]): #Go through the list in reverse because of how the window expects data
                if ww == 0:
                    #Initialize values
                    window.x_data = win.window_values.copy()
                    window.x_labels = win.raw_data.columns[1:]
                    window.x_times = np.array([win.nth_window]*len(win.genes))
                else:
                    window.x_data = np.hstack((window.x_data, win.window_values))
                    window.x_labels = np.append(window.x_labels, win.raw_data.columns[1:])
                    window.x_times = np.append(window.x_times, np.array([win.nth_window]*len(win.genes)))
            if window.x_data is None:
                window.include_window = False
        return

    def compile_roller_edges(self, self_edges=False):
        """
        Edges across all windows will be compiled into a single edge list
        :return:
        """
        print "Compiling all model edges...",
        df = None
        for ww, window in enumerate(self.window_list):
            if window.include_window:
                if df is None:
                    df = window.make_edge_table()
                else:
                    df = df.append(window.make_edge_table(), ignore_index=True)

        if not self_edges:
            df = df[df.Parent != df.Child]

        df['Edge'] = zip(df.Parent, df.Child)
        df['Lag'] = df.C_window - df.P_window
        self.full_edge_list = df.copy()
        print "[DONE]"
        return

    def make_static_edge_dict(self, true_edges, self_edges=False, lag_method='max_median'):
        """
        Make a dictionary of edges
        :return:
        """
        print "Lumping edges...",
        df = self.full_edge_list.copy()
        #Only keep edges with importance > 0
        df = df[df['Importance']>0]
        if not self_edges:
            df = df[df.Parent != df.Child]
        edge_set = list(set(df.Edge))
        full_edge_set = set(self.make_possible_edge_list(self.gene_list, self.gene_list, self_edges=self_edges))
        edge_diff = full_edge_set.difference(edge_set)
        self.edge_dict = {}
        lag_importance_score, lag_lump_method = lag_method.split('_')
        score_method = eval('np.'+lag_importance_score)
        lump_method = eval('np.'+lag_lump_method)
        for edge in full_edge_set:
            if edge in edge_diff:
                self.edge_dict[edge] = {"dataframe": None, "mean_importance": 0, 'real_edge': (edge in true_edges),
                                        "max_importance": 0, 'max_edge': None, 'lag_importance': 0,
                                        'lag_method': lag_method}
                continue
            current_df = df[df.Edge==edge]
            max_idx = current_df['Importance'].idxmax()
            lag_set = list(set(current_df.Lag))
            lag_imp = score_method([lump_method(current_df.Importance[current_df.Lag==lag]) for lag in lag_set])
            self.edge_dict[edge] = {"dataframe":current_df, "mean_importance":np.mean(current_df.Importance),
                                    'real_edge':(edge in true_edges), "max_importance":current_df.Importance[max_idx],
                                    'max_edge':(current_df.P_window[max_idx], current_df.C_window[max_idx]),
                                    'lag_importance': lag_imp, 'lag_method':lag_method}
        print "[DONE]"
        if edge_diff:
            message = 'The last %i edges had no meaningful importance score' \
                      ' and were placed at the bottom of the list' %len(edge_diff)
            warnings.warn(message)
        return

    def make_sort_df(self, df, sort_by='mean'):
        """
        Calculate the mean for each edge
        :param df: dataframe
        :return: dataframe
        """
        sort_field = sort_by+"_importance"
        print "Calculating %s edge importance..." %sort_by,
        temp_dict = {edge:df[edge][sort_field] for edge in df.keys()}
        sort_df = pd.DataFrame.from_dict(temp_dict, orient='index')
        sort_df.columns = [sort_field]
        sort_df.sort(sort_field, ascending=False, inplace=True)
        #sort_df['mean_importance'] = stats.zscore(sort_df['mean_importance'], ddof=1)
        sort_df.index.name = 'regulator-target'
        sort_df=sort_df.reset_index()
        print "[DONE]"
        return sort_df

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
        print "Scoring model...",
        if gold_standard_file is None:
            current_gold_standard = self.file_path.replace("timeseries.tsv","goldstandard.tsv")
        else:
            current_gold_standard = gold_standard_file
        evaluator = Evaluator(current_gold_standard, '\t')
        tpr, fpr, auroc = evaluator.calc_roc(sorted_edge_list)
        auroc_dict = {'tpr':np.array(tpr), 'fpr':np.array(fpr), 'auroc': np.array(auroc)}
        precision, recall, aupr = evaluator.calc_pr(sorted_edge_list)
        aupr_random = [len(evaluator.gs_flat)/float(len(evaluator.full_list))]*len(recall)
        aupr_dict = {"precision": np.array(precision), "recall": np.array(recall), "aupr": np.array(aupr),
                     "aupr_random": np.array(aupr_random)}
        print "[DONE]"
        return auroc_dict, aupr_dict

    def plot_scoring(self, roc_dict, pr_dict):
        f, axarr = plt.subplots(2, figsize=(7,10))
        axarr[0].plot(roc_dict['fpr'], roc_dict['tpr'])
        axarr[0].plot(roc_dict['fpr'], roc_dict['fpr'])
        axarr[0].legend(['tdRoller', 'Random'], loc='best')
        axarr[0].set_xlabel('False Positive Rate')
        axarr[0].set_ylabel('True Positive Rate')
        title = 'ROC Curve (AUROC = %0.3f)' %roc_dict['auroc'][-1]
        axarr[0].set_title(title)
        axarr[1].plot(pr_dict['recall'], pr_dict['precision'])
        axarr[1].plot(pr_dict['recall'], pr_dict['aupr_random'])
        axarr[1].legend(['tdRoller', 'Random'], loc='best')
        axarr[1].set_xlabel('Recall')
        axarr[1].set_ylabel('Precision')
        title = 'PR Curve (AUPR = %0.3f)' %pr_dict['aupr'][-1]
        axarr[1].set_title(title)
        plt.tight_layout()
        plt.show()
        return

