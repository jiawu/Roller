__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'
import pdb
from Roller import Roller
from RFRWindow import tdRFRWindow
from LassoWindow import tdLassoWindow
from DionesusWindow import tdDionesusWindow
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
        if self.window_type == "RandomForest":
            window_obj = tdRFRWindow(dataframe, window_info_dict, self.norm_data)
        elif self.window_type =="Lasso":
            window_obj = tdLassoWindow(dataframe, window_info_dict, self.norm_data)

        elif self.window_type =="Dionesus":
            window_obj = tdDionesusWindow(dataframe, window_info_dict, self.norm_data)

        return window_obj

    def augment_windows(self, min_lag=0, max_lag=None):
        """
        Window data is augmented to include data from previous time points and labeled accordingly
        :param min_lag: int
        :param max_lag: int or None
            if None, all earlier windows will included
        :return:
        """
        new_window_list = []

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
                    #try not to assign these outside the classfile. I could not figure out where these were assigned (they weren't in the file)
                    window.x_labels = win.data.columns[1:]
                    window.x_times = np.array([win.nth_window]*len(win.genes))

                else:
                    window.x_data = np.hstack((window.x_data, win.window_values))
                    window.x_labels = np.append(window.x_labels, win.data.columns[1:])
                    window.x_times = np.append(window.x_times, np.array([win.nth_window]*len(win.genes)))
            if window.x_data is not None:
                new_window_list.append(window)
            else:
                window.include_window = False
        self.window_list = new_window_list
        return

    def compile_roller_edges(self, self_edges=False, mse_adjust=True):
        """
        Edges across all windows will be compiled into a single edge list
        :return:
        """
        print "Compiling all model edges...",
        df = None
        for ww, window in enumerate(self.window_list):
            if window.include_window:
                # Get the edges and associated values in table form
                current_df = window.make_edge_table()

                # Only retain edges if the p_value is below the threshold
                #current_df = current_df[current_df['p_value'] <= 0.05]

                # Only retain edges if the MSE_diff is negative
                if mse_adjust:
                    current_df = current_df[current_df['MSE_diff'] < 0]


                current_df['adj_imp'] = current_df['Importance']*(1-current_df['p_value'])#*current_df['MSE_diff']

                current_df.sort(['adj_imp'], ascending=False, inplace=True)
                current_df['Rank'] = np.arange(len(current_df))

                if df is None:
                    df = current_df.copy()
                else:
                    df = df.append(current_df.copy(), ignore_index=True)
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
        print
        df = self.full_edge_list.copy()

        # Only keep edges with importance > 0. Values below 0 are not helpful for model building
        df = df[df['Importance'] > 0]

        # Ignore self edges if desired
        if not self_edges:
            df = df[df.Parent != df.Child]
        edge_set = list(set(df.Edge))

        # Calculate the full set of potential edges
        full_edge_set = set(self.make_possible_edge_list(self.gene_list, self.gene_list, self_edges=self_edges))

        # Identify edges that could exist, but do not appear in the inferred list
        edge_diff = full_edge_set.difference(edge_set)

        self.edge_dict = {}
        lag_importance_score, lag_lump_method = lag_method.split('_')
        score_method = eval('np.'+lag_importance_score)
        lump_method = eval('np.'+lag_lump_method)
        for edge in full_edge_set:
            if edge in edge_diff:
                self.edge_dict[edge] = {"dataframe": None, "mean_importance": 0, 'real_edge': (edge in true_edges),
                                        "max_importance": 0, 'max_edge': None, 'lag_importance': 0,
                                        'lag_method': lag_method, 'rank_importance':np.nan}
                continue
            current_df = df[df.Edge == edge]
            max_idx = current_df['Importance'].idxmax()
            lag_set = list(set(current_df.Lag))
            lag_imp = score_method([lump_method(current_df.Importance[current_df.Lag == lag]) for lag in lag_set])
            lag_rank = score_method([lump_method(current_df.Rank[current_df.Lag == lag]) for lag in lag_set])
            self.edge_dict[edge] = {"dataframe":current_df, "mean_importance":np.mean(current_df.Importance),
                                    'real_edge':(edge in true_edges), "max_importance":current_df.Importance[max_idx],
                                    'max_edge':(current_df.P_window[max_idx], current_df.C_window[max_idx]),
                                    'lag_importance': lag_imp, 'lag_method':lag_method,
                                    'rank_importance': lag_rank}
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
        if sort_by.lower() == 'rank':
            sort_df.sort(sort_field, ascending=True, inplace=True)
        else:
            sort_df.sort(sort_field, ascending=False, inplace=True)
        #sort_df['mean_importance'] = stats.zscore(sort_df['mean_importance'], ddof=1)
        sort_df.index.name = 'regulator-target'
        sort_df = sort_df.reset_index()
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

    def get_samples(self):
        df=pd.read_csv(self.file_path,sep='\t')
        node_list = df.columns.tolist()
        node_list.pop(0)
        return(node_list)

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

        evaluator = Evaluator(current_gold_standard, '\t', node_list=self.get_samples())
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

