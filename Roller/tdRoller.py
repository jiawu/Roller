__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
from RFRWindow import tdRFRWindow
import numpy as np
import pandas as pd
from scipy import stats
import sys
import random
import matplotlib.pyplot as plt
from util.Evaluator import Evaluator

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

    def create_windows(self, random_time=False):
        """
        Create window objects for the roller to use

        Called by:
            pipeline

        :return:
        """
        window_list = [self.get_window_object(self.get_window_raw(index, random_time), index,
                                              {"time_label": self.time_label+'_'+str(index),
                                               "gene_start": self.gene_start,
                                               "gene_end": self.gene_end,
                                               "nth_window": index}) if (
        index + self.window_width <= self.overall_width) else '' for index in range(self.get_n_windows())]
        self.window_list = window_list

    def get_window_object(self, dataframe, index, window_info_dict):
        """
        Return a window object from a data-frame

        Called by:
            create_windows

        :param dataframe: data-frame
        :param window_info_dict: dict
            Dictionary containing information needed for window initialization
        :return:
        """

        dataframe.columns = [col+'_'+str(index) for col in dataframe.columns]
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
        # todo: should only allow for regression against earlier timepoints? (Testing seems to indicate no, just correct afterward)
        # Enumerate all of the windows except the first
        for ww, window in enumerate(self.window_list):
            if ww == 0:
                window.x_data = window.window_values.copy()
                window.x_labels = window.raw_data.columns[1:]
            else:
                window.x_data = np.hstack((window.window_values, self.window_list[ww-1].x_data))
                window.x_labels = np.append(window.raw_data.columns[1:], self.window_list[ww-1].x_labels)
            window.augmented_edge_list = window.possible_edge_list(window.x_labels, window.raw_data.columns[1:])

def point_slope(x1,y1, x2,y2):
    slope = (y2-y1)/float(x2-x1)
    return slope

def elbow_criteria(x,y):
    x = np.array(x)
    y = np.array(y)
    # Slope between elbow endpoints
    m1 = point_slope(x[0], y[0], x[-1], y[-1])
    # Intercept
    b1 = y[0] - m1*x[0]

    # Slope for perpendicular lines
    m2 = -1/m1

    # Calculate intercepts for perpendicular lines that go through data point
    b_array = y-m2*x
    x_perp = (b_array-b1)/(m1-m2)
    y_perp = m1*x_perp+b1

    # Calculate where the maximum distance to a line connecting endpoints is
    distances = np.sqrt((x_perp-x)**2+(y_perp-y)**2)
    index_max = np.where(distances==np.max(distances))[0][0]
    elbow_x = x[index_max]
    elbow_y = y[index_max]
    return elbow_x, elbow_y

if __name__ == "__main__":
    file_path = "../data/dream4/yeast_size100_1_timeseries.tsv"
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    #pd.set_option('display.width', 1000)

    np.random.seed(8)

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.zscore_all_data()
    tdr.set_window(5)
    tdr.create_windows()
    tdr.augment_windows()
    tdr.fit_windows(n_trees=5000)
    #todo: functionalize and speed up this part
    full_edge_list = []
    full_edge_importance = []
    print 'lumping edges'
    for window in tdr.window_list[1:]:
        full_edge_importance += list(window.edge_importance.flatten())
        full_edge_list += window.augmented_edge_list
    print len(full_edge_importance)
    print len(full_edge_list)
    print 'reverse zip'
    parents, children = zip(*full_edge_list)
    df = pd.DataFrame([list(parents), list(children), full_edge_importance], index=['Parent', 'Child', 'Importance']).T
    print "sort"
    df.sort(columns='Importance', ascending=False, inplace=True)
    print 'split'
    a = df['Parent'].str.split('_').apply(pd.Series,1)
    print 'split'
    b = df['Child'].str.split('_').apply(pd.Series,1)
    df['Parent'] = a.iloc[:,0]
    df['P_window'] = a.iloc[:,1]
    df['Child'] = b.iloc[:,0]
    df['C_window'] = b.iloc[:,1]
    #df = df[df.P_window != df.C_window]
    df = df[df.Parent != df.Child]
    df['Edge'] = zip(df.Parent, df.Child)
    #print df
    #print len(df)
    edge_set = list(set(df.Edge))
    x, y = elbow_criteria(range(len(df.Importance)), df.Importance.values.astype(np.float64))
    #df = df[df.Importance>y]

    print 'calc edge imp'
    edge_imp_vals = {edge:df.Importance[df.Edge==edge] for edge in edge_set}
    edge_mean_importance = [np.mean(df.Importance[df.Edge==edge]) for edge in edge_set]
    #edge_mean_importance = [np.std(df.Importance[df.Edge==edge], ddof=1) for edge in edge_set]

    df2 = pd.DataFrame([edge_set, edge_mean_importance],
                       index=['regulator-target', 'mean_imp']).T
    print 'second sort'
    df2.sort(columns='mean_imp', ascending=False, inplace=True)
    #ranked_edge_list = df2['Edge'].tolist()
    current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
    evaluator = Evaluator(current_gold_standard, '\t')
    tpr, fpr, auroc = evaluator.calc_roc(df2)
    print "mean", auroc[-1]#+(1-fpr[-1])
    #plt.plot(fpr,tpr)
    #plt.plot(fpr, fpr)
    #plt.show()
    sys.exit()
    #plt.figure()
    nbins = 15

    top_n = len(df2)
    #print df2.head(top_n)
    ks_val = []
    gene_list = []
    f, axarr = plt.subplots(top_n)
    for edge_n in range(top_n):
        #axarr[edge_n].hist(edge_imp_vals[df2['regulator-target'].iloc[edge_n]].values, bins=nbins, alpha=0.5, label='->'.join(df2['regulator-target'].iloc[edge_n]))
        data = edge_imp_vals[df2['regulator-target'].iloc[edge_n]].values.astype(np.float64)
        x = np.linspace(0, np.max(data))
        params = stats.weibull_min.fit(data)
        #Parameters returned are shape, location, scale
        #print params
        #print np.mean(data)
        #print df2['regulator-target'].iloc[edge_n]
        # print 2*params[1]*np.sqrt(2/np.pi)
        #print stats.kstest(data, 'norm', args=params)
        gene_list.append(df2['regulator-target'].iloc[edge_n])
        ks_val.append(stats.kstest(data, 'weibull_min', args=params)[0])
        #axarr[edge_n].plot(x, stats.norm.pdf(x, *params), lw=3)
        #axarr[edge_n].legend(loc='upper right')
    df3 =pd.DataFrame([gene_list, ks_val], index=['regulator-target', 'KS']).T
    df3.sort(columns='KS', inplace=True)
    print df3.head(10)
    tpr, fpr, auroc = evaluator.calc_roc(df3)
    print "mean", auroc[-1]
    # plt.hist(edge_imp_vals[df2['regulator-target'].iloc[1]].values, bins=nbins, alpha=0.5, label='->'.join(df2['regulator-target'].iloc[1]))
    # plt.hist(edge_imp_vals[df2['regulator-target'].iloc[2]].values, bins=nbins, alpha=0.5, label='->'.join(df2['regulator-target'].iloc[2]))
    # plt.hist(edge_imp_vals[df2['regulator-target'].iloc[3]].values, bins=nbins, alpha=0.5, label='->'.join(df2['regulator-target'].iloc[3]))
    # plt.hist(edge_imp_vals[df2['regulator-target'].iloc[4]].values, bins=nbins, alpha=0.5, label='->'.join(df2['regulator-target'].iloc[4]))

    #print df2['regulator-target'].iloc[0]
    #plt.show()

