import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import unittest
import numpy as np
import Roller
import pdb
import Roller.util.Grapher as Grapher
import pandas as pd
from Roller.util.Evaluator import Evaluator

class TestRFGrapher(unittest.TestCase):
    def setUp(self):
        file_path = "../../data/dream4/insilico_size10_1_timeseries.tsv"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        self.roller = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator,window_type = "RandomForest")
        roller_result_path = "sample_RF_roller_results.pickle"
        self.roller_results = pd.read_pickle(roller_result_path)
    def test_heatmap(self):
        self.roller.set_window(20)
        self.roller.create_windows_no_next()
        raw_data=self.roller.raw_data
        figure,ax1 = Grapher.generate_heatmap_from_df(raw_data)

    def test_combo_graphs(self):
        roller = self.roller_results
        raw_data = roller.raw_data
        window_size = roller.window_width
        window_start_list = []
        auroc_list = []
        for window in roller.window_list:
            window_size = roller.window_width
            min_time = min(window.raw_data['Time'].unique())
            window_start_list.append(min_time)
            gold_standard = roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
            evaluator = Evaluator("../../"+gold_standard,sep="\t")
            sorted_edge_list = window.results_table
            sorted_edge_list.sort(['p_value'], ascending=[True], inplace=True)
            sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p_value'])]
            tpr,fpr, auroc = evaluator.calc_roc(sorted_edge_list)
            auroc_list.append(auroc[-1])
            
        pdb.set_trace()
        f = plt.figure()
        axarr1 = plt.subplot2grid((3,1),(0,0), rowspan=2)
        axarr2 = plt.subplot2grid((3,1),(2,0), sharex = axarr1)
        f.set_tight_layout(True)
        #figure, heatmap = Grapher.generate_heatmap_from_df(self.roller_results.raw_data)
        #axarr1 = figure
        time_vector = raw_data['Time'].unique()
        nrow,ncol=raw_data.shape
        nrepeats = nrow/len(time_vector)
        #create group assignment
        groups = [x for x in range(0,nrepeats) for i in range(0,len(time_vector))]
        raw_data['Group'] = groups
        sorted = raw_data.sort(['Time','Group'])
        heatmap_values = pd.DataFrame()
        col_labels = map(int,time_vector)
        partial_row_labels = raw_data.drop(['Time','Group'],1).columns.values.tolist()
        row_labels = [i for x in range(0,nrepeats) for i in partial_row_labels]
        for time in time_vector: 
            target_set = sorted[sorted['Time']==time].drop(['Time','Group'],1)
            vector = target_set.values.flatten()
            heatmap_values[time] = vector
        my_axis =axarr1.imshow(heatmap_values,interpolation='nearest',cmap=matplotlib.cm.RdBu)
        my_axi = my_axis.get_axes()
        axarr1.set_yticks(np.arange(heatmap_values.shape[0]))
        axarr1.yaxis.set_ticks_position('left')
        axarr1.set_yticklabels(row_labels)
        axarr1.set_xticks(np.arange(heatmap_values.shape[1]))
        axarr1.xaxis.set_ticks_position('top')
        xlabelsL = axarr1.set_xticklabels(col_labels)
        for label in xlabelsL:
            label.set_rotation(90)
            label.set_fontsize(8)
        for label in (axarr1.get_yticklabels()):
            label.set_fontsize(8)
        for l in axarr1.get_xticklines() + axarr1.get_yticklines():
            l.set_markersize(0)
        
        
        
        
        heatmapGS = matplotlib.gridspec.GridSpec(2,1,height_ratios=[2,0.25])
        axarr2.scatter(map(int,window_start_list), auroc_list)
        xlabels = axarr2.get_xticklabels()
        ylabels = axarr2.get_yticklabels()

        for label in xlabels:
            label.set_rotation(90)
            label.set_fontsize(8)
        for label in ylabels:
            label.set_fontsize(8)
        plt.tight_layout()
        #f.subplots_adjust(wspace=0.001,hspace=0.001,top=None,bottom=None)
        f.savefig("/home/jjw036/Roller/yeast_heatmap2.png")
if __name__ == '__main__':
    unittest.main()
