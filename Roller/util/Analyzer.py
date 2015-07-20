import sys, os
import Roller
import pandas as pd
import pdb
from Roller.util.Evaluator import Evaluator


class Analyzer:
    
    """ Analyzer is a object that analyzes groups of Rollers. Figures out the completeness if your experiment. Checks for errored Rollers. Helps open and sort pickle files. It also has methods that ranks Rollers by their model and cragging scores."""

    def __init__(self, pickle_path_folder):
        self.sorted_edge_lists = []
        self.total_files_unpickled = []
        self.error_list = []

        for pickle_path in os.listdir(pickle_path_folder):
          self.current_pickle_path = pickle_path
          self.current_roller = pd.read_pickle(pickle_path_folder+"/"+ pickle_path)
          self.is_ranked = self.check_ranked_list(self.current_roller)

          self.total_files_unpickled.append(pickle_path)


    def get_window_tag(self):
        window_size = self.current_roller.window_width
        tag = self.current_pickle_path + "Width: " + str(window_size)
        return(tag)
    
    def check_ranked_list(self,roller_obj):
        for index,window in enumerate(roller_obj.window_list):
            try:
                sorted_edge_list = window.results_table
                #check if the sorted edge list actually has importance/ranking values. if it doesn't, raise an error
                if len(sorted_edge_list.columns) < 2:
                    raise AttributeError
                pdb.set_trace()
                gold_standard = self.current_roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
                evaluator = Evaluator("/projects/p20519/Roller/" +gold_standard,sep="\t")
                sorted_edge_list.sort(['importance'], ascending=[False], inplace=True)
                tpr,fpr, auroc = evaluator.calc_roc(sorted_edge_list)
            except AttributeError:
                window_tag = self.get_window_tag()
                self.error_list.append(window_tag + "Window Index " + str(index) + " : No results table")
                
                

        
