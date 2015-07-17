import sys
import Roller
import pandas as pd
import pdb

class Analyzer:
    
    """ Analyzer is a object that analyzes groups of Rollers. Figures out the completeness if your experiment. Checks for errored Rollers. Helps open and sort pickle files. It also has methods that ranks Rollers by their model and cragging scores."""

    def __init__(self, pickle_path_folder):
        self.sorted_edge_lists = []
        self.total_files_unpickled = []
        self.error_list = []

        for pickle_path in os.filedir(pickle_path_folder):
          self.current_pickle_path = pickle_path
          self.current_roller = pd.read_pickle(pickle_path)
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
                pdb.set_trace()
            except AttributeError:
                self.error_list.append(window_tag + "Window Index " + str(index) + " : No results table")
                
                

        
