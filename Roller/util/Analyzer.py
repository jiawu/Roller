import sys, os
import Roller
import pandas as pd
import numpy as np
import pdb
from Roller.util.Evaluator import Evaluator
import warnings


class Analyzer:
    
    """ Analyzer is a object that analyzes groups of Rollers. Figures out the completeness if your experiment. Checks for errored Rollers. Helps open and sort pickle files. It also has methods that ranks Rollers by their model and cragging scores."""

    def __init__(self, pickle_path_folder):
        self.sorted_edge_lists = []
        self.total_files_unpickled = []
        self.error_list = []
        self.overall_df = pd.DataFrame()
        for pickle_path in os.listdir(pickle_path_folder):
            self.current_pickle_path = pickle_path
            try:
                
                self.current_roller = pd.read_pickle(pickle_path_folder+"/"+ pickle_path)
                df = self.aggregate_ranked_list(self.current_roller)
                self.overall_df = self.overall_df.append(df)

            except KeyError:
                continue

            self.total_files_unpickled.append(pickle_path)

        

        

        ### correlation between ###

        ### guess best window ###

        ### check if best window is better than status quo ###
    
    def get_best_window(self):
        best_row = self.overall_df.loc[self.overall_df['window_width'].idxmax()]
        return(best_row)

    def get_max_window(self):
        ### identify status quo ###
        max_row = self.overall_df.loc[self.overall_df['window_width'].idxmax()]
        max_width = self.current_roller.overall_width
        if max_row['window_width'] != max_width:
            max_width = max_row['window_width']
            #find max window size
            warnings.warn("Roller with all timepoints is not present. Using Roller with a maximum width of %s as comparison window" % (max_width))
        return(max_row)


    def get_window_tag(self):
        window_size = self.current_roller.window_width
        tag = self.current_pickle_path + "Width: " + str(window_size)
        return(tag)
    
    def aggregate_ranked_list(self,roller_obj):
        #generate a dataframe that aggregates the window stats for each window/roller
        df = pd.DataFrame()
        pickle_paths = []
        network_paths = []
        auroc_list = []
        aupr_list = []
        window_index_list = []
        crag_mse_average_list = []
        crag_r2_average_list = []
        crag_ev_average_list =[]
        
        crag_mse_median_list = []
        crag_r2_median_list = []
        crag_ev_median_list =[]
        
        crag_mse_max_list = []
        crag_r2_max_list = []
        crag_ev_max_list = []

        window_width_list = []

        for index,window in enumerate(roller_obj.window_list):
            pickle_paths.append(self.current_pickle_path)
            network_paths.append(self.current_roller.file_path)
            window_width_list.append(roller_obj.window_width)

    def check_ranked_list(self,roller_obj):
        for index,window in enumerate(roller_obj.window_list):
            try:
                sorted_edge_list = window.results_table
                #check if the sorted edge list actually has importance/ranking values. if it doesn't, raise an error
                if len(sorted_edge_list.columns) < 2:
                    raise AttributeError
                gold_standard = self.current_roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
                evaluator = Evaluator("/projects/p20519/Roller/" +gold_standard,sep="\t")
                sorted_edge_list.sort(['importance'], ascending=[False], inplace=True)
                tpr,fpr,auroc = evaluator.calc_roc(sorted_edge_list)
                precision,recall,aupr = evaluator.calc_pr(sorted_edge_list)
                print(aupr.values[-1])
                print(auroc.values[-1])
                
                auroc_list.append(auroc.values[-1])
                aupr_list.append(aupr.values[-1])

                model_crag=[{  'ev': 0,
                              'mse':0,
                              'r2':0
                          }]
                if roller_obj.window_width != 34:
                    crag_iterations = len(window.test_scores)/window.n_genes
                    cragging_scores = []
                    for i in range(0,crag_iterations):
                        cragging_scores.append(window.test_scores[i*window.n_genes:(i+1)*window.n_genes])
                    # unfortunately, get_coeffs is also called by the null model, so the cragging function also evaluates null models and appends them to window.training_scores. The first indices are the cragging scores for the model.

                    model_crag = cragging_scores[0]
                crag_ev_average_list.append(self.average_dict(model_crag,'ev'))
                crag_mse_average_list.append(self.average_dict(model_crag,'mse'))
                crag_r2_average_list.append(self.average_dict(model_crag,'r2'))
                
                crag_ev_median_list.append(self.median_dict(model_crag,'ev'))
                crag_mse_median_list.append(self.median_dict(model_crag,'mse'))
                crag_r2_median_list.append(self.median_dict(model_crag,'r2'))
                
                crag_ev_max_list.append(self.max_dict(model_crag,'ev'))
                crag_mse_max_list.append(self.max_dict(model_crag,'mse'))
                crag_r2_max_list.append(self.max_dict(model_crag,'r2'))

                window_index_list.append(index+1)



            except (AttributeError,IndexError):
                window_tag = self.get_window_tag()
                self.error_list.append(window_tag + "Window Index " + str(index) + " : No results table")
        
        if auroc_list:
          df = pd.DataFrame( {'pickle_paths':pickle_paths,'network_paths':network_paths,'auroc':auroc_list,'aupr':aupr_list,'window_index':window_index_list,'crag_mse_average':crag_mse_average_list,'crag_ev_average':crag_ev_average_list,'crag_r2_average':crag_r2_average_list,'crag_mse_median':crag_mse_median_list,'crag_ev_median':crag_ev_median_list,'crag_r2_median':crag_r2_median_list,'crag_ev_max':crag_ev_max_list,'crag_mse_max':crag_mse_max_list,'crag_r2_max':crag_r2_max_list,'window_width':window_width_list})
        return(df) 



    def average_dict(self,total,key): 
        return( (sum(d[key] for d in total))/len(total))

    def median_dict(self,total,key):
        aggr = [x[key] for x in total]
        return(np.median(aggr))

    def max_dict(self,total,key):
        aggr=[x[key] for x in total]
        return(np.max(aggr))
