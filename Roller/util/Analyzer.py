import sys, os
import Roller
import pandas as pd
import numpy as np
import pdb
from Roller.util.Evaluator import Evaluator
import warnings
import Roller.util.utility_module as Rutil


class Analyzer:
    
    """ Analyzer is a object that analyzes groups of Rollers. Figures out the completeness if your experiment. Checks for errored Rollers. Helps open and sort pickle files. It also has methods that ranks Rollers by their model and cragging scores."""

    def __init__(self, pickle_path_folder):
        self.sorted_edge_lists = []
        self.total_files_unpickled = []
        self.error_list = []
        self.overall_df = pd.DataFrame()
        self.pickle_paths = os.listdir(pickle_path_folder)
        self.pickle_folder = pickle_path_folder
        self.max_width_results = None
        self.min_width_results = None
        self.all_ranked_lists = []
        

        for pickle_path in self.pickle_paths:
            self.current_pickle_path = pickle_path
            try:
                self.current_roller = pd.read_pickle(pickle_path_folder+"/"+ pickle_path)
                # saving rollers

                df = self.aggregate_ranked_list(self.current_roller)
                self.overall_df = self.overall_df.append(df)
            except KeyError:
                continue
            self.total_files_unpickled.append(pickle_path)
    
    def get_ranked_list(self, pickle_path, target_window_index):
        roller_obj = pd.read_pickle(pickle_path)
        for index,window in enumerate(roller_obj.window_list):
            if index == target_window_index:
                if not roller_obj.file_path.startswith("/"):
                    roller_obj.file_path = roller_obj.file_path.replace("data/","/projects/p20519/Roller/data/")
                    sorted_edge_list = window.results_table 
                    sorted_edge_list.sort(['importance'], ascending=[False], inplace=True)
        return(sorted_edge_list)
        

    def aggregate_best_windows(self, top_percentile=10):
        corr_list = []
        width_list = []
        agg_auroc = []
        agg_auroc_table = []
        for target_width in range(4,22):
            width_list.append(target_width)
            target_df = self.overall_df[self.overall_df['window_width']==target_width]
            corr_list.append(target_df.corr()['auroc'])

            sorted = self.overall_df.sort('crag_mse_average',ascending=True)
            
            target_sorted = sorted[(sorted['window_width']<target_width) &(sorted['window_width']>target_width-3)]
            n_rows = round(top_percentile*.01*len(target_sorted))
            print(n_rows," windows incorporated")
            top_windows = target_sorted.head(n_rows)
            
            ## aggregate ranked lists from a list of pickle paths and indices
            top_ranked_lists = []
            for index, row in top_windows.iterrows():
                full_path = self.pickle_folder + row['pickle_paths']
                window_ranked_list=self.get_ranked_list(full_path, row['window_index'])
                top_ranked_lists.append(window_ranked_list)
            agg_auroc_table.append(top_ranked_lists)
            ## change the value into a rank
            for ranked_list in top_ranked_lists:
                ranked_list['importance_rank'] = ranked_list.rank(axis=0, ascending=False)['importance']
            
            
            
            top_auroc=self.overall_df.sort('auroc', ascending=False).head()
            top_ranked_lists2 = []
            for index, row in top_auroc.iterrows():
                full_path = self.pickle_folder + row['pickle_paths']
                window_ranked_list=self.get_ranked_list(full_path, row['window_index'])
                top_ranked_lists2.append(window_ranked_list)
            
            ## change the value into a rank
            for ranked_list in top_ranked_lists2:
                ranked_list['importance_rank'] = ranked_list.rank(axis=0, ascending=False)['importance']
            
            
            
            averaged_lists = Rutil.average_rank(top_ranked_lists,col_string='importance_rank')
            gold_standard = self.current_roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
            averaged_lists.sort('mean-rank',ascending=True,inplace=True)
            evaluator = Evaluator(gold_standard,sep="\t")
            tpr,fpr,auroc=evaluator.calc_roc(averaged_lists)
            precision, recall, aupr = evaluator.calc_pr(averaged_lists)
            agg_auroc.append(auroc.tolist()[-1])
        
        my_r = zip(width_list, agg_auroc)
        pdb.set_trace()

    def predict_best_window(self):
        #max_value = self.overall_df['crag_mse_average'].max()
        max_value = 0
        counter = 1
        while max_value == 0:
            max_value = self.overall_df['crag_mse_average'].nsmallest(counter).values[-1]
            counter += 1

        best_row = self.overall_df[self.overall_df['crag_mse_average']== max_value]
        return(best_row)

    def get_correlation(self):
        return(self.overall_df.corr())

        
    def load_list(self,csv_file_path):
        self.overall_df = pd.read_csv(csv_file_path)
        return(df)

    def get_best_window(self):
        best_row = self.overall_df.loc[self.overall_df['window_width'].idxmax()]
        return(best_row)

    def get_max_window(self):
        ### identify status quo ###
        max_row = self.overall_df.loc[self.overall_df['window_width'].idxmax()]
        max_width = self.current_roller.overall_width
        """
        if max_row['window_width'] != max_width:
            max_width = max_row['window_width']
            #find max window size
            warnings.warn("Roller with all timepoints is not present. Using Roller with a maximum width of %s as comparison window" % (max_width))
        """
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
            if not self.current_roller.file_path.startswith("/"):
                self.current_roller.file_path = self.current_roller.file_path.replace("data/","/projects/p20519/Roller/data/")
            pickle_paths.append(self.current_pickle_path)
            network_paths.append(self.current_roller.file_path)
            window_width_list.append(roller_obj.window_width)

            try:
                sorted_edge_list = window.results_table
                #check if the sorted edge list actually has importance/ranking values. if it doesn't, raise an error
                if len(sorted_edge_list.columns) < 2:
                    raise AttributeError
                ## replace relative file paths with absolute file paths

                gold_standard = self.current_roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
               
                evaluator = Evaluator(gold_standard,sep="\t")
                #sorted_edge_list.sort(['p_value'], ascending=[True], inplace=True)
                sorted_edge_list.sort(['importance'], ascending=[False], inplace=True)
                self.all_ranked_lists.append(sorted_edge_list)
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
                if roller_obj.window_width != roller_obj.overall_width:
                    if self.max_width_results:
                        if auroc.values[-1] > self.max_width_results['auroc'].values[-1]:
                            self.min_width_results = {'tpr':tpr, 'fpr':fpr,'auroc':auroc,'precision':precision,'recall':recall,'aupr':aupr}

                    crag_iterations = len(window.test_scores)/window.n_genes
                    cragging_scores = []
                    for i in range(0,crag_iterations):
                        cragging_scores.append(window.test_scores[i*window.n_genes:(i+1)*window.n_genes])
                    # unfortunately, get_coeffs is also called by the null model, so the cragging function also evaluates null models and appends them to window.training_scores. The first indices are the cragging scores for the model.

                    model_crag = cragging_scores[0]
                else:
                    self.max_width_results = {'tpr':tpr, 'fpr':fpr,'auroc':auroc,'precision':precision,'recall':recall,'aupr':aupr}
                    
                crag_ev_average_list.append(self.average_dict(model_crag,'ev'))
                crag_mse_average_list.append(self.average_dict(model_crag,'mse'))
                crag_r2_average_list.append(self.average_dict(model_crag,'r2'))
                
                crag_ev_median_list.append(self.median_dict(model_crag,'ev'))
                crag_mse_median_list.append(self.median_dict(model_crag,'mse'))
                crag_r2_median_list.append(self.median_dict(model_crag,'r2'))
                
                crag_ev_max_list.append(self.max_dict(model_crag,'ev'))
                crag_mse_max_list.append(self.max_dict(model_crag,'mse'))
                crag_r2_max_list.append(self.max_dict(model_crag,'r2'))

                window_index_list.append(index)

            except (AttributeError,IndexError):
                window_tag = self.get_window_tag()
                self.error_list.append(window_tag + "Window Index " + str(index) + " : No results table")
        
        if auroc_list:
          if roller_obj.window_width == roller_obj.overall_width:
              window_index_list = [0]
              crag_mse_average_list = [0]
              crag_r2_average_list = [0]
              crag_ev_average_list =[0]
              
              crag_mse_median_list = [0]
              crag_r2_median_list = [0]
              crag_ev_median_list =[0]
              
              crag_mse_max_list = [0]
              crag_r2_max_list = [0]
              crag_ev_max_list = [0]
          df = pd.DataFrame( {'pickle_paths':pickle_paths,
                              'network_paths':network_paths,
                              'auroc':auroc_list,
                              'aupr':aupr_list,
                              'window_index':window_index_list,
                              'crag_mse_average':crag_mse_average_list,
                              'crag_ev_average':crag_ev_average_list,
                              'crag_r2_average':crag_r2_average_list,
                              'crag_mse_median':crag_mse_median_list,
                              'crag_ev_median':crag_ev_median_list,
                              'crag_r2_median':crag_r2_median_list,
                              'crag_ev_max':crag_ev_max_list,
                              'crag_mse_max':crag_mse_max_list,
                              'crag_r2_max':crag_r2_max_list,
                              'window_width':window_width_list})
        return(df) 



    def average_dict(self,total,key): 
        return( (sum(d[key] for d in total))/len(total))

    def median_dict(self,total,key):
        aggr = [x[key] for x in total]
        return(np.median(aggr))

    def max_dict(self,total,key):
        aggr=[x[key] for x in total]
        return(np.max(aggr))

