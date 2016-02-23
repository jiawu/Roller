import pdb
import sys
import Pipelines as pl
import pandas as pd
from datetime import datetime
import numpy as np
import time

# saving the models for the iteration tests:
# to save the models for the iteration tests, we will save a dataframe (in the form of the final dataframe from Analyzer...) instead of a full model, because it is too computationally expensive, and as of this day, we are running out of room on QUEST.

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def main(data_folder, output_path, target_dataset, my_iterating_param, param_test_style, param_tests, n_trials):

    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    default_params = {'data_folder':data_folder, 'file_path':target_dataset, 'td_window':10,'min_lag':1,'max_lag':3,'n_trees':500,'permutation_n':5, 'lag_method':'mean_mean', 'calc_mse':False, 'bootstrap_n':50,'n_trials':n_trials, 'run_time':current_time, 'sort_by': 'rank','iterating_param':my_iterating_param, 'filter_noisy':False}

    overall_df = pd.DataFrame()

    #**kwargs allows me to change the iterating parameter very easily

    trial_times = []
    for current_param_value in param_tests:
        for trial in range(0,n_trials):
            trial_start = time.time()
            run_params = default_params.copy()
            run_params[my_iterating_param]=current_param_value

            if (run_params['td_window'] == 21) or (('_sampling' in run_params['file_path']) and run_params['td_window'] == 7) or (('dream8/insilico' in run_params['file_path']) and run_params['td_window'] == 11) or(('dream8/invitro' in run_params['file_path']) and run_params['td_window'] == 7):
                run_params['min_lag'] = 0
                run_params['max_lag'] = 0
            
            if 'community' in data_folder:
                roc,pr, tdr = pl.get_td_community(**run_params)
            else:
                roc,pr, tdr = pl.get_td_stats(**run_params)
            run_params['auroc']=roc
            run_params['aupr']=pr

            trial_end = time.time()
            run_params['trial_time'] = trial_end-trial_start
            run_result=pd.Series(run_params)
            overall_df = overall_df.append(run_result, ignore_index=True)
            print(run_result)
    overall_df.to_csv(output_path+current_time+'.tsv', index=False, sep='\t')

if __name__ == "__main__":
    """
    Runs a number of trials with tdRoller with specified params, saves a dataframe with the run params and result

    Example call: python run_tdSwing_scan.py data_folder output_path target_dataset iterating_param param_test_style
    :param data_folder: /projects/p20519/roller_output/optimizing_window_size/RandomForest/janes
    :param output_path: /projects/p20519/roller_output/stability_analysis/RandomForest/janes_ntrees_
        the output will be a tsv named janes_ntrees_<currenttime>.tsv
    :param target_dataset: /projects/p20519/Swing/data/invitro/janes_timeseries.tsv or/projects/p20519/Swing/data/dream4/ecoli_timeseries.tsv
    :param my_iterating_param: = n_trees

    :param param_test_style: str that defines either logarithmic 10,100,1000 or specify a min/max or string
    

    """

    data_folder = str(sys.argv[1])
    output_path = str(sys.argv[2])
    target_dataset = str(sys.argv[3])
    my_iterating_param = str(sys.argv[4])
    param_test_style = str(sys.argv[5])

    if param_test_style == "log":
        param_tests = [10,100,500,1000]

    elif param_test_style == "minmax":
        param_min = int(sys.argv[6])
        param_max = int(sys.argv[7])
        param_tests = [i for i in range(param_min, param_max+1)]
    
    elif param_test_style == "num":
        param_tests = [int(x) for x in sys.argv[6:]]
        
    elif param_test_style == "string":
        param_tests = [str(x) for x in sys.argv[6:]]
    
    elif param_test_style == "boolean":
        param_tests = [str2bool(x) for x in sys.argv[6:]]
        
    n_trials = 2

    #always save the full parameter list and date in the dataframe for each test. for posterity!

    main(data_folder, output_path, target_dataset, my_iterating_param, param_test_style, param_tests, n_trials)


    
