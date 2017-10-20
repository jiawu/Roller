import pdb
import sys
import Pipelines as pl
import pandas as pd
from datetime import datetime
import numpy as np
import time
import tempfile
import os

# saving the models for the iteration tests:
# to save the models for the iteration tests, we will save a dataframe (in the form of the final dataframe from Analyzer...) instead of a full model, because it is too computationally expensive, and as of this day, we are running out of room on QUEST.

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def main(data_folder, output_path, target_dataset, my_iterating_param, param_test_style, param_tests, n_trials):

    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S:%f')
    if 'Dionesus' in data_folder:
        n_trials = 2
    default_params = {'data_folder':data_folder, 'file_path':target_dataset, 'td_window':15,'min_lag':1,'max_lag':3,'n_trees':500,'permutation_n':5, 'lag_method':'mean_mean', 'calc_mse':False, 'bootstrap_n':5,'n_trials':n_trials, 'run_time':current_time, 'sort_by': 'rank','iterating_param':my_iterating_param, 'filter_noisy':False, 'alpha': None}

    overall_df = pd.DataFrame()

    #**kwargs allows me to change the iterating parameter very easily

    trial_times = []
    for current_param_value in param_tests:
        for trial in range(0,n_trials):
            trial_start = time.time()
            run_params = default_params.copy()

            if param_test_style == "pair":
                run_params[my_iterating_param[0]] = current_param_value[0]
                run_params[my_iterating_param[1]] = current_param_value[1]
            
            elif param_test_style == "triplet":
                run_params[my_iterating_param[0]] = current_param_value[0]
                run_params[my_iterating_param[1]] = current_param_value[1]
                run_params[my_iterating_param[2]] = current_param_value[2]
            
            else:
                run_params[my_iterating_param]=current_param_value
            
            # Check max_lag restriction
            if 'size10' in run_params['data_folder']:
                max_window = 21
            if 'high_sampling' in run_params['data_folder']:
                if 'even' in run_params['file_path']:
                    max_window = 7
                else:
                    interval = run_params['file_path'].split('/')[-1].split('_')[1]
                    max_window = int(1000/int(interval)+1)
            elif 'gardner_out' in run_params['data_folder']:
                interval = run_params['file_path'].split('/')[-1].split('_')[1]
                max_window = int(round(14/int(interval)))
            else:
                max_window = 21
            lag_gap = max_window-run_params['td_window']
            # max window = 1, td window = 1
            # lag gap = 1-1 =0 
            # if lag gap (0) <= max_lag = 1
            if lag_gap <= run_params['max_lag']:
                run_params['max_lag'] = lag_gap
            if run_params['max_lag'] >= max_window:
                run_params['max_lag'] = max_window - 1
            
            if run_params['min_lag'] > run_params['max_lag']:
                run_params['min_lag'] = run_params['max_lag']
            
            if 'community' in data_folder:
                roc,pr, tdr, _ = pl.get_td_community(**run_params)
            else:
                roc,pr, tdr, _ = pl.get_td_stats(**run_params)
            run_params['auroc']=roc
            run_params['aupr']=pr

            trial_end = time.time()
            run_params['trial_time'] = trial_end-trial_start
            for key in run_params.keys():
                if run_params[key] is None:
                    run_params[key] = "None"
            run_result=pd.Series(run_params)
            overall_df = overall_df.append(run_result, ignore_index=True)
            print(run_result)
    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S:%f')
    full_path = output_path+current_time
    directory = os.path.dirname(full_path)
    _, filename = os.path.split(full_path)
    with tempfile.NamedTemporaryFile(prefix=filename, suffix='.tsv', dir=directory, delete=False) as temp:
        overall_df.to_csv(temp.name, index=False, sep='\t')
        print(temp.name)
        temp.close()

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
    param_test_combo = str(sys.argv[5])

    param_test_style = param_test_combo.split("_")[0]

    if param_test_style == "log":
        param_tests = [10,100,500,1000]

    elif param_test_style == "minmax":
        param_min = param_test_combo[1]
        param_max = param_test_combo[2]
        param_tests = [i for i in range(param_min, param_max+1)]
    
    elif param_test_style == "num":
        param_tests = [int(x) for x in param_test_combo.split("_")[1:]]
        
    elif param_test_style == "string":
        param_tests = [str(x) for x in param_test_combo.split("_")[1:]]

    elif param_test_style == "boolean":
        param_tests = [False, True]

    elif param_test_style == "pair":
        pli =param_test_combo.split("_")
        param_tests = list(zip( map(int, pli[1::2]), map(int, pli[2::2])))
        my_iterating_param = my_iterating_param.split("^")
    
    elif param_test_style == "triplet":
        pli =param_test_combo.split("_")
        param_tests = list(zip( map(int, pli[1::3]), map(int, pli[2::3]), map(int, pli[3::3]) ) )
        my_iterating_param = my_iterating_param.split("^")
        
    n_trials = 50

    #always save the full parameter list and date in the dataframe for each test. for posterity!

    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S:%f')

    overall_df = pd.DataFrame()

    #**kwargs allows me to change the iterating parameter very easily


    #always save the full parameter list and date in the dataframe for each test. for posterity!

    main(data_folder, output_path, target_dataset, my_iterating_param, param_test_style, param_tests, n_trials)
    
