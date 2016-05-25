import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/jjw036/Roller/pipelines')
import Pipelines as pl
import numpy as np
import pandas as pd
from datetime import datetime
import time

from Swing.util.Evaluator import Evaluator
import pickle
import pdb

""" 
Sample calls:
- to run low sampling dream4 data with n_trials:
  submit script:
    msub -t dream4.[1-10] submit_analyze_promotion_plots_dream4.py
  call inside submit:
    python analyze_promotion_plot_dream4.py ${param_set} 100

- to run high sampling data:
  submit script:
    msub -t dream4.[11-55] submit_analyze_promotion_plots_dream4.py
  call inside submit:
    python analyze_promotion_plot_dream4.py ${param_set} 1
    

"""
def parse_method(method_string, max_window):
    """
    Parameters:
    method_string is a str that contains abbreviated information about parameter settings
    ie Dionesus-td10 indicates that dionesus will be run with a windowsize of 10

    """
    min_lag = 1
    max_lag = 3
    td_window = 15

    inf_method = method_string.split('-')[0]

    misc = method_string.split('-')[1]

    if "td" in misc:
        td_window = int(misc.split('_')[1])
        if td_window == max_window:
            min_lag = 0
            max_lag = 0
        elif td_window + max_lag > max_window:
            max_lag = max_window - td_window

    elif "ml" in misc:
        case = str(misc.split('_')[1])
        min_lag = int(case[0])
        max_lag = int(case[1])
        if td_window > max_window:
            td_window = 3
    return(inf_method, td_window, min_lag, max_lag)

def parse_method_hs(method_string, max_window):
    """
    Parameters:
    method_string is a str that contains abbreviated information about parameter settings
    ie Dionesus-td10 indicates that dionesus will be run with a windowsize of 10

    """
    min_lag = 25
    max_lag = 75
    td_window = 375

    inf_method = method_string.split('-')[0]

    misc = method_string.split('-')[1]

    if "td" in misc:
        td_window = int(misc.split('_')[1])
        if td_window == max_window:
            min_lag = 0
            max_lag = 0
        elif td_window + max_lag > max_window:
            max_lag = max_window - td_window

    elif "ml" in misc:
        case = str(misc.split('_')[1])
        min_lag = int(case[0:2])
        max_lag = int(case[2:4])
        if td_window > max_window:
            td_window = 375
    return(inf_method, td_window, min_lag, max_lag)

def parse_job_index(job_index):
    """ Job Indices:
    1-5: size 10
    5-10: size 100
    11-15: high sampling size 11 - first method
    x # of methods
    11-55

    """
    if job_index < 100:
        inf_method = 'Dionesus'
    elif job_index < 200:
        inf_method = 'RandomForest'
    elif job_index < 300:
        inf_method = 'Lasso'
    elif job_index < 400:
        inf_method = 'community'

    organism_index = job_index%100
    
    if organism_index < 6:
        organism = 'dream4_insilico_size_10'
    elif organism_index < 11:
        organism = 'dream4_insilico_size_100'
    else:
        organism = 'dream4_insilico_size_10_high_sampling'

    network_index = organism_index%5
    if network_index == 0:
        network_index = 5
    
    if organism is 'dream4_insilico_size_10_high_sampling':
        file_path = "/home/jjw036/Roller/data/dream4/high_sampling/insilico_size10_%d_timeseries.tsv" % (network_index)
        data_folder = "/projects/p20519/roller_output/ranks/%s/%s_%d_" % (inf_method,organism,network_index)
    elif 'size_100' in organism:
        file_path = "/home/jjw036/Roller/data/dream4/insilico_size100_%d_timeseries.tsv" % (network_index)
        data_folder = "/projects/p20519/roller_output/ranks/%s/%s_%d_" % (inf_method,organism,network_index)
    else:
        file_path = "/home/jjw036/Roller/data/dream4/insilico_size10_%d_timeseries.tsv" % (network_index)
        data_folder = "/projects/p20519/roller_output/ranks/%s/%s_%d_" % (inf_method,organism,network_index)

    return(file_path, data_folder, inf_method)
    
def main(job,n_trials=1):
    """
    Prints a series of text files with the parameters: denoted by ! and the resulting ranked list for each trial, denoted by a time-stamp
    Parameters:
        job - an int that corresponds to a method/network combination
            ** each job runs a series of methods for one network, each method produces its own rank file
    """
    for i in range(n_trials):
        current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

        default_params = {'data_folder':None, 'file_path':None, 'td_window':0,'min_lag':0,'max_lag':0,'n_trees':100,'permutation_n':5, 'lag_method':'mean_mean', 'calc_mse':False, 'bootstrap_n':50,'n_trials':n_trials, 'run_time':current_time, 'sort_by': 'rank','iterating_param':'promotion', 'filter_noisy':False, 'alpha': None, 'trial_time': 0, 'auroc': 0.0, 'aupr': 0.0}

        run_params = default_params.copy()
        file_path,data_folder, inf_method = parse_job_index(job)

        run_params['file_path'] = file_path
        run_params['data_folder'] = data_folder

        methods_of_interest = ['-td_10', '-td_15','-td_21', '-ml_01', '-ml_11', '-ml_02','-ml_03','-ml_12', '-ml_22', '-ml_33', '-ml_34']
        max_window = 21
        if 'high_sampling' in file_path:
            methods_of_interest = ['-td_501','-td_375','-td_400', '-ml_0510', '-ml_2550', '-ml_2575','-ml_2590', '-ml_9090']
            ### High sampling data takes a long time to process, so we split up the method based on the job index that is passed in. Therefore, even though there are only 5 networks, the job indices to run high sampling data is from 10 to 50
            method_indices = [x for x in range(len(methods_of_interest)) for y in range(5)]
            
            method_index = method_indices[(job-11)%100]
            methods_of_interest = [methods_of_interest[method_index]]
            
            ## 11 is net1 method 1 = 0
            ## 15 is net5 method 1 
            ## 16 is net1 method 2 = 1
            ## 21 is net1 method 3 = 2

            max_window = 501
        method_strings = [inf_method + x for x in methods_of_interest]
        
        for method_string in method_strings:
            #update current_time
            current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
            run_params['run_time'] = current_time
            
            trial_start = time.time()
        
            if 'high_sampling' in file_path:
                inf_method, td_window, min_lag, max_lag = parse_method_hs(method_string, max_window)
            else:
                inf_method, td_window, min_lag, max_lag = parse_method(method_string, max_window)
            run_params['td_window'] = td_window
            run_params['min_lag'] = min_lag
            run_params['max_lag'] = max_lag
            print(run_params)
            
            if 'community' in data_folder:
                roc,pr, tdr, rank_table = pl.get_td_community(**run_params)
            else:    
                roc, pr, tdr = pl.get_td_stats_custom(**run_params)
            
            run_params['auroc']=roc
            run_params['aupr']=pr

            trial_end = time.time()
            run_params['trial_time'] = trial_end-trial_start

            if 'community' in data_folder:
                result_table = rank_table
            else:
                result_table = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
                result_table['rank_importance'] = np.arange(len(result_table))
            
            output_file = data_folder+current_time+".csv"
            
            with open(output_file,'a') as output:
                for key, value in run_params.items():
                    output.write('!%s,%s\n' % (key, value))
                result_table.to_csv(output, header=True, index=False, sep='\t')

if __name__ == '__main__':
    job_index = int(sys.argv[1])    
    if len(sys.argv) >= 3:
        n_trials = int(sys.argv[2])
    else:
        n_trials = 1
    main(job_index, n_trials)
