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


def parse_job_index(job_index):
    ## Job Indices:
    ## 1-20: ecoli 1 to 20, 10 node networks, Dionesus
    ## 21-40: yeast 1 to 20, 10 node networks, Dionesus
    ## 41-60: ecoli 1 to 20, 100 node networks, Dionesus
    # 61-80: yeast 1 to 20, 100 node networks, Dionesus 
    ## 81 omranian Dionesus

    ## 101-120: ecoli 1 to 20, 10 node networks, RF
    ## 121-140: yeast 1 to 20, 10 node networks, RF
    ## 141-160: ecoli 1 to 20, 100 node networks, RF
    ## 161-180: yeast 1 to 20, 100 node networks, RF
    ## 181 omranian RF

    ## 201-220: ecoli 1 to 20, 10 node networks, LASSO
    ## 221-240: yeast 1 to 20, 10 node networks, LASSO          
    ## 241-260: ecoli 1 to 20, 100 node networks, LASSO
    ## 261-280: yeast 1 to 20, 100 node networks, LASSO
    ## 281 omranian lasso

    ## 301-320: ecoli 1 to 20, 10 node networks, community
    ## 321-340: yeast 1 to 20, 10 node networks, community
    ## 341-360: ecoli 1 to 20, 100 node networks, community
    ## 361-380: yeast 1 to 20, 100 node networks, community
    ## 381 omranian community

    if job_index < 100:
        inf_method = 'Dionesus'
    elif job_index < 200:
        inf_method = 'RandomForest'
    elif job_index < 300:
        inf_method = 'Lasso'
    elif job_index < 400:
        inf_method = 'community'

    organism_index = job_index%100
    
    if organism_index < 21:
        organism = 'Ecoli'
    elif organism_index < 41:
        organism = 'Yeast'
    elif organism_index < 61:
        organism = 'Ecoli100'
    elif organism_index < 81:
        organism = 'Yeast100'
    elif organism_index == 81:
        organism = 'omranian'
    elif organism_index >= 82:
        organism = 'dream5'

    network_index = organism_index%20
    if network_index == 0:
        network_index = 20
    
    if organism is 'omranian':
        file_path = "/home/jjw036/Roller/data/invitro/omranian_parsed_timeseries.tsv"
        data_folder = "/projects/p20519/roller_output/ranks/%s/%s_" % (inf_method,organism)

    elif organism is 'dream5':
        file_path = "/home/jjw036/Roller/data/dream5/insilico_timeseries.tsv"
        data_folder = "/projects/p20519/roller_ouput/ranks/%s/%s_" % (inf_method,organism)

    else:
        file_path = "/home/jjw036/Roller/data/gnw_insilico/network_data/%s/%s-%d_timeseries.tsv" % (organism,organism,network_index)
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
        # 11 methods of interest. so it goes from 82 to 93
        max_window = 21
        if 'omranian' in file_path:
            methods_of_interest = ['-td_6','-td_5','-td_4', '-ml_01', '-ml_11', '-ml_02','-ml_12', '-ml_22', '-ml_23', '-ml_13']
            max_window = 6
        elif 'dream5' in file_path:
            methods_of_interest = methods_of_interest[(job-82)%100]
        
        method_strings = [inf_method + x for x in methods_of_interest]
        
        for method_string in method_strings:
            #update current_time
            current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
            run_params['run_time'] = current_time
            
            trial_start = time.time()
            
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
