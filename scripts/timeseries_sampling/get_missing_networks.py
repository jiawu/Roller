import pandas as pd
import pdb
import os
import time

def read_tdr_results(folder_list, folder_str):
    agg_df = pd.DataFrame()
    for input_folder in folder_list:
        for file_path in os.listdir(input_folder):
            if folder_str in file_path:
              try:
                  df = pd.read_csv(input_folder+file_path,sep='\t', engine='python')
              except pd.io.common.EmptyDataError:
                  continue
              agg_df = agg_df.append(df)
    return(agg_df)

output_path = "/home/jjw036/"

input_folder_list = ["/projects/p20519/roller_output/high_sampling/Lasso/"]  
#input_folder_list = ["/projects/p20519/roller_output/gnw/RandomForest/", "/projects/p20519/roller_output/gnw/Lasso/", "/projects/p20519/roller_output/gnw/Dionesus/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "sampling_comparison"
n_trials = 100

#datasets = ["_"]

#datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']
start = time.time()
agg_df = read_tdr_results(input_folder_list, folder_str = "2017-09")

result_counts = agg_df['file_path'].dropna().value_counts()
completed_files = result_counts[result_counts > 99].index.tolist()

job_file = pd.read_csv('/home/jjw036/Roller/pipelines/job_params_high_sampling.txt', sep = ' ', header=None)

job_file.columns = ['data_path','data_path', 'input_file', 'iterating_param', 'iterating_style']

mask = job_file['input_file'].str.contains('|'.join(completed_files))
new_jobs = job_file[~mask]
n_jobs = len(new_jobs)
new_jobs.to_csv('/home/jjw036/Roller/pipelines/job_params_high_sampling_missing.txt', sep = ' ', header=None, index=False)
print("Added {} missing jobs".format(n_jobs))
