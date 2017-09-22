import pandas as pd
import pdb

RF_sampling = pd.read_pickle('D_high_sampling_networks.pkl')

pdb.set_trace()
completed_files = RF_sampling['file_path'].unique().tolist()

job_file = pd.read_csv('/home/jjw036/Roller/pipelines/job_params_high_sampling.txt', sep = ' ', header=None)

job_file.columns = ['data_path','data_path', 'input_file', 'iterating_param', 'iterating_style']

mask = job_file['input_file'].str.contains('|'.join(completed_files))

new_jobs = job_file[~mask]

new_jobs.to_csv('/home/jjw036/Roller/pipelines/job_params_high_sampling_missing.txt', sep = ' ', header=None, index=False)
print("Added missing jobs")
