import os
import pandas as pd
import pdb
from Swing import Swing
import pickle

gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None

agg_df = pd.DataFrame()


input_folder_list = ["/home/jjw036/Roller/pipelines/lasso_alpha_scans/output/"]
for input_folder in input_folder_list:
    for file_path in os.listdir(input_folder):
        df = pd.read_csv(input_folder+file_path,sep=',|\t', engine='python')
        agg_df = agg_df.append(df)

my_group = agg_df.groupby(['file_path', 'min_lag', 'max_lag', 'td_window'])
groups = []
group_names = []
group_alphas = []
for name, group in my_group:
    group_names.append(name)
    groups.append(group)

for param in group_names:
    file_path = param[0]
    min_lag = int(param[1])
    max_lag = int(param[2])
    td_window = int(param[3])
    tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag = min_lag,max_lag=max_lag,window_type='Lasso')
    tdr.zscore_all_data()
    tdr.set_window(td_window)
    tdr.create_windows()
    tdr.optimize_params()
    selected_alphas = [x.alpha for x in tdr.window_list]
    group_alphas.append(selected_alphas)
pickle_name = 'selected_alphas.pkl'
pickle.dump(group_alphas, open(pickle_name,'wb'))
