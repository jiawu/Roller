import os
import pandas as pd
import pdb
from parse_cluster_info import main as pinf
def parse_summary(fp, lag_thresh=4, fill_na=False, median_thresh=2):
    df = pd.read_csv(fp, sep='\t')
    aupr_cols = [col for col in df.columns if 'swing_aupr' in col]
    auroc_cols = [col for col in df.columns if 'swing_auroc' in col]
    norm_aupr = df[aupr_cols].sub(df['baseline_aupr'],axis=0)
    norm_auroc = df[auroc_cols].sub(df['baseline_auroc'],axis=0)
    norm_aupr.columns = ["norm_" + col for col in aupr_cols]
    norm_auroc.columns = ["norm_" + col for col in auroc_cols]
    final_df = df.join(norm_aupr).join(norm_auroc)
    final_df.sort('cluster_id', inplace=True)
    return(final_df)

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

df_list = []
prev_df = None
iseq = True
for filepath in os.listdir():
    if 'cluster_summary_within_c7_' in filepath and 'swp' not in filepath:
        final_df = parse_summary(filepath)
        # each summary is the result of running the tdrollers on a certain cluster
        # check if final df cluster ids is equal to prev df
        df_list.append(final_df)
## I want to get a positive correlation of the norm auroc columns


exp =df_list[0]
parsed_df_list = []
for df in df_list:
    if exp['cluster_id'].equals(df['cluster_id']):
        parsed_df_list.append(df)


big_df = pd.concat(parsed_df_list)
mean_df=big_df.groupby(level=0).mean()

t_list = [0,1,2,3,4,5,6,7,8,9,10]
m_list = [1,2]
param_list = []
for m in m_list:
    for t in t_list:
        info = pinf(lag_thresh=t, fill_na=True, median_thresh=m)
        inf_col = info.columns.tolist()
        temp_df = mean_df.merge(info, on='cluster_id')
        norm_cols = [col for col in mean_df.columns if 'norm' in col]
        norm_cor = temp_df.corr(method='spearman')['percent_lagged_y'].loc[norm_cols]
        if any(norm_cor>0):
            print('found one! t is ' + str(t))
            param_list.append((m,t))
        test_stat = 'lag_median'
        temp_df.sort(test_stat, inplace=True)
        thresh = temp_df[test_stat].describe()[5]
        not_lagged = temp_df[temp_df[test_stat] < thresh]
        lagged = temp_df[temp_df[test_stat] >= thresh]
        diff = lagged[norm_cols].mean() - not_lagged[norm_cols].mean()
        print(diff)
        if any(diff>0):
            print('yes! threshold',thresh)
            param_list.append((m,t))


print(param_list)
mean_df.to_csv('mean_cluster_summary_within_c7.csv', sep = '\t', index = False)
mean_df.corr()
mean_df.corr().to_csv('mean_corr_cluster_summary_within_c7.csv', sep = '\t', index = False)
