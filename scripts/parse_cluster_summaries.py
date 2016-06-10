import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import pandas as pd
import pdb
from parse_cluster_info import main as pinf
import numpy as np
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


def plot_lag(df, idx):
    # separate into two groups, between and within groups
    within_clusters = df[df['parent_cluster'] == df['child_cluster']]
    between_clusters = df[df['parent_cluster'] != df['child_cluster']]
    
    plt.figure(1)
    #within_clusters_values = within_clusters['lag_median'].fillna(0).values
    #between_clusters_values = between_clusters['lag_median'].fillna(0).values
    within_clusters_values = within_clusters['Lag'].dropna().values
    between_clusters_values = between_clusters['Lag'].dropna().values

    within_weights= np.ones_like(within_clusters_values)/len(within_clusters_values)
    between_weights= np.ones_like(between_clusters_values)/len(between_clusters_values)
    ax_clust = fig.add_subplot(10,1,idx)
    bins = [x for x in np.arange(0,3,0.25)]
    ax_clust.hist(within_clusters_values,bins = bins, alpha = 0.5, label = 'Within', weights=within_weights, color = 'red')
    ax_clust.hist(between_clusters_values, bins=bins,alpha = 0.5, label = 'Between', weights=between_weights, color = 'blue')
    plt.legend(shadow=True, fancybox=True)


CLUSTER=26
df_list = []

prev_df = None
iseq = True
directory ='/projects/p20519/roller_output/cluster_summaries/' 
for filepath in os.listdir(directory):
    if 'cluster_summary_within_community_c'+str(CLUSTER)+'_' in filepath and 'swp' not in filepath:
        try:
            final_df = parse_summary(directory+filepath)
            # each summary is the result of running the tdrollers on a certain cluster
            # check if final df cluster ids is equal to prev df
            df_list.append(final_df)
        except ValueError:
            continue
## I want to get a positive correlation of the norm auroc columns


exp =df_list[0]
parsed_df_list = []
for df in df_list:
    if exp['cluster_id'].equals(df['cluster_id']):
        parsed_df_list.append(df)


big_df = pd.concat(parsed_df_list)
mean_df=big_df.groupby(level=0).mean()

t_list = [0]
#t_list = [5,7,8]
m_list = [10]
#t_list = [6]
#m_list = [2]
param_list = []

fig = plt.figure(1,figsize=(20,20))

counter = 1
for m in m_list:
    for t in t_list:
        info, lag_df = pinf(lag_thresh=t, fill_na=False, median_thresh=m, CLUSTER=CLUSTER, img_append=str(t))
        inf_col = info.columns.tolist()
        temp_df = mean_df.merge(info, on='cluster_id')
        plot_lag(lag_df, counter)
        
        norm_cols = [col for col in mean_df.columns if 'norm' in col]
        test_stat = 'percent_lagged_y'
        #test_stat = 'lag_mean'
        norm_cor = temp_df.corr(method='spearman')[test_stat].loc[norm_cols]
        
        print('norm_correlation :', norm_cor)
        print('t,m ',t,m)
        param_list.append((m,t))
        temp_df.sort(test_stat, inplace=True)
        #thresh = temp_df[test_stat].describe()[5]
        thresh = 0.1
        not_lagged = temp_df[temp_df[test_stat] < thresh]
        lagged = temp_df[temp_df[test_stat] >= thresh]
        diff = lagged[norm_cols].mean() - not_lagged[norm_cols].mean()
        print('diff improvement: ', diff)
        print('lagged not lagged',len(lagged), len(not_lagged))
        print('threshold',thresh)
        counter +=1


print(param_list)

info, lag_df = pinf(lag_thresh=t, fill_na=False, median_thresh=m, CLUSTER=CLUSTER)
inf_col = info.columns.tolist()
mean_df = mean_df.merge(info, on='cluster_id')

mean_df.to_csv('mean_cluster_summary_within_c'+str(CLUSTER)+'.csv', sep = '\t', index = False)
mean_df.corr()
mean_df.corr().to_csv('mean_corr_cluster_summary_within_c'+str(CLUSTER)+'.csv', sep = '\t', index = False)
plt.figure(1)
plt.savefig('CLUSTER_'+str(CLUSTER)+'_hist_between_within.png')
pdb.set_trace()
