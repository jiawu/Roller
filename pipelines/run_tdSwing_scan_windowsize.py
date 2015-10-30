import Pipelines as pl

# td_window from min to max

# min_lag from min to max

# max_lag from min to max

# n_trees from 10-100-1000

# permutation_n from 10-100-1000

# all different lag methods


# optimization prioritization: let's do least sensitive to most sensitive.
# start with n_trees, then permutation_n, then td_window, then min/max lag

#we will identify the change in AUROC and AUPR in 100 iterations
#we will plot a bar graph showing mean, median, standard deviation/quartiles for each parameter

# saving the models for the iteration tests:
# to save the models for the iteration tests, we will save a dataframe (in the form of the final dataframe from Analyzer...) instead of a full model, because it is too computationally expensive, and as of this day, we are running out of room on QUEST.

data_folder = "/projects/p20519/roller_output/optimizing_window_size/RandomForest/ecoli"
output_path = "/home/jjw036/Swing/ecoli"
#target_dataset = "/projects/p20519/Swing/data/invitro/janes_timeseries.tsv"
target_dataset = "/projects/p20519/Swing/data/dream4/ecoli_timeseries.tsv"
my_statistic = 'aupr'
roc,pr = get_td_stats(target_dataset, 5)
pdb.set_trace()
td_score = pr

main(data_folder, output_path, target_dataset, 'auroc', roc)
main(data_folder, output_path, target_dataset, 'aupr', pr)



pl.get_td_stats


def get_td_stats(file_path, td_window = 6, min_lag=1, max_lag=4, n_trees=10, permutation_n=10, lag_method='median_median'): 

