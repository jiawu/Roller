import pandas as pd
import scipy.stats as ss

def create_3D_linked_list(labels, numpy_array_3D, value_label):
    """returns a panel with interaction (x-axis) - value (y axis) - time (Z axis)"""
    windows_n = numpy_array_3D.shape[2]
    linked_list_3D ={}

    for i in xrange(windows_n):
        target_2D_array = numpy_array_3D[:,:,i]
        linked_list = create_linked_list(labels, target_2D_array, value_label)
        linked_list_3D[i] = linked_list
    return pd.Panel(linked_list_3D)

def create_linked_list(labels, numpy_array_2D, value_label):
    """labels and array should be in row-major order"""
    linked_list = pd.DataFrame({'regulator-target':labels, value_label:numpy_array_2D.flatten()})
    return linked_list

def average_rank(ranked_result_list, col_string):
    """finds the average rank and standard deviation throughout time"""
    aggregate_ranks = []
    for nth_window in ranked_result_list:
        aggregate_ranks.append(nth_window[[col_string, 'regulator-target']])
    #now merge the panels in an interesting loop. The merge function insures the keys are always matched up correctly. Zip would've worked too.
    left_df = aggregate_ranks[0] #initialize the left_df.
    left_df.columns= [col_string+"_0", 'regulator-target']
    for window_index in range(1,len(aggregate_ranks)):
        right_df = aggregate_ranks[window_index]
        right_df.columns= [col_string+"_"+str(window_index), 'regulator-target']
        left_df = left_df.merge(right_df,on = 'regulator-target')

    aggr_ranks = left_df.drop(['regulator-target'], axis = 1)
    #assign to temporary variables to prevent the calc columns to be involved in other calculations
    range_col = zip(aggr_ranks.min(axis = 1), aggr_ranks.max(axis = 1))
    mean_col = aggr_ranks.mean(axis = 1)
    median_col = aggr_ranks.median(axis = 1)
    sd_col = aggr_ranks.std(axis = 1, ddof=1)

    aggr_ranks['range'] = range_col
    aggr_ranks['mean-rank'] = mean_col
    aggr_ranks['median-rank'] = median_col
    aggr_ranks['sd-rank'] = sd_col
    aggr_ranks['regulator-target'] = left_df['regulator-target']
    return(aggr_ranks)

def rank_results_3D(result_list, col_string, ascending = True):
    """input: list of result pandas dfs, column name. output: each time window is sorted by column name, most significant to least"""
    rank_column_name = col_string + "-rank"
    for nth_window in result_list:
        nth_window[rank_column_name] = nth_window[col_string].rank(method="dense", ascending = ascending)
    return result_list

def rank_index(vector):
        return [vector.index(x) for x in sorted(range(vector), key=vector.__getitem__)]
