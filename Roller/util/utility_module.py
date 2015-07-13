import pandas as pd
import scipy.stats as ss
import sklearn.metrics as skmet
import numpy as np

def get_test_set(window_raw_data, roller_raw_data):
    roller_vec = roller_raw_data['Time'].unique()
    window_vec = window_raw_data['Time'].unique()
    test_set_vec = np.setdiff1d(roller_vec,window_vec)
    test_data = roller_raw_data.loc[roller_raw_data['Time'].isin(test_set_vec)].drop('Time',1)
    return(test_data)


def get_cragging_scores(model, predictor, response_true):
    response_pred = model.predict(predictor)
    scores = {}
    scores['ev'] = skmet.explained_variance_score(response_true, response_pred)
    scores['mae'] = skmet.mean_absolute_error(response_true, response_pred)
    scores['mse'] = skmet.mean_squared_error(response_true, response_pred)
    #scores['medae'] = skmet.median_absolute_error(response_true, response_pred)
    scores['r2'] = skmet.r2_score(response_true, response_pred)
    return(scores)

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

def point_slope(x1,y1, x2,y2):
    slope = (y2-y1)/float(x2-x1)
    return slope

def elbow_criteria(x,y):
    x = np.array(x)
    y = np.array(y)
    # Slope between elbow endpoints
    m1 = point_slope(x[0], y[0], x[-1], y[-1])
    # Intercept
    b1 = y[0] - m1*x[0]

    # Slope for perpendicular lines
    m2 = -1/m1

    # Calculate intercepts for perpendicular lines that go through data point
    b_array = y-m2*x
    x_perp = (b_array-b1)/(m1-m2)
    y_perp = m1*x_perp+b1

    # Calculate where the maximum distance to a line connecting endpoints is
    distances = np.sqrt((x_perp-x)**2+(y_perp-y)**2)
    index_max = np.where(distances==np.max(distances))[0][0]
    elbow_x = x[index_max]
    elbow_y = y[index_max]
    return elbow_x, elbow_y