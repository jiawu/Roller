import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import pandas as pd
from Swing_old.util.Evaluator import Evaluator
import pdb
import numpy as np
import kdpee
from scipy.stats.stats import pearsonr

def aggregate_scores(dict_list):
  aggr_dict = {}
  for key in dict_list[0].keys():    
    new_list = [x[key] for x in dict_list]
    average = sum(new_list)/len(new_list)
    aggr_dict[key] = average

  return(aggr_dict)

#get all pickle files
path="/projects/p20519/roller_output/optimizing_window_size/RandomForest/ecoli/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Swing/aggregated"


#todo: calculate distance between two edge lists
#average models: calculate AUROC between averaged models
#single line graph: window size, auroc
gp_left = 0.2
gp_bottom = 0.1
gp_width = 0.7
gp_height = 0.2

padding = 0.01
numTFs = 10
dm_left = 0.2
dm_bottom = 0.1+0.2+2*0.01
dm_width = 0.7
dm_height = 0.03*10

target_dataset =  "data/dream4/ecoli_timeseries.tsv"
img_suffix = "2"
#dataset_list = ["5"]
dataset_list = ["2"]


for dataset_counter in dataset_list:
  window_start_list = []
  auroc_list = []
  rate_list = []
  average_list = []
  aggr_training_list = []
  aggr_test_list = []
  checked_windows = []
  img_suffix = str(dataset_counter)
  target_dataset =  "data/dream4/ecoli_timeseries.tsv"
  for file in filenames:
      full_path = path + file
      print(full_path)
      try:
          roller_obj = pd.read_pickle(full_path)
      except EOFError:
          continue
      attributes = dir(roller_obj)
      if any("file_path" in attribute for attribute in attributes):
        counter += 1
        print(str(counter) + " out of " + str(nfiles))
        key = roller_obj.file_path
        if key == target_dataset:
          window_size = roller_obj.window_width
          target_roller = roller_obj
          if window_size not in checked_windows:
            if window_size !=0:
              checked_windows.append(window_size)
              for window in roller_obj.window_list:
                min_time = min(window.raw_data['Time'].unique())
                window_start_list.append(min_time)
                gold_standard = target_roller.file_path.replace("timeseries.tsv","goldstandard.tsv")
                evaluator = Evaluator(gold_standard,sep="\t")
                sorted_edge_list = window.results_table
                pdb.set_trace()
                sorted_edge_list.sort(['importance'], ascending=[False], inplace=True)
                sorted_edge_list =sorted_edge_list[np.isfinite(sorted_edge_list['p_value'])]
                #sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p-means'])]
                #pdb.set_trace()
                #sorted_edge_list.sort(['stability'], ascending=[False], inplace=True)
                tpr,fpr, auroc = evaluator.calc_roc(sorted_edge_list)
                auroc_list.append(auroc[-1])

                #posthoc analysis of rates and averages
                training_scores = window.training_scores[0:10]
                test_scores = window.test_scores[0:10]
                aggr_training= aggregate_scores(training_scores)
                aggr_test= aggregate_scores(test_scores)
                aggr_training_list.append(aggr_training)
                aggr_test_list.append(aggr_test)
                rate_dict = window.get_rate_analysis(1)
                rate_list.append(rate_dict)

                average_list.append(window.get_average())

  #scatterplot
  aggr_mae = [rate['mae'].mean() for rate in aggr_test_list] 
  aggr_mse = [rate['mse'].mean() for rate in aggr_test_list] 
  aggr_ev = [rate['ev'].mean() for rate in aggr_test_list] 
  #aggr_r2 = [rate['r2'].mean() for rate in aggr_test_list] 
  aggrt_mae = [rate['mae'].mean() for rate in aggr_training_list] 
  aggrt_mse = [rate['mse'].mean() for rate in aggr_training_list] 
  aggrt_ev = [rate['ev'].mean() for rate in aggr_training_list] 
  #aggrt_r2 = [rate['r2'].mean() for rate in aggr_training_list] 
  axis_font = {'size':'28'}
  
  pdb.set_trace()
  f = plt.figure(figsize=(10,10))
  ax = f.add_subplot(1,1,1)
  from pylab import rc, gca
  ax = gca()
  ax.scatter(auroc_list, aggr_mae, s = 50, c="blue")
  #ax.axhline(linewidth=4)
  #ax.axvline(linewidth=4)
  plt.xlabel("AUROC", **axis_font)
  plt.ylabel("MAE", **axis_font)
  plt.xticks(fontsize = 16, fontweight='bold')
  plt.yticks(fontsize = 16, fontweight='bold')
  rc('axes', linewidth=2)
  fontsize=16
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
  plt.tight_layout()
  #ax.set_xlim([0.5,0.9])
  #ax.set_ylim([-1.0,1.0])
  image_save = "scatter_mae.png"
  f.savefig(image_save,format="png")
  print(pearsonr(auroc_list, aggr_mae))
  
  ## get the statistics for each gene
  row_labels = window.genes
  row_labels = np.append(row_labels, 'Aggr')
  heatmap_values = pd.DataFrame()

  col_labels = ['Rate (mean)', 'Rate (median)', 'Rate (min)', 'Rate (max)', 'Rate (std)']

  ## dataframe labels: mean values, 
  for index, gene in enumerate(window.genes):
      #correlate auroc with gene mean
      #auroc of a window is in auroc_list
      #auroc_list should have the same length as the rate_dict_list
  
      #get the mean for gene 1
      item_mean = [rate_dict['mean'][index] for rate_dict in rate_list]
      item_min = [rate_dict['min'][index] for rate_dict in rate_list]
      item_max = [rate_dict['max'][index] for rate_dict in rate_list]
      item_median = [rate_dict['median'][index] for rate_dict in rate_list]
      item_std = [rate_dict['std'][index] for rate_dict in rate_list]

      
      mean_c = pearsonr(item_mean, auroc_list)[0]
      median_c = pearsonr(item_median, auroc_list)[0]
      min_c = pearsonr(item_min, auroc_list)[0]
      max_c = pearsonr(item_max, auroc_list)[0]
      std_c = pearsonr(item_std, auroc_list)[0]
      
      vector = [mean_c, median_c, min_c, max_c, std_c]
      heatmap_values[gene] = vector
  """all_means = [rate['all_rates'].mean() for rate in rate_list] 
  all_medians = [np.median(rate['all_rates']) for rate in rate_list] 
  all_mins = [rate['all_rates'].min() for rate in rate_list] 
  all_maxs = [rate['all_rates'].max() for rate in rate_list] 
  all_std = [rate['all_rates'].std(ddof=1) for rate in rate_list]
  all_means_c = pearsonr(all_means, auroc_list)[0]
  all_medians_c = pearsonr(all_medians, auroc_list)[0]
  all_mins_c = pearsonr(all_mins, auroc_list)[0]
  all_maxs_c = pearsonr(all_maxs, auroc_list)[0]
  all_std_c = pearsonr(all_std, auroc_list)[0]
  """
  all_means = [rate['mae'].mean() for rate in aggr_test_list] 
  all_medians = [rate['mse'].mean() for rate in aggr_test_list] 
  all_mins = [rate['all_rates'].min() for rate in rate_list] 
  all_maxs = [rate['all_rates'].max() for rate in rate_list] 
  all_std = [rate['all_rates'].std(ddof=1) for rate in rate_list]
  all_means_c = pearsonr(all_means, auroc_list)[0]
  all_medians_c = pearsonr(all_medians, auroc_list)[0]
  all_mins_c = pearsonr(all_mins, auroc_list)[0]
  all_maxs_c = pearsonr(all_maxs, auroc_list)[0]
  all_std_c = pearsonr(all_std, auroc_list)[0]

  all_vector = [all_means_c, all_medians_c, all_mins_c, all_maxs_c,all_std_c]
  #rates_list = window.get_rates(1)
  
  heatmap_values['all'] = all_vector
  f = plt.figure(figsize=(10,10))
  axarr1 = f.add_axes([dm_left, dm_bottom, dm_width, dm_height])
  my_axis = axarr1.matshow(heatmap_values,cmap=matplotlib.cm.RdBu,aspect='auto', vmin=-1, vmax=1)
  my_axi = my_axis.get_axes()
  axarr1.set_xlabel('Genes')
  axarr1.xaxis.set_label_coords(0.5,1.15)
  axarr1.yaxis.set_ticks_position('left')
  axarr1.set_yticklabels(col_labels)
  
  axarr1.xaxis.set_ticks_position('top')
  xlabelsL = axarr1.set_xticklabels(row_labels)
  for label in xlabelsL:
      label.set_rotation(90)
      label.set_fontsize(8)
  for label in (axarr1.get_yticklabels()):
      label.set_fontsize(8)
  axarr1.set_xticks(np.arange(heatmap_values.shape[1]))
  axarr1.set_yticks(np.arange(heatmap_values.shape[0]))
  axarr1.grid("off")
  for l in axarr1.get_xticklines() + axarr1.get_yticklines():
      l.set_markersize(0)


  image_save = "test_heat_map"+dataset_counter+".png"
  f.savefig(image_save,format="png")

