import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import os
import pandas as pd
from Roller.util.Evaluator import Evaluator
import pdb
import numpy as np
import kdpee

#get all pickle files
path="/projects/p20519/Roller_outputs_RF/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Roller/aggregated"

target_dataset =  "data/dream4/insilico_size10_1_timeseries.tsv"

#todo: calculate distance between two edge lists
#average models: calculate AUROC between averaged models
#single line graph: window size, auroc


dataset_dict = {}
best_alpha_list = []
aupr_list2 = []
auroc_list2 = []
tpr_list2 = []
fpr_list2 = []
precision_list2 = []
recall_list2 = []
entropies = []
window_start_list = []
window_end_list = []
window_auroc = []

a_aupr_list = []
a_auroc_list = []
a_tpr_list = []
a_fpr_list = []
a_precision_list = []
a_recall_list = []
a_window_size_list =[]

window_size_list = []
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
colors = np.linspace(0,1,22)

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
      averaged_edge_list = roller_obj.average_rank(rank_by='p_value', ascending = False)
      averaged_edge_list.sort(['mean-rank'], ascending=[False], inplace=True)
      averaged_edge_list = averaged_edge_list[np.isfinite(averaged_edge_list['mean-rank'])]
      gold_standard = roller_obj.file_path.replace("timeseries.tsv","goldstandard.tsv")
      evaluator = Evaluator(gold_standard, sep="\t")
      pdb.set_trace()
      a_precision, a_recall, a_aupr = evaluator.calc_pr(averaged_edge_list)
      a_tpr, a_fpr, a_auroc = evaluator.calc_roc(averaged_edge_list)
      
      a_window_size_list.append(window_size)
      a_precision_list.append(a_precision[14])
      a_recall_list.append(a_recall[14])
      a_aupr_list.append(a_aupr[-1])
      a_tpr_list.append(a_tpr[14])
      a_fpr_list.append(a_fpr[14])
      a_auroc_list.append(a_auroc[-1])

      # create barplots for each Roller of a certain window size
      # window size = bar plot
      # then group all barplots into a 3D bar plot
      if window_size not in window_size_list:
        for nwindow,window in enumerate(roller_obj.window_list):
          min_time = min(window.raw_data['Time'].unique())
          window_start_list.append(min_time)
          max_time = max(window.raw_data['Time'].unique())
      #other graphs correlating window size to statistics
            #have to redefine gold standard here because it wasn't saved in the object
          edge_cutoff = len(evaluator.gs_flat)
          sorted_edge_list = window.results_table
          sorted_edge_list.sort(['p_value'], ascending=[True], inplace=True)
          sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p_value'])]
          precision, recall, aupr = evaluator.calc_pr(sorted_edge_list)
          tpr, fpr, auroc = evaluator.calc_roc(sorted_edge_list)
          
          window_auroc.append(auroc[-1])

          window_size_list.append(window_size)
          precision_list2.append(precision[14])
          recall_list2.append(recall[14])
          aupr_list2.append(aupr[-1])
          tpr_list2.append(tpr[14])
          fpr_list2.append(fpr[14])
          auroc_list2.append(auroc[-1])

          window_transposed = window.df.values.transpose()
          window_entropy = kdpee.kdpee(window_transposed)
          entropies.append(window_entropy)
        print(max_time-min_time)
        print(window_size)
        print(window_auroc)
        print(window_start_list)
        #ax.bar(window_start_list, window_auroc, zs=window_size,zorder=window_size, zdir='y',alpha=0.8,width=(max_time-min_time), color = plt.cm.RdYlBu(colors[window_size]))


print(window_size_list)

ax.set_xlabel('Time')
ax.set_ylabel('Window Size')
ax.set_zlabel('AUROC')
image_save = image_file_path + "_3D_RF.png"
fig.savefig(image_save)

beforenorm = pd.read_csv(roller_obj.file_path, sep='\t')
beforenorm = beforenorm.dropna(axis=0,how='all')
ax3 = beforenorm.plot(x='Time',y="G1",kind="scatter")
fig3 = ax3.get_figure()
image_save = image_file_path + "_time_series.png"
fig3.savefig(image_save)

ax2 = roller_obj.raw_data.plot(x='Time', y="G1", kind="scatter")
fig2 = ax2.get_figure()
image_save = image_file_path + "_time_series_zscore.png"
fig2.savefig(image_save)

result_titles = ["precision","recall","aupr","tpr","fpr","auroc", "entropy"]

result_list2 = [precision_list2, recall_list2, aupr_list2, tpr_list2, fpr_list2,
    auroc_list2, entropies]

for count,result in enumerate(result_list2):
  print(window_size_list)
  print(result)
  fig = plt.figure(6+count)
  plt.scatter(window_size_list, result)
  title_string = "Window Size, RF, One Dataset"
  plt.title(title_string)
  plt.xlabel('Window Size')
  plt.ylabel(result_titles[count])
  image_save = image_file_path + "_windowsize_RF_" + str(result_titles[count]) + ".png"
  fig.savefig(image_save)

