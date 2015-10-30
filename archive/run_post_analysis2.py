import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import pandas as pd
from Swing_old.util.Evaluator import Evaluator
#get all pickle files
path="/projects/p20519/Roller_outputs/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Swing/aggregated"

dataset_dict = {}

for file in filenames:
  full_path = path + file
  print(full_path)
  roller_obj = pd.read_pickle(full_path)
  attributes = dir(roller_obj)
  if any("file_path" in attribute for attribute in attributes):
    counter += 1
    print(str(counter) + " out of " + str(nfiles))
    obj_list.append(roller_obj)
    
    key = roller_obj.file_path
    if key not in dataset_dict:
      dataset_dict[key] = []

    dataset_dict[key].append(roller_obj)
    
print(len(obj_list))
print(dataset_dict.keys())
#optimum alpha. hypothesis: best aupr and best auroc near optimum alpha.
#distance from alpha, aupr auroc

# x axis is distance from best alpha
# y axis is aupr, auroc, tpr, fpr

best_alpha_list = []
aupr_list = []
auroc_list = []
tpr_list = []
fpr_list = []
precision_list = []
recall_list = []

window_size_list = []


for dataset in dataset_dict.keys():
  #print aupr and auroc
  #same window size? 
  if dataset == dataset_dict.keys()[1]:
    for roller_obj in dataset_dict[dataset]:    
      window_size = roller_obj.window_width
      for window in roller_obj.window_list:
    #other graphs correlating window size to statistics
        if window_size == 20:    
          window_size_list.append(window_size)
          current_alpha = window.alpha
          alpha_table = window.cv_table["alpha"]
          best_a_index = window.cv_table["positive_q2"].idxmax(1)
          best_alpha = alpha_table[best_a_index]
          
          alpha_dist = best_alpha - current_alpha
          
          best_alpha_list.append(alpha_dist)

          #have to redefine gold standard here because it wasn't saved in the object
          gold_standard = roller_obj.file_path.replace("timeseries.tsv","goldstandard.tsv")
          evaluator = Evaluator(gold_standard, sep="\t")
          edge_cutoff = len(evaluator.gs_flat)
          sorted_edge_list = window.results_table
          precision, recall, aupr = evaluator.calc_pr(sorted_edge_list[0:edge_cutoff+1])
          tpr, fpr, auroc = evaluator.calc_roc(sorted_edge_list[0:edge_cutoff+1])
          precision_list.append(precision[-1])
          recall_list.append(recall[-1])
          aupr_list.append(aupr[-1])
          tpr_list.append(tpr[-1])
          fpr_list.append(fpr[-1])
          auroc_list.append(auroc)

print(window_size_list)

result_list = [precision_list, recall_list, aupr_list, tpr_list, fpr_list, auroc_list]
result_titles = ["precision","recall","aupr","tpr","fpr","auroc"]
fig = plt.figure(100)
plt.scatter(window_size_list, precision_list)
title_string = "Window Size 21, One Dataset"
plt.title(title_string)
plt.xlabel('Window Size')
plt.ylabel(result_titles[0])
image_save = image_file_path + "_windowsize_" + str(result_titles[0]) + ".png"
fig.savefig(image_save)


for count,result in enumerate(result_list):
  fig = plt.figure(count)
  plt.scatter(best_alpha_list, result)
  title_string = "Window Size 21, One Dataset"
  plt.title(title_string)
  plt.xlabel('Approx Distance from Best Alpha')
  plt.ylabel(result_titles[count])
  image_save = image_file_path + "_" + str(result_titles[count]) + ".png"
  fig.savefig(image_save)

for count,result in enumerate(result_list):
  fig = plt.figure(6+count)
  plt.scatter(window_size_list, result)
  title_string = "Window Size 21, One Dataset"
  plt.title(title_string)
  plt.xlabel('Window Size')
  plt.ylabel(result_titles[count])
  image_save = image_file_path + "_windowsize_" + str(result_titles[count]) + ".png"
  fig.savefig(image_save)

