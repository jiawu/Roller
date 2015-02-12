import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import pandas as pd
from Roller.util.Evaluator import Evaluator
#get all pickle files
path="/projects/p20519/sorted_Rollers/data/dream4/insilico_size10_1_timeseries/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Roller/aggregated"

target_dataset =  "data/dream4/insilico_size10_1_timeseries.tsv"

dataset_dict = {}
best_alpha_list = []
aupr_list = []
auroc_list = []
tpr_list = []
fpr_list = []
precision_list = []
recall_list = []
aupr_list2 = []
auroc_list2 = []
tpr_list2 = []
fpr_list2 = []
precision_list2 = []
recall_list2 = []

window_size_list = []

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
      for nwindow,window in enumerate(roller_obj.window_list):
        print(nwindow)
    #other graphs correlating window size to statistics
        if window_size == 21:    
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
          
        if alpha_dist < 0.05:
          window_size_list.append(window_size)
          precision_list2.append(precision[-1])
          recall_list2.append(recall[-1])
          aupr_list2.append(aupr[-1])
          tpr_list2.append(tpr[-1])
          fpr_list2.append(fpr[-1])
          auroc_list2.append(auroc)

print(window_size_list)

result_list = [precision_list, recall_list, aupr_list, tpr_list, fpr_list, auroc_list]
result_titles = ["precision","recall","aupr","tpr","fpr","auroc"]

result_list2 = [precision_list2, recall_list2, aupr_list2, tpr_list2, fpr_list2, auroc_list2]

for count,result in enumerate(result_list):
  fig = plt.figure(count)
  plt.scatter(best_alpha_list, result)
  title_string = "Window Size 21, One Dataset"
  plt.title(title_string)
  plt.xlabel('Approx Distance from Best Alpha')
  plt.ylabel(result_titles[count])
  image_save = image_file_path + "_best_alpha_" + str(result_titles[count]) + ".png"
  fig.savefig(image_save)

for count,result in enumerate(result_list2):
  fig = plt.figure(6+count)
  plt.scatter(window_size_list, result)
  title_string = "Window Size, Optimal Alpha, One Dataset"
  plt.title(title_string)
  plt.xlabel('Window Size')
  plt.ylabel(result_titles[count])
  image_save = image_file_path + "_windowsize_" + str(result_titles[count]) + ".png"
  fig.savefig(image_save)

