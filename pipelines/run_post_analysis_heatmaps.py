import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#from mpl_toolkits.mplot3d import Axes3D

import os
import pandas as pd
from Roller.util.Evaluator import Evaluator
import pdb
import numpy as np
import kdpee

#get all pickle files
path="/projects/p20519/roller_output/optimizing_window_size/RandomForest/janes/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Roller/janes"

target_dataset =  '/projects/p20519/Roller/data/invitro/janes_timeseries.tsv'
img_suffix = "1"
gp_left = 0.2
gp_bottom = 0.1
gp_width = 0.7
gp_height = 0.2

padding = 0.01
numTFs = 200
dm_left = gp_left
dm_bottom = gp_bottom+gp_height+2*padding
dm_width = gp_width
box_height = 0.03
dm_height = box_height*numTFs

tableau20 = [ (152,223,138),(31, 119, 180), (174, 199, 232), (255, 127, 14), 
              (255, 187, 120),  (44, 160, 44), (255, 152, 150),  
              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
              (214,39,40)]
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.) 

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
entropies = []
window_start_list = []
window_end_list = []
window_auroc = []

window_size_list = []
fig = plt.figure()
colors = np.linspace(0,1,22)
raw_data_list = []

for file in filenames:
  full_path = path + file
  #print(full_path)
  try:
    roller_obj = pd.read_pickle(full_path)
  except EOFError:
    continue
  attributes = dir(roller_obj)
  if any("file_path" in attribute for attribute in attributes):
    counter += 1
    #print(str(counter) + " out of " + str(nfiles))
    key = roller_obj.file_path
    if key == target_dataset:
      window_size = roller_obj.window_width
      raw_data_list.append(roller_obj.raw_data)
      raw_data = roller_obj.raw_data
      #roller_obj.average_rank(rank_by='stability', ascending = False)
#      roller_obj.average_rank(rank_by='importance', ascending = True)
      
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
          gold_standard = roller_obj.file_path.replace("timeseries.tsv","goldstandard.tsv")
          evaluator = Evaluator(gold_standard, sep="\t")
          edge_cutoff = len(evaluator.gs_flat)
          sorted_edge_list = window.results_table
          #sorted_edge_list.sort(['stability'], ascending=[True], inplace=True)
          sorted_edge_list.sort(['p_value'], ascending=[True], inplace=True)
          #sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p-means'])]
          sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p_value'])]
          precision, recall, aupr = evaluator.calc_pr(sorted_edge_list)
          precision = precision.tolist()
          recall = recall.tolist()
          aupr = aupr.tolist()
          tpr, fpr, auroc = evaluator.calc_roc(sorted_edge_list)
          auroc = auroc.tolist()
          tpr = tpr.tolist()
          fpr = fpr.tolist()

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
        #print(max_time-min_time)
        #print(window_size)
        #print(window_auroc)
        #print(window_start_list)
        #ax.bar(window_start_list, window_auroc, zs=window_size,zorder=window_size, zdir='y',alpha=0.8,width=(max_time-min_time), color = plt.cm.RdYlBu(colors[window_size]))


#print(window_size_list)
size_to_auroc = zip(window_size_list,window_start_list, window_auroc)
size_to_aupr = zip(window_size_list,window_start_list, aupr_list2)
pdb.set_trace()
f = plt.figure(figsize=(10,10))
axarr2 = f.add_axes([gp_left, gp_bottom, gp_width, gp_height])
axarr1 = f.add_axes([dm_left, dm_bottom, dm_width, dm_height])

#figure, heatmap = Grapher.generate_heatmap_from_df(self.roller_results.raw_data)
#axarr1 = figure
time_vector = raw_data['Time'].unique()
nrow,ncol=raw_data.shape
nrepeats = nrow/len(time_vector)
#create group assignment
groups = [x for x in range(0,nrepeats) for i in range(0,len(time_vector))]
raw_data['Group'] = groups
sorted = raw_data.sort(['Time','Group'])
heatmap_values = pd.DataFrame()
col_labels = map(int,time_vector)
partial_row_labels = raw_data.drop(['Time','Group'],1).columns.values.tolist()
row_labels = [i for x in range(0,nrepeats) for i in partial_row_labels]
for time in time_vector: 
    target_set = sorted[sorted['Time']==time].drop(['Time','Group'],1)
    vector = target_set.values.flatten()
    heatmap_values[time] = vector
my_axis =axarr1.matshow(heatmap_values,cmap=matplotlib.cm.RdBu,aspect='auto')
my_axi = my_axis.get_axes()
axarr1.set_xlabel('Time')
axarr1.xaxis.set_label_coords(0.5, 1.08)

axarr1.set_yticks(np.arange(heatmap_values.shape[0]))
axarr1.yaxis.set_ticks_position('left')
axarr1.set_yticklabels(row_labels)
axarr1.set_xticks(np.arange(heatmap_values.shape[1]))
axarr1.xaxis.set_ticks_position('top')
xlabelsL = axarr1.set_xticklabels(col_labels)
for label in xlabelsL:
    label.set_rotation(90)
    label.set_fontsize(8)
for label in (axarr1.get_yticklabels()):
    label.set_fontsize(8)
for l in axarr1.get_xticklines() + axarr1.get_yticklines():
    l.set_markersize(0)
unique_window_sizes = list(set(window_size_list))
for color_index,window_size in enumerate(unique_window_sizes):
    if window_size == roller_obj.overall_width:  
      if window_size == max(unique_window_sizes):
          #add horizontal line instead of adding dot plot
          line_start_list = [min(time_vector), max(time_vector)]
          line_auroc = [item[2] for item in size_to_auroc if item[0] == window_size]
          line_auroc = [line_auroc[0],line_auroc[0]]
          axarr2.plot(map(int,line_start_list), line_auroc, linestyle='--',color =
              tableau20[color_index], label = "WS "+str(window_size))
    else:
        line_start_list = [item[1] for item in size_to_auroc if item[0] == window_size]
        line_auroc = [item[2] for item in size_to_auroc if item[0] == window_size]
        axarr2.plot(map(int,line_start_list), line_auroc, 'o', linestyle='-',color =
            tableau20[color_index], label = "WS "+str(window_size))
axarr2.xaxis.set_ticks_position('bottom')

pdb.set_trace()
axarr2.get_xaxis().get_major_formatter().labelOnlyBase = False
line_ticks = np.arange(roller_obj.time_vec[0],roller_obj.time_vec[-1],200)
axarr2.xaxis.set_ticks(line_ticks)

#axarr2.xaxis.set_ticks(np.arange(0,1050,50))
xlabels = axarr2.get_xticklabels()
ylabels = axarr2.get_yticklabels()
for label in xlabels:
    label.set_rotation(90)
    label.set_fontsize(8)
for label in (axarr2.get_yticklabels()):
    label.set_fontsize(8)
for l in axarr2.get_xticklines() + axarr2.get_yticklines():
    l.set_markersize(0)
box = axarr2.get_position()
axarr2.set_position([box.x0, box.y0 + box.height * 0.2,
                   box.width, box.height * 0.8])

axarr2.legend(fontsize=8,bbox_to_anchor=(0.5, -0.2), loc='upper center',ncol=7,fancybox=True, shadow=True)
axarr1.set_ylabel('Genes (multiple perturbations)')
axarr2.set_ylabel('AUROC')
image_save = image_file_path + "_windowsize_RF_heatmap_" + img_suffix + ".png"
f.savefig(image_save,format="png")
print("AUROC:")
size_to_auroc.sort(key=lambda tup:tup[0],reverse=True)
print(size_to_auroc[0:5])
size_to_auroc.sort(key=lambda tup:tup[2],reverse=True)
print(size_to_auroc[0:5])
print("AUPR:")

size_to_aupr.sort(key=lambda tup:tup[0],reverse=True)
print(size_to_aupr[0:5])
size_to_aupr.sort(key=lambda tup:tup[2],reverse=True)
print(size_to_aupr[0:5])

