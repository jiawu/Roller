import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#from mpl_toolkits.mplot3d import Axes3D

import os
import pandas as pd
from Swing.util.Evaluator import Evaluator
import pdb
import numpy as np
import kdpee
#import seaborn
#get all pickle files
path="/projects/p20519/roller_output/optimizing_window_size/RandomForest/janes/"
filenames = next(os.walk(path))[2]
nfiles = len(filenames)
#organize pickled objects by dataset analyzed
obj_list = []
counter = 0
image_file_path = "/home/jjw036/Swing/aggregated_janes_RF"


#todo: calculate distance between two edge lists
#average models: calculate AUROC between averaged models
#single line graph: window size, auroc
gp_left = 0.2
gp_bottom = 0.1
gp_width = 0.7
gp_height = 0.2

padding = 0.01
numTFs = 20
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



fig = plt.figure()

colors = np.linspace(0,1,22)

overall_df = pd.DataFrame()

target_dataset =  "data/invitro/janes_timeseries.tsv"
img_suffix = "1"
#dataset_list = ["1"]

dataset_list = ["1"]
for dataset_counter in dataset_list:
  edge_list = []
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
  window_percentage_optimal = []
  #percentage of timepoints in the optimal region

  a_aupr_list = []
  a_auroc_list = []
  a_tpr_list = []
  a_fpr_list = []
  a_precision_list = []
  a_recall_list = []
  a_window_size_list =[]

  window_size_list = []
  img_suffix = str(dataset_counter)
  target_dataset =  "data/invitro/janes_timeseries.tsv"
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
          pdb.set_trace()
          target_roller = roller_obj

          window_size = roller_obj.window_width
          averaged_edge_list = roller_obj.average_rank(rank_by='importance',
              ascending = False)
          averaged_edge_list.sort(['mean-rank'], ascending=[False], inplace=True)
          averaged_edge_list = averaged_edge_list[np.isfinite(averaged_edge_list['mean-rank'])]
          gold_standard = roller_obj.file_path.replace("timeseries.tsv","goldstandard.tsv")
          evaluator = Evaluator(gold_standard, sep="\t")
          a_precision, a_recall, a_aupr = evaluator.calc_pr(averaged_edge_list)
          a_tpr, a_fpr, a_auroc = evaluator.calc_roc(averaged_edge_list)
          
          a_window_size_list.append(window_size)
          a_precision_list.append(a_precision[14])
          a_recall_list.append(a_recall[14])
          a_aupr_list.append(a_aupr[-1])
          a_tpr_list.append(a_tpr[14])
          a_fpr_list.append(a_fpr[14])
          a_auroc_list.append(a_auroc[-1])

          # create barplots for each Swing of a certain window size
          # window size = bar plot
          # then group all barplots into a 3D bar plot
          optimal_time = target_roller.raw_data['Time'].unique()[11:]
          if window_size not in window_size_list:
            for nwindow,window in enumerate(roller_obj.window_list):
              time = window.raw_data['Time'].unique()
              #find points in optimal region
              points_optimal_region = np.intersect1d(optimal_time, time)
              #percentage in optimal region
              percent_optimal_region = float(len(points_optimal_region))/len(time)
              window_percentage_optimal.append(percent_optimal_region)
              min_time = min(window.raw_data['Time'].unique())
              start_time = roller_obj.raw_data['Time'].unique()[window.nth_window]
              window_start_list.append(start_time)
              max_time = max(window.raw_data['Time'].unique())
          #other graphs correlating window size to statistics
                #have to redefine gold standard here because it wasn't saved in the object
              edge_cutoff = len(evaluator.gs_flat)
              sorted_edge_list = window.results_table
              #sorted_edge_list.sort(['stability'], ascending=[], inplace=True)
              sorted_edge_list = sorted_edge_list[np.isfinite(sorted_edge_list['p_value'])]
              edge_list.append(sorted_edge_list)
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

  network_df = pd.DataFrame({"average_auroc":a_auroc_list,"window_size": a_window_size_list})
  network_df['Network'] = "Network"+dataset_counter
  overall_df = overall_df.append(network_df)
  a_auroc_list = []
  a_window_size_list = []
  print(window_size_list)
  """ax.set_xlabel('Time')
  ax.set_ylabel('Window Size')
  ax.set_zlabel('AUROC')
  image_save = image_file_path + "_3D_lasso.png"
  fig.savefig(image_save)
  """
  figure_networks = plt.figure()
  #figure_networks=seaborn.lmplot("window_size","average_auroc",overall_df,col="Network",
  #    fit_reg=False)
  image_save = image_file_path + "_lasso_averages.png"
  figure_networks.savefig(image_save)

  optimal_correlation,ax = plt.subplots()
  ax.scatter(window_percentage_optimal, auroc_list2,marker='o', color='blue',edgecolor='none')
  #fig_corr=optimal_correlation.get_figure()
  optimal_correlation.tight_layout()
  image_save = image_file_path + "_correlation.png"
  optimal_correlation.savefig(image_save)
  
  pdb.set_trace()
  optimal_time = target_roller.raw_data['Time'].unique()[11:]
  
  #plot GvsG
  
  #first get all the best edges from the best windows
  auroc_edges = zip(auroc_list2, edge_list)
  auroc_edges.sort(key = lambda t: t[0], reverse=True)
  #take 5 edges from top 10 windows
  best_edges = [edge_list[1]['regulator-target'][0:6].tolist() for edge_list in auroc_edges[0:10]]
  best_edges = [item for sublist in best_edges for item in sublist]
  seen = set()
  unique_edges = [item for item in best_edges if item[1] not in seen and not seen.add(item[1])]
  opt_series = target_roller.raw_data[target_roller.raw_data['Time']>500]
  subopt_series = target_roller.raw_data[target_roller.raw_data['Time']<550]
  for tuple in unique_edges:
      regulator = tuple[0]
      target = tuple[1]
      #get values from regulator
      reg_o = opt_series[regulator].tolist()
      targ_o = opt_series[target].tolist()
      reg_s = subopt_series[regulator].tolist()
      targ_s = subopt_series[target].tolist()
      gvg,ax = plt.subplots()
      opt = ax.scatter(reg_o, targ_o,marker='o', color='blue',edgecolor='none',
          label="Optimal Time Window")
      sub = ax.scatter(reg_s, targ_s,marker='o', color='red',edgecolor='none',
          label="Suboptimal Time Window")
      ax.set_xlabel(regulator+" Expression (z-scored)")
      ax.set_ylabel(target+" Expression (z-scored)")
      ax.set_title(regulator+" vs " + target)
      ax.legend([opt, sub], ['Optimal Time Window', 'Suboptimal Time Window'])
      #fig_corr=optimal_correlation.get_figure()
      gvg.tight_layout()
      image_save = image_file_path + "_"+ dataset_counter + "_gvg_"+str(regulator)+"_"+str(target)+".png"
      gvg.savefig(image_save)
  target_roller.raw_data['Time']
  
  beforenorm = pd.read_csv(target_roller.file_path, sep='\t')
  beforenorm = beforenorm.dropna(axis=0,how='all')
  ax3 = beforenorm.plot(x='Time',y="G1",kind="scatter")
  fig3 = ax3.get_figure()
  image_save = image_file_path + "_time_series.png"
  fig3.savefig(image_save)

  ax2 = target_roller.raw_data.plot(x='Time', y="G1", kind="scatter")
  fig2 = ax2.get_figure()
  image_save = image_file_path + "_time_series_zscore.png"
  fig2.savefig(image_save)

  result_titles = ["precision","recall","aupr","tpr","fpr","auroc", "entropy"]

  result_list2 = [a_precision_list, a_recall_list, a_aupr_list, a_tpr_list, a_fpr_list,
      a_auroc_list]

  """for count,result in enumerate(result_list2):
      print(a_window_size_list)
      print(result)
      fig = plt.figure(6+count)
      plt.plot(a_window_size_list, result,'o',)
      title_string = "Window Size, lasso, Network " +img_suffix
      plt.title(title_string)
      plt.xlabel('Window Size')
      plt.ylabel(result_titles[count])
      image_save = image_file_path + "_windowsize_lasso_" +img_suffix+"_"+ str(result_titles[count]) + ".png"
      fig.savefig(image_save)

  print(window_size_list)
  """

  size_to_auroc = zip(window_size_list,window_start_list, window_auroc)
  f = plt.figure(figsize=(10,10))
  axarr2 = f.add_axes([gp_left, gp_bottom, gp_width, gp_height])
  axarr1 = f.add_axes([dm_left, dm_bottom, dm_width, dm_height])

  #figure, heatmap = Grapher.generate_heatmap_from_df(self.roller_results.raw_data)
  #axarr1 = figure

  raw_data= target_roller.raw_data
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
  #axarr1.grid("off")
  #axarr2.grid("off")
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
  axarr2.xaxis.set_ticks(np.arange(0,1050,50))
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
  image_save = image_file_path + "_windowsize_lasso_heatmap_" + dataset_counter + ".png"
  f.savefig(image_save,format="png")

