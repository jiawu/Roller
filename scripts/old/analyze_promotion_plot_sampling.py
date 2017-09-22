import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.append('../pipelines')
import pdb
import Pipelines as pl
import numpy as np
import pandas as pd

from Swing.util.Evaluator import Evaluator
import pickle

def parse(method_string):
    
    min_lag = 25
    max_lag = 50
    td_window = 10
    preinf_method = method_string.split('-')[0]
    inf_method = inf_method_key[preinf_method]

    misc = method_string.split('-')[1]

    if "td" in misc:
        td_window = int(misc.split('_')[1])
        if td_window == 21:
            min_lag = 0
            max_lag = 0

    elif "ml" in misc:
        case = int(misc.split('_')[1])
        if case == 0:
            min_lag = 0
            max_lag = 25
        elif case == 1:
            min_lag = 15
            max_lag = 25
        elif case == 2:
            min_lag = 25
            max_lag = 35

        elif case == 3:
            min_lag = 35
            max_lag = 45

        elif case == 4:
            min_lag = 45
            max_lag = 55

        elif case == 5:
            min_lag = 55
            max_lag = 65    
        elif case == 6:
            min_lag = 65
            max_lag = 75    
        elif case == 7:
            min_lag = 75
            max_lag = 85    
        elif case == 8:
            min_lag = 25
            max_lag = 45    
        elif case == 9:
            min_lag = 25
            max_lag = 70    
        elif case == 10:
            min_lag = 10
            max_lag = 50    

    return(inf_method, td_window, min_lag, max_lag)

td_edge_lists_lag = []

run_params = {}
inf_list = ["Dionesus", "RandomForest", "Lasso"]
inf_list2 = ["Dionesus-td10", "RandomForest-td10", "Lasso-td10", "Dionesus-td15", "RandomForest-td15", "Lasso-td15"]

methods_of_interest = ['Dionesus-td_2', 'Dionesus-td_10', 'Dionesus-td_15','Dionesus-td_21', 'Dionesus-ml_0', 'Dionesus-ml_1', 'Dionesus-ml_2','Dionesus-ml_3','Dionesus-ml_4', 'Dionesus-ml_5','Dionesus-ml_6','Dionesus-ml_7','Dionesus-ml_8','Dionesus-ml_9','Dionesus-ml_10','RF-td_2', 'RF-td_10', 'RF-td_15','RF-td_21', 'RF-ml_0', 'RF-ml_1', 'RF-ml_2','RF-ml_3','RF-ml_4', 'RF-ml_5','RF-ml_6','RF-ml_7','RF-ml_8','RF-ml_9','RF-ml_10', 'Lasso-td_10','Lasso-td_15']
 

inf_method_key = {'RF':'RandomForest', 'Dionesus': 'Dionesus', 'Lasso':'Lasso'}

organisms = ['high_sampling']

organism_results = []
for organism in organisms:    
    network_results = []
    for network_index in range(1,6):
        method_results = []
        for method in methods_of_interest:
        
            inf_method, td_window, min_lag, max_lag = parse(method)
            td_edge_lists_lag.append([])

            run_params['file_path'] = "/home/jjw036/Roller/data/dream4/"+str(organism)+"/insilico_size10_"+str(network_index)+"_timeseries.tsv"
            run_params['data_folder'] = "/projects/p20519/roller_output/sampling/" + str(inf_method) +"/" + organism + "-insilico_size10_" + str(network_index)
            run_params['td_window'] = td_window
            run_params['min_lag'] = min_lag
            run_params['max_lag'] = max_lag
            run_params['n_trees'] = 500
            run_params['bootstrap_n'] = 50
            run_params['filter_noisy'] = False
            print(run_params)
            
            roc, pr, tdr = pl.get_td_stats(**run_params)
            result_table = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
            result_table['rank_importance'] = np.arange(len(result_table))
            method_results.append(result_table)

        for i in range(0,len(methods_of_interest)):
            method_results[i] = method_results[i].rename(columns = {"rank_importance":"rank_importance_" + methods_of_interest[i]})
        current_results = method_results[0]
        for i, frame in enumerate(method_results[1:], 2):
            current_results = current_results.merge(frame, on='regulator-target')
        
        current_gold_standard = run_params['file_path'].replace("timeseries.tsv","goldstandard.tsv")
        node_list = ['G'+ str(x) for x in range(1,11)]
        evaluator = Evaluator(current_gold_standard, '\t', node_list = node_list)
        true_edges = evaluator.gs_flat.tolist()

        true_only = current_results[current_results['regulator-target'].isin(true_edges)]
        
        for method in methods_of_interest:
            diff_name = 'Dionesus-td_21-' + method
            true_only[diff_name] = true_only['rank_importance_Dionesus-td_21']-true_only['rank_importance_'+method]
        pickle_name = '/projects/p20519/roller_output/pickles/' + organism + '_net' + str(network_index) + '_promotion.pkl'
        pickle.dump(current_results, open(pickle_name, 'wb'))

        #network_results.append(current_results)
    #organism_results.append(network_results)
"""
##pickle.dump(organism_results, open("promotion_plot_results.p","wb"))

mm = pickle.load(open("promotion_plot_results.p","rb"))
mm = mm[0][0]

current_gold_standard = run_params['file_path'].replace("timeseries.tsv","goldstandard.tsv")
node_list = ['G'+ str(x) for x in range(1,11)]
evaluator = Evaluator(current_gold_standard, '\t', node_list = node_list)
true_edges = evaluator.gs_flat.tolist()

true_only = mm[mm['regulator-target'].isin(true_edges)]
true_only['diff'] = true_only['rank_importance_RF-td_21']-true_only['rank_importance_RF-td_15']

print(true_only)
print('promoted',len(true_only[true_only['diff']>0]))
print('demoted',len(true_only[true_only['diff']<0]))
f = plt.figure(figsize = (10,10) )
axes = f.gca()


for idx,row in true_only.iterrows():
    name = row['regulator-target']
    y1 = row['rank_importance_RF-td_21']
    y2 = row['rank_importance_RF-td_15']
    if row['diff'] > 0:
        color = 'blue'
    elif row['diff'] < 0:
        color = 'red'
    elif row['diff'] == 0:
        color = 'gray'
    y = [y1,y2]
    x = [1,2]
    axes.plot(x,y,color=color, label = name, linewidth = 2.0)

axes.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off


axes.invert_yaxis()    
f.savefig('ecoli-5-promotion-plot.ps', format = "ps", bbox_inches='tight')


# for each inference method, get results
# list of things to compare:
#  DIONESUS
#  LASSO
#  RF
#  td dionesus, lag = 1
#  td lasso
#  td RF
#  lag = 2
#  lag = 3
"""
