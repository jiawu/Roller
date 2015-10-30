__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Swing.tdSwing import tdSwing
from Swing.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Loop through all 5 in silico networks
for ii in range(1,6):
        print ii
        print "=============================================="
        file_path = "../data/dream4/insilico_size10_%i_timeseries.tsv"%ii
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None
        current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
        evaluator = Evaluator(current_gold_standard, '\t')
        true_edges = evaluator.gs_flat.tolist()

        np.random.seed(8)

        tdr = tdSwing(file_path, gene_start_column, gene_end, time_label, separator)
        tdr.zscore_all_data()
        w_widths = range(2,21)
        auroc_list = []
        aupr_list = []
        save_name = file_path.split('/')[-1].split('.')[0]+'.csv'
        save_testing = 'mean_mean_lag_'
        lag_method = save_testing.replace('_lag_', '')
        save_path = '../output/tdRoller_testing/'+save_testing[:-1]+"/"
        for w in range(2,21):
            tdr.set_window(w)
            tdr.create_windows()
            tdr.augment_windows()
            tdr.fit_windows(n_trees=10, show_progress=False)
            tdr.compile_roller_edges(self_edges=True)
            #tdr.full_edge_list = tdr.full_edge_list[tdr.full_edge_list.Lag>1]
            tdr.make_static_edge_dict(true_edges, lag_method=lag_method)
            df2 = tdr.make_sort_df(tdr.edge_dict, 'lag')
            roc_dict, pr_dict = tdr.score(df2)
            auroc_list.append(roc_dict['auroc'][-1]+(1-roc_dict['fpr'][-1]))
            aupr_list.append(pr_dict['aupr'][-1]+(1-pr_dict['recall'][-1]))
        result_df = pd.DataFrame([w_widths, auroc_list, aupr_list], index=['W', 'AUROC', 'AUPR']).T
        result_df.to_csv(save_path+save_testing+save_name)
        plt.plot(result_df.W, result_df.AUROC, '.-', result_df.W, result_df.AUPR, '.-')
        plt.legend(['AUROC', 'AUPR'], loc='best')
        plt.xlabel('Window Size')
        plt.ylabel('Score')
        plt.yticks(np.linspace(0,1,11))
        save_name = save_name.replace('csv', 'png')
        plt.savefig(save_path+save_testing+save_name)
        plt.close()