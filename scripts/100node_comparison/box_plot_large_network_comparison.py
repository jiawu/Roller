import matplotlib
import seaborn as sns

# from Swing.util.BoxPlot import BoxPlot
# from matplotlib.backends.backend_pdf import PdfPages

import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from scipy.stats import mannwhitneyu

pd.set_option('display.width', 2000)
"""
Script that loads data from a dataframe and generates boxplots

"""


def read_tdr_results(folder_list):
    agg_df = pd.DataFrame()
    for input_folder in folder_list:
        for file_path in os.listdir(input_folder):
            if ".tsv" in file_path:
                # print(file_path)
                df = pd.read_csv(input_folder + file_path, sep=',|\t', engine='python')

            # Correct weird inputs:
            # Omranian data has a messed up header
            if 'omranian' in file_path:
                new_cols = df.columns[1:]
                # Delete last column
                del df[df.columns[-1]]
                df.columns = new_cols

            agg_df = agg_df.append(df)
    return (agg_df)


def parse_tdr_results(agg_df, test_statistic, datasets):
    label_list = []
    auroc_list = []

    ## Analyze:
    # nonuniform
    # uniform
    # for all networks 1 2 3 4 5
    # parsing for windows = 7, windows = 4

    for dataset in datasets:
        current_df = agg_df[agg_df['file_path'].str.contains(dataset)]
        RF = current_df[(current_df['td_window'] == 21)]
        SWING_RF = current_df[(current_df['td_window'] == 15)]

        comparisons = [RF, SWING_RF]

        for category in comparisons:
            auroc_list.append(category[test_statistic][0:n_trials].tolist())
        label_list.append("Dionesus")
        label_list.append("SWING Dionesus")

        # label_list.append("Dionesus")
        # label_list.append("SWING Dionesus")

    return ((label_list, auroc_list))


output_path = "./"
methods = ['Dionesus', 'Lasso', 'RandomForest']

input_folder_list = ["/Users/jfinkle/Downloads/RandomForest/"]
test_statistic = ['aupr', 'auroc']
save_tag = "Dionesus_Yeast100_11-20"
n_trials = 100

datasets = ["Yeast100-" + str(index) + "_" for index in range(1, 21)]
# datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']
agg_df = read_tdr_results(input_folder_list)
g = agg_df.groupby('file_path')
net_stats = []
for name, group in g:
    # Quick and dirty parsing
    if 'Yeast' in name or 'Ecoli' in name:
        net_name = os.path.basename(name).split('_')[0]
    else:
        net_name = os.path.basename(name).split('_t')[0]
    system = net_name.split('100')[0].split('_')[0]
    subgroup = group.groupby('td_window')
    sg_mean = subgroup.mean()
    sg_std = subgroup.std()
    try:
        current_stats = [net_name, system, sg_mean.loc[21, 'aupr'], sg_std.loc[21, 'aupr'],
                         sg_mean.loc[15, 'aupr'], sg_std.loc[15, 'aupr'],
                         sg_mean.loc[21, 'auroc'], sg_std.loc[21, 'auroc'],
                         sg_mean.loc[15, 'auroc'], sg_std.loc[15, 'auroc']]
        net_stats.append(current_stats)
    except:
        pass

net_stats = pd.DataFrame(net_stats, columns=['network', 'system', 'base_aupr_mean', 'base_aupr_std', 'swing_aupr_mean',
                                             'swing_aupr_std', 'base_auroc_mean', 'base_auroc_std', 'swing_auroc_mean',
                                             'swing_auroc_std'])
net_stats['aupr_diff'] = (net_stats['swing_aupr_mean']-net_stats['base_aupr_mean'])/net_stats['base_aupr_mean']*100
net_stats['auroc_diff'] = (net_stats['swing_auroc_mean']-net_stats['base_auroc_mean'])/net_stats['base_auroc_mean']*100

sns.boxplot(data=net_stats, x='system', y='auroc_diff')
plt.ylabel('AUROC Increase (%)')
plt.savefig('./RF_100_node_auroc_increase.pdf', format='pdf')
plt.close()

sns.boxplot(data=net_stats, x='system', y='aupr_diff')
plt.ylabel('AUPR Increase (%)')
plt.savefig('./RF_100_node_aupr_increase.pdf', format='pdf')
for n, g in net_stats.groupby('system'):
    print(n, mannwhitneyu(g.base_aupr_mean, g.swing_aupr_mean), mannwhitneyu(g.base_auroc_mean, g.swing_auroc_mean))
sys.exit()

with PdfPages(output_path + save_tag + '.pdf') as pdf:
    for test in test_statistic:
        label_list, auroc_list = parse_tdr_results(agg_df, test, datasets)
        bp_data = auroc_list
        bp = BoxPlot()

        bp.plot_box(bp_data, label_list)
        title = save_tag
        bp.add_formatting(title, y_label=test.upper())
        pdf.savefig(bp.f)




        # auroc_1 = df['auroc'].values
        # auroc_2 = df['auroc'].values

        # bp_data = [auroc_1,auroc_2]

        # bp = BoxPlot()

        # bp.plot_box(bp_data, ['n_trees = 10', 'n_trees = 20'])


        # bp.save_plot(output_path, save_tag)



        # grouped.get_group((2,2)).mean()['aupr']
