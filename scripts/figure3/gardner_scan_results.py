import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
from palettable.colorbrewer.qualitative import Set1_3

df = pd.read_pickle('gardner_RF_scan_1000trees.pkl')
df['minmax'] = list(zip(df.min_lag, df.max_lag))
avg = df.groupby(['td_window', 'minmax']).mean().reset_index()
auroc = ((avg.pivot('minmax', 'td_window', 'auroc')-.756)/.756*100)
aupr = ((avg.pivot('minmax', 'td_window', 'aupr')-.286)/.286*100)

minmax = df.groupby('minmax')
pdata = minmax.get_group((0, 1))
pdata = pdata[pdata.td_window != 14]
baseline = df[(df.td_window==14) & (df.min_lag==0) & (df.max_lag == 0)]
baseline_mean = np.mean(baseline.aupr)
pdata = pd.concat([pdata, baseline], axis=0)
# Color selector
colors = []
for w in range(2, 14):
    test_data = pdata[pdata.td_window == w]
    diff = np.mean(test_data.aupr)-baseline_mean
    if stats.ttest_rel(test_data.aupr, baseline.aupr).pvalue < 0.05:
        if diff > 0:
            colors.append(Set1_3.mpl_colors[0])
        else:
            colors.append(Set1_3.mpl_colors[1])
    else:
        colors.append('grey')
# Append for baseline
colors.append('k')

sns.boxplot(data=pdata, x='td_window', y='aupr', palette=colors)
# plt.plot([-1, 12], [baseline_mean, baseline_mean], '-k', zorder=0)
plt.show()

#
# # print(auroc[10.0])
# ax = sns.heatmap(auroc)
#
# plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
# plt.show()
