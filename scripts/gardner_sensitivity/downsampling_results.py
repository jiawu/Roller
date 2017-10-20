import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel
import sys
pd.set_option('display.width', 2000)
# Results from sensitivity analysis
data = pd.read_pickle('../../output/sensitivity_analysis/gardner_downsample_results_2017-10-19.pkl') # type: pd.DataFrame
results_by_dataset = data.groupby('file_path')
for g, df in results_by_dataset:
    grouped = df.groupby(['method', 'min_lag', 'max_lag', 'td_window'])
    print(grouped.mean())
sys.exit()

print(grouped.mean())
base = grouped.get_group(('RandomForest', 0, 0, 7))
swing = grouped.get_group(('RandomForest', 0, 0, 7))
# print(grouped.get_group(('Lasso', 'ecoli', '13', 1, 3, 15)).mean())
    # print(grouped.get_group(('Lasso', 'ecoli', '13', 0, 0, 21)).mean())
idx = pd.IndexSlice
print(grouped.mean())
avg = grouped.mean().loc[:, ['aupr', 'auroc']].stack().reset_index()
print(avg)
sys.exit()
avg = avg[((avg.td_window==5) & (avg.min_lag==1) & (avg.max_lag==3))| (avg.td_window==21)]
avg.columns = avg.columns[:-2].tolist()+['score', 'value']
g = sns.factorplot(x='td_window', y='value', data=avg, col='score', row='model', hue='method', kind='box')
plt.show()
sys.exit()