import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
pd.set_option('display.width', 2000)
# Results from sensitivity analysis
data = pd.read_pickle('../../output/sensitivity_analysis/gardner_downsample_results_2017-10-19.pkl') # type: pd.DataFrame
grouped = data.groupby(['method', 'min_lag', 'max_lag', 'td_window'])
idx = pd.IndexSlice
avg = grouped.mean().loc[:, ['aupr', 'auroc']].stack().reset_index()
print(avg)
sys.exit()
avg = avg[((avg.td_window==5) & (avg.min_lag==1) & (avg.max_lag==3))| (avg.td_window==21)]
avg.columns = avg.columns[:-2].tolist()+['score', 'value']
g = sns.factorplot(x='td_window', y='value', data=avg, col='score', row='model', hue='method', kind='box')
plt.show()
sys.exit()