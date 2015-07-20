__author__ = 'jfinkle'

"""
This is a basic script to produce figures as requested by Dr. Bagheri.

Hopefully it produces some reusable code.
"""

import sys
import itertools
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd

# Read in the data
filepath = '../output/tdRoller_testing/20150719_roller_tdroller_comparison.xlsx'
aupr = pd.read_excel(filepath, 'AUPR', skip_footer=2) # Ignore last two lines with comments
aupr_std = pd.Series([np.std(row, ddof=1) for row in aupr.values.T], index=aupr.columns)
aupr_mean = pd.Series([np.mean(row) for row in aupr.values.T], index=aupr.columns)
aupr_std.name='std'

auroc = pd.read_excel(filepath, 'AUROC', skip_footer=2)
auroc_std = pd.Series([np.std(row, ddof=1) for row in auroc.values.T], index=auroc.columns)
auroc_mean = pd.Series([np.mean(row) for row in auroc.values.T], index=aupr.columns)
auroc_std.name='std'


idx = ['All_Data', 'Roller', 'tdRFR_optimized']

temp_mean = auroc_mean[idx].values
temp_std = auroc_std[idx].values

width = 0.9

# Pull the formatting out here
bar_kwargs = {'width':width,'color':'c','zorder':5, 'align':'center'}
err_kwargs = {'zorder':0,'fmt':'none','lw':2,'ecolor':'k'}

idx = np.arange(len(temp_mean))
plt.bar(idx, temp_mean,  **bar_kwargs)
plt.errs = plt.errorbar(idx, temp_mean, yerr=temp_std, **err_kwargs)
plt.show()


sys.exit()
# This stuff didn't work out so well
ttest_vals=[1, 1, 1]
ttest_vals[0] = round(stats.ttest_ind(auroc['All_Data'], auroc['Roller'])[1], 4)
ttest_vals[1] = round(stats.ttest_ind(auroc['All_Data'], auroc['tdRFR_optimized'])[1], 4)
ttest_vals[2] = round(stats.ttest_ind(auroc['Roller'], auroc['tdRFR_optimized'])[1], 4)
ind  = np.arange(len(idx))    # the x locations for the groups
width= 0.85
labels = (['All Data', 'Roller', 'tdRFR'])

# Pull the formatting out here
bar_kwargs = {'width':width,'color':'c','zorder':5}
err_kwargs = {'zorder':0,'fmt':'none','lw':2,'ecolor':'k'}

X = ind+width/2
print X

fig, ax = plt.subplots()
ax.p1 = plt.bar(ind, temp_mean, **bar_kwargs)
ax.errs = plt.errorbar(X, temp_mean, yerr=temp_std, **err_kwargs)


# Custom function to draw the diff bars

def label_diff(i,j,text,X,Y):
    x = (X[i]+X[j])/2
    y = 1.1*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'<->', 'shrinkA':40,'shrinkB':20,'lw':1}
    print text, X[i], X[j], y

    ax.annotate(text, xy=(X[i],y), xytext=(x,y+0.1), zorder=10, ha='center')
    ax.annotate('', xy=(X[i],y-0.1), xytext=(X[j],y-0.1), arrowprops=props, ha='center')

# Call the function
diff1 = 'p = %0.4f' %ttest_vals[0]
diff2 = 'p = %0.4f' %ttest_vals[1]
diff3 = 'p = %0.4f' %ttest_vals[2]
label_diff(0,1,diff1,X,temp_mean)
label_diff(0,2,diff2,X,temp_mean)
#label_diff(1,2,diff3,X,temp_mean)



plt.ylim(ymax=1.2)
plt.xticks(X, labels, color='k')
plt.show()
