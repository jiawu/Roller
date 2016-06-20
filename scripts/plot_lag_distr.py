import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdb
from matplotlib.table import Table

def stext(ax,text):
    ax.text(0.9, 0.9, text, horizontalalignment='right', verticalalignment='top',transform=ax.transAxes)
  
def main():
    data = pd.DataFrame()
    # overall lag distribution, with total number of edges + filename
    #fig2 = plt.figure(2, figsize=(20,20))
    # lag counts vs # of edges above certain threshold, percentage, total number of edges
         
    for ind in range(2,5):
        fig1 = plt.figure(figsize=(20,20))
        filename = 'lag_df2_parse_biocyc_%s.pkl' % (ind)
        lag_df = pd.read_pickle(filename)

        lag_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in lag_df['Lag'].tolist()]

        ax1 = fig1.add_subplot(2,2,1)
        bins = [x for x in range(0,30)]
        
        ax1.hist(lag_df['lag_counts'], normed=True, bins = bins, color='salmon', alpha=0.8)
        stext(ax1, 'Counts')      
        ax2 = fig1.add_subplot(2,2,2)
        ax2.hist(lag_df['lag_std'].dropna(), normed=True, color='salmon', alpha=0.8)
        stext(ax2, 'Stdev')      
        ax3 = fig1.add_subplot(2,2,3)
        ax3.hist(lag_df['lag_mean'].dropna(), normed=True, color='salmon', alpha=0.8)
        stext(ax3, 'Mean')      
        ax4 = fig1.add_subplot(2,2,4)
        ax4.hist(lag_df['lag_median'].dropna(), normed=True, color='salmon', alpha=0.8)
        stext(ax4, 'Median')      

        plt.savefig('lag_stats_biocyc_%s.png'%(ind))

        fig2 = plt.figure(figsize=(20,20))
        # get percentage of edges and number over thresholds

        lag_thresh = [0,1,2]
        count_thresh = [0,5,10,15,20,25]

        data = np.zeros([len(lag_thresh),len(count_thresh)])
        for idx,lag in enumerate(lag_thresh):
            for jdx,count in enumerate(count_thresh):
                items = lag_df[(lag_df['lag_counts'] >= count) & (lag_df['lag_mean'] >=lag)]
                data[idx][jdx] = len(items)

        p_data = data/len(lag_df)
        data = pd.DataFrame(data, index = lag_thresh, columns = count_thresh)
        p_data = pd.DataFrame(p_data, index = lag_thresh, columns = count_thresh)
        checkerboard_table(data)
        plt.savefig('lag_table_biocyc_%s.png'%(ind))
        checkerboard_table(p_data)
        plt.savefig('lag_table_p_biocyc_%s.png'%(ind))

        # get the distribution of lag counts (histogram)
        # the total number of edges

        # lag counts vs # of edges above certain threshold


def checkerboard_table(data, fmt='{:.2f}', bkg_colors=['yellow', 'white']):
    fig, ax = plt.subplots()
    ax.set_axis_off()
    tb = Table(ax, bbox=[0,0,1,1])

    nrows, ncols = data.shape
    width, height = 1.0 / ncols, 1.0 / nrows

    # Add cells
    for (i,j), val in np.ndenumerate(data):
        # Index either the first or second item of bkg_colors based on
        # a checker board pattern
        idx = [j % 2, (j + 1) % 2][i % 2]
        color = bkg_colors[idx]

        tb.add_cell(i, j, width, height, text=fmt.format(val), 
                    loc='center', facecolor=color)

    # Row Labels...
    for i, label in enumerate(data.index):
        tb.add_cell(i, -1, width, height, text=label, loc='right', 
                    edgecolor='none', facecolor='none')
    # Column Labels...
    for j, label in enumerate(data.columns):
        tb.add_cell(-1, j, width, height/2, text=label, loc='center', 
                           edgecolor='none', facecolor='none')
    ax.add_table(tb)
    plt.ylabel('Lag Threshold')
    plt.xlabel('Count Threshold')
    return fig

if __name__ == '__main__':
    main()


