import pandas as pd
import pdb


fp = '/home/jjw036/Roller/data/gnw_insilico/high_sampling/Yeast10'
current_file = '/Yeast-1_dream4_timeseries.tsv'
rate = 333

rates = [10,30,50,100,200,333,500]

for network in range(1,21):
    current_file = '/Yeast-{}_dream4_timeseries.tsv'.format(network)
    table = pd.read_csv(fp+current_file, sep='\t')
    for rate in rates:
        interval = range(0,1001, rate)
        cropped_table = table[table['Time'].isin(interval)]
        current_file = '/Yeast-{}_{}_timeseries.tsv'.format(network,rate)
        cropped_table.to_csv(fp+current_file, header=True, index=False, sep='\t')

# regular sampling
#1, 10, 30, 50, 100, 200, 333, 500

# uneven sampling 0, 15, 30, 60, 120, 240, 480, 720
