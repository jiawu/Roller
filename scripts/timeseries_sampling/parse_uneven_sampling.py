import pandas as pd
import pdb


fp = '/home/jjw036/Roller/data/gnw_insilico/high_sampling/Ecoli10'
current_file = '/Yeast-1_dream4_timeseries.tsv'
rate = 333

rates = [10,30,50,100,200,333,500]
uneven_rate = [0, 15, 30, 60, 120, 240, 480, 720]
# uneven sampling 0, 15, 30, 60, 120, 240, 480, 720

for network in range(1,21):
    current_file = '/Ecoli-{}_dream4_timeseries.tsv'.format(network)
    table = pd.read_csv(fp+current_file, sep='\t')
    cropped_table = table[table['Time'].isin(uneven_rate)]
    current_file = '/Ecoli-{}_uneven_timeseries.tsv'.format(network)
    cropped_table.to_csv(fp+current_file, header=True, index=False, sep='\t')

# regular sampling
#1, 10, 30, 50, 100, 200, 333, 500

# uneven sampling 0, 15, 30, 60, 120, 240, 480, 720
