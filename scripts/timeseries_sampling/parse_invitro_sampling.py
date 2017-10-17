import pandas as pd
import pdb


fp = '/home/jjw036/Roller/data/invitro'
current_file = '/gardner_timeseries.tsv'
rate = 333

rates = [2, 3, 4]

table = pd.read_csv(fp+current_file, sep='\t')
for rate in rates:
    interval = range(1,14, rate)
    cropped_table = table[table['Time'].isin(interval)]
    current_file = '/gardner_{}_timeseries.tsv'.format(rate)
    cropped_table.to_csv(fp+current_file, header=True, index=False, sep='\t')

# regular sampling
#1, 10, 30, 50, 100, 200, 333, 500

# uneven sampling 0, 15, 30, 60, 120, 240, 480, 720
