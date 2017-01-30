import os
import re

import pdb
from shutil import copyfile

directory = "/home/jjw036/Roller/data/gnw_insilico/network_data/Yeast100/"
regexp = re.compile(r'Yeast100-\d+_timeseries.tsv')

for file in os.listdir(directory):
    if regexp.search(file) is not None:
        filename = directory+file
        with open(filename,'r') as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line)
        final_lines = lines[:211]
        nf = file.replace('Yeast100-','Yeast100-10p-')
        new_file_name = directory + nf
        with open(new_file_name,'w') as f_out:
            f_out.write( '\n'.join( final_lines ) )
        gs_file = file.replace('_timeseries.tsv','_goldstandard.tsv')
        new_gs_file = gs_file.replace('Yeast100-','Yeast100-10p-')
        copyfile(directory+gs_file, directory+new_gs_file)

