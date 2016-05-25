import pandas as pd
import pdb

# Run the following scripts in order:

# preprocessing_omranian_gs.py
# preprocessing_omranian.py
# preprocessing_omranian_gs_part2.py

gs = pd.read_csv('../data/invitro/omranian_parsed_goldstandard.tsv', sep = '\t', header = None)
signed_gs = pd.read_csv('../data/invitro/omranian_signed_parsed_goldstandard.tsv', sep = '\t', header= None)


# remove df in signed_gs not in gs
edges = tuple(zip(gs[0],gs[1]))
signed_gs['edge'] = tuple(zip(signed_gs[0],signed_gs[1]))

new_signed_gs = signed_gs[signed_gs['edge'].isin(edges)]
new_signed_gs = new_signed_gs[[0,1,2]] 
new_signed_gs.to_csv('../data/invitro/omranian_signed_parsed_goldstandard.tsv', sep='\t',index=False, header=False)
