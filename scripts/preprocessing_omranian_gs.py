import pandas as pd
import pdb
raw_gs = pd.read_csv('../data/invitro/omranian_goldstandard.tsv', sep = '\t')
raw_gs['effect'] = 1

tf_list = raw_gs['TF'].unique().tolist()
target_list = raw_gs['gene'].unique().tolist()

with open('../data/invitro/omranian_tf_list.tsv', 'w') as outfile:
    outfile.write("\n".join(tf_list))
with open('../data/invitro/omranian_target_list.tsv', 'w') as outfile:
    outfile.write("\n".join(target_list))

# there are 61 TFs that are not regulated by other TFs...

raw_gs.to_csv('../data/invitro/omranian_parsed_goldstandard.tsv', sep='\t',index=False, header=False)



