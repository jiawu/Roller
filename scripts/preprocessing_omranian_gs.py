import pandas as pd
import pdb
raw_gs1 = pd.read_csv('../data/invitro/regulon_tf_tf.tsv', sep = '\t')

raw_gs2 = pd.read_csv('../data/invitro/regulon_tf_gene.tsv', sep = '\t')

combined_gs = raw_gs1.append(raw_gs2)
combined_gs['parent'] = combined_gs['parent'].str.lower()
combined_gs['child'] = combined_gs['child'].str.lower()
# Only get strong connections
#combined_gs = combined_gs[combined_gs['strength'] == 'Strong']
combined_gs = combined_gs.drop_duplicates(['parent','child', 'reg'])

# Drop self-edges
combined_gs = combined_gs[combined_gs['parent'] != combined_gs['child']]

# Filter duplicated values

duplicated_values = combined_gs[combined_gs.duplicated(['parent','child'])]
for idx,row in duplicated_values.iterrows():
  current_set = combined_gs[(combined_gs['parent'] == row['parent']) & (combined_gs['child'] == row['child'])]
  if ('+' in current_set['reg'].unique()) and ('-' in current_set['reg'].unique()) or ('+-' in current_set['reg'].unique()):
    current_set['reg'] = '+-'
    combined_gs = combined_gs.append(current_set)

# Check the "?" values and change them to '+-'
question_values = combined_gs[combined_gs['reg'] == '?']
for idx,row in question_values.iterrows():
  current_set = combined_gs[(combined_gs['parent'] == row['parent']) & (combined_gs['child'] == row['child'])]
  current_set['reg'] = '+-'
  combined_gs = combined_gs.append(current_set)

no_dup = combined_gs.drop_duplicates(['parent','child'])

no_dup.to_csv('../data/invitro/regulon_parsed.tsv', sep='\t',index=False, header=False)

no_dup['edge'] = 1

final_gs = no_dup[['parent','child', 'edge']]
signed_final_gs = no_dup[['parent','child','reg']]

final_gs.to_csv('../data/invitro/omranian_parsed_goldstandard.tsv', sep='\t',index=False, header=False)
signed_final_gs.to_csv('../data/invitro/omranian_signed_parsed_goldstandard.tsv', sep='\t',index=False, header=False)


tf_list = no_dup['parent'].unique().tolist()
target_list = no_dup['child'].unique().tolist()

with open('../data/invitro/omranian_tf_list.tsv', 'w') as outfile:
    outfile.write("\n".join(tf_list))
with open('../data/invitro/omranian_target_list.tsv', 'w') as outfile:
    outfile.write("\n".join(target_list))

# there are 61 TFs that are not regulated by other TFs...
"""
raw_gs.to_csv('../data/invitro/omranian_parsed_goldstandard.tsv', sep='\t',index=False, header=False)

"""

