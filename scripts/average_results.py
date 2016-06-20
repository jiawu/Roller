
# get a list of DFs by reading the csv, append

# average list of DFs based off of 'regulator-target'

big_df = pd.concat(parsed_df_list)
mean_df=big_df.groupby(level=0).mean()


