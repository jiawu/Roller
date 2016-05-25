from datetime import datetime
import pandas as pd
import os
import sys
import time

import pdb

import pickle

def get_results(input_folder):
    agg_df = pd.DataFrame()
    rank_data = {}
    print("Processing ", input_folder)
    print("There are %d files in this folder." % (len(os.listdir(input_folder))))
    for counter,file_path in enumerate(os.listdir(input_folder)):
        if os.stat(input_folder+file_path).st_size > 0:
            with open(input_folder+file_path, 'r') as f_input:
                param_dict = {}
                for line in f_input:
                    if line.startswith('!'):
                        (key,value) = line[1:].rstrip().split(',')
                        param_dict[key] = value
            rank_df = pd.read_csv(input_folder+file_path, sep='\t',comment='!')
            rank_data[input_folder+file_path] = rank_df
            
            param_dict['result_path'] = input_folder+file_path
            agg_df = agg_df.append(pd.Series(param_dict), ignore_index=True)
            if counter%1000 == 0:
                print("Done with %d out of %d files" % (counter,len(os.listdir(input_folder))))

    return(agg_df, rank_data)

def main(file_index):
    start_time = time.time()
    agg_df = pd.DataFrame()
    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    input_folder_list = ["/projects/p20519/roller_output/ranks/community/","/projects/p20519/roller_output/ranks/Lasso/","/projects/p20519/roller_output/ranks/Dionesus/","/projects/p20519/roller_output/ranks/RandomForest/"]
    input_folder_list = [input_folder_list[file_index-1]]
    for input_folder in input_folder_list:
        agg_df, rank_data = get_results(input_folder)
        substr = input_folder.split('/')[-2]
        agg_df.to_pickle(substr+'_agg_rank.pkl')
        with open(substr+'_rank_data.pkl', 'wb') as fb:
            pickle.dump(rank_data,fb)

        my_group = agg_df.groupby(['file_path', 'min_lag', 'max_lag', 'td_window'])

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    file_index = int(sys.argv[1])
    main(file_index)
