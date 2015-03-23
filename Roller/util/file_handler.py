__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import os
import sys
import pandas as pd
import numpy as np


def load_roller_pickles(pickle_list):
    #organize pickled objects by dataset analyzed
    obj_list = []
    counter = 0

    dataset_dict = {}

    for file in pickle_list:
      full_path = path + file
      print(full_path)
      roller_obj = pd.read_pickle(full_path)
      attributes = dir(roller_obj)
      if any("file_path" in attribute for attribute in attributes):
        counter += 1
        print(str(counter) + " out of " + str(nfiles))
        obj_list.append(roller_obj)

        key = roller_obj.file_path
        if key not in dataset_dict:
          dataset_dict[key] = []

        dataset_dict[key].append(roller_obj)

    print(len(obj_list))
    print(dataset_dict.keys())

if __name__ == '__main__':
    path="../../output/Roller_outputs_RF_moretrees/"
    filenames = next(os.walk(path))[2]

    nfiles = len(filenames)
    load_roller_pickles(filenames)