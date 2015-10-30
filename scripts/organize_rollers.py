import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import pandas as pd
from Swing.util.Evaluator import Evaluator
import shutil
import sys

if __name__ == "__main__":
  file_start = int(sys.argv[1])
  #get all pickle files
  path="/projects/p20519/Roller_outputs/"
  filenames = next(os.walk(path))[2]
  nfiles = len(filenames)
  #organize pickled objects by dataset analyzed
  obj_list = []
  counter = 0
  image_file_path = "/home/jjw036/Swing/aggregated"

  target_dataset =  "data/dream4/insilico_size10_1_timeseries.tsv"

  target_directory = "/projects/p20519/sorted_Rollers/"


  #currently it's really hard to unpickle all these objects at once. This script
  #tries to mv files into different directories according to dataset.

  file_end = file_start + 170

  if file_end > nfiles:
    file_end = nfiles -1

  for file in filenames[file_start:file_end]:
    full_path = path + file
    print(full_path)
    roller_obj = pd.read_pickle(full_path)
    attributes = dir(roller_obj)
    if any("file_path" in attribute for attribute in attributes):
      counter += 1
      key = roller_obj.file_path
      print(str(counter) + " out of " + str(nfiles))
      directory = target_directory + key
      directory = directory[:-4]
      if not os.path.exists(directory):
        os.makedirs(directory)
      shutil.copy(full_path, directory)
