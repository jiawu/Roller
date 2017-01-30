import matplotlib
matplotlib.use('Agg')
from Swing.util.BoxPlot import BoxPlot
from matplotlib.backends.backend_pdf import PdfPages

import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

"""
Script that loads data from a dataframe and generates boxplots

"""
input_file = "cluster_summary_all_c26_2016-06-07_15-46-49.tsv"

df = pd.read_csv(input_file, sep='\t')
df['diff_auroc'] = df['swing_auroc'] - df['baseline_auroc']
df['diff_aupr'] = df['swing_aupr'] - df['baseline_aupr']
lagged_modules = df[df['percent_lagged'] > 0.1]
nlagged_modules = df[df['percent_lagged'] < 0.1]

lagged_b_auroc = lagged_modules['baseline_auroc'].tolist()
nlagged_b_auroc = nlagged_modules['baseline_auroc'].tolist()

lagged_s_auroc = lagged_modules['swing_auroc'].tolist()
nlagged_s_auroc = nlagged_modules['swing_auroc'].tolist()

lagged_d_auroc = lagged_modules['diff_auroc'].tolist()
nlagged_d_auroc = nlagged_modules['diff_auroc'].tolist()

lagged_b_aupr = lagged_modules['baseline_aupr'].tolist()
nlagged_b_aupr = nlagged_modules['baseline_aupr'].tolist()

lagged_s_aupr = lagged_modules['swing_aupr'].tolist()
nlagged_s_aupr = nlagged_modules['swing_aupr'].tolist()

lagged_d_aupr = lagged_modules['diff_aupr'].tolist()
nlagged_d_aupr = nlagged_modules['diff_aupr'].tolist()

with PdfPages('regulon_db_network_comparison.pdf') as pdf:

  bp = BoxPlot()

  labels = ['NTD','NTD','TD','TD']
  bp.plot_box([nlagged_b_auroc, nlagged_s_auroc, lagged_b_auroc, lagged_s_auroc], labels)
  bp.add_formatting('AUROC', y_label='AUROC')        
  pdf.savefig(bp.f)

  bp = BoxPlot()
  
  bp.plot_box([nlagged_b_aupr,nlagged_s_aupr, lagged_b_aupr, lagged_s_aupr], labels)
  bp.add_formatting('AUPR', y_label='AUPR')        
  pdf.savefig(bp.f)



