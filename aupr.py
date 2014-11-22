__author__ = 'jfinkle'

import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve
import sys
import matplotlib.pyplot as plt

# Load adjacency matrices

def create_link_list(df, w):
    parent_names = df.index.values
    child_names = df.columns.values
    parent_index = range(len(parent_names))
    child_index = range(len(child_names))
    a, b = np.meshgrid(parent_index, child_index)
    parents = parent_names[a.flatten()]
    children = child_names[b.flatten()]
    directed_edges = df.values.flatten()
    all_edges = np.abs(directed_edges)
    ll_array = np.array([parents, children, zip(parents, children), directed_edges, all_edges, weights])
    link_list = pd.DataFrame(ll_array).transpose()
    link_list.columns = ['Parent', 'Child', 'Edge', 'Directed_Edge', 'Edge_Exists', 'W']
    #link_list.sort(columns='Edge_Exists', ascending=False, inplace=True)
    return link_list

def calc_roc(ref, pred):
    # True Positive Rate (TPR) = TP/(TP+FN)
    # False Positive Rate (FPR) = FP/(FP+TN)
    ref.sort(columns='Edge', inplace=True)
    pred.sort(columns='Edge', inplace=True)
    if not np.array_equal(ref.Edge.values, pred.Edge.values):
        print 'Not same edges'
        return

    pred.sort(columns='W', ascending=False, inplace=True)
    ref_edge_list = ref.Edge.values.tolist()
    num_edges = len(ref_edge_list)

    tp = 0.0
    fp = 0.0
    fn = 0.0
    tn = 0.0
    precision = []
    recall = []

    for ii, row in enumerate(pred.iterrows()):
        pred_edge = row[1].Edge
        predicted = row[1].Edge_Exists
        ref_idx = ref_edge_list.index(pred_edge)
        real = ref.Edge_Exists.values[ref_idx]
        if real and predicted:
            tp+=1
            #print real, predicted, 'tp'
        elif not real and predicted:
            fp +=1
            #print real, predicted, 'fp'
        elif real and not predicted:
            fn +=1
            #print real, predicted, 'fn'
        elif not real and not predicted:
            tn+=1
            #print real, predicted, 'tn'
        cur_fn = num_edges-ii+1+fn
        if tp ==0 and cur_fn ==0:
            recall.append(0.0)
        else:
            recall.append(tp/(tp+cur_fn))
        if fp ==0 and tp ==0:
            precision.append(0.0)
        else:
            precision.append(tp/(tp+fp))

    return precision, recall

xls = pd.ExcelFile('goldbetter_model/adjacency_matrix.xlsx')
df = xls.parse()

xls2 = pd.ExcelFile('goldbetter_model/test_matrix.xlsx')
df2 = xls2.parse()

np.random.seed(8)
weights = np.random.random(len(df2)**2)
reference = create_link_list(df, weights)
random = np.array([np.sum(reference.Edge_Exists.values)/256.0]*256)
prediction = create_link_list(df2, weights)
p, r = calc_roc(reference, prediction)
plt.plot(r, p, r, random, 'r')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend(['Test', 'Random'])
plt.show()