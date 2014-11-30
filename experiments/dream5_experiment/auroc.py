__author__ = 'jfinkle'

import pandas as pd
import numpy as np

# Load adjacency matrices

def create_link_list(df):
    parent_names = df.index.values
    child_names = df.columns.values
    parent_index = range(len(parent_names))
    child_index = range(len(child_names))
    a, b = np.meshgrid(parent_index, child_index)
    parents = parent_names[a.flatten()]
    children = child_names[b.flatten()]
    directed_edges = df.values.flatten()
    all_edges = np.abs(directed_edges)
    ll_array = np.array([parents, children, directed_edges, all_edges])
    link_list = pd.DataFrame(ll_array).transpose()
    link_list.columns = ['Parent', 'Child', 'Directed_Edge', 'Edge_Exists']
    link_list.sort(columns='Edge_Exists', ascending=False, inplace=True)
    return link_list

def calc_roc(ref, pred):
    # True Positive Rate (TPR) = TP/(TP+FN)
    # False Positive Rate (FPR) = FP/(FP+TN)

    if not np.array_equal(ref.Parent.values, pred.Parent.values) or not np.array_equal(ref.Child.values, pred.Child.values):
        print 'Not same edges'
        return

    tp = 0.0
    fp = 0.0
    fn = 0.0
    tn = 0.0
    tpr = []
    fpr = []
    real_edges = ref.Edge_Exists.values
    predicted_edges = prediction.Edge_Exists.values
    for real, predicted in zip(real_edges, predicted_edges):
        if real and predicted:
            tp+=1
        elif not real and predicted:
            fp +=1
        elif real and not predicted:
            fn +=1
        elif not real and not predicted:
            tn+=1

        if tp ==0 and fn ==0:
            tpr.append(0.0)
        else:
            tpr.append(tp/(tp+fn))
        if fp ==0 and tn ==0:
            fpr.append(0.0)
        else:
            fpr.append(fp/(fp+tn))

    return tpr, fpr

xls = pd.ExcelFile('goldbetter_model/adjacency_matrix.xlsx')
df = xls.parse()

reference = create_link_list(df)
prediction = create_link_list(df)
t_rate, f_rate = calc_roc(reference, prediction)
print t_rate
print f_rate