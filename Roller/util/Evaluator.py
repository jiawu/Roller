__author__ = 'jfinkle'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def possible_edges(parents, children):
    """
    Create a list of all the possible edges between parents and children

    :param parents: array
        labels for parents
    :param children: array
        labels for children
    :return: array, length = parents * children
        array of parent, child combinations for all possible edges
    """
    parent_index = range(len(parents))
    child_index = range(len(children))
    a, b = np.meshgrid(parent_index, child_index)
    parent_list = parents[a.flatten()]
    child_list = children[b.flatten()]
    possible_edge_list = np.array(zip(parent_list, child_list))
    return possible_edge_list

def create_link_list(df, w):
    parent_names = df.index.values
    child_names = df.columns.values
    edges = possible_edges(parent_names, child_names)
    parents = edges[:, 0]
    children = edges[:, 1]
    directed_edges = df.values.flatten()
    all_edges = np.abs(directed_edges)
    ll_array = np.array([parents, children, zip(parents, children), directed_edges, all_edges, w])
    link_list = pd.DataFrame(ll_array).transpose()
    link_list.columns = ['Parent', 'Child', 'Edge', 'Directed_Edge', 'Edge_Exists', 'W']
    #link_list.sort(columns='Edge_Exists', ascending=False, inplace=True)
    return link_list

def calc_pr(ref, pred):
    # True Positive Rate (TPR) = TP/(TP+FN)
    # False Positive Rate (FPR) = FP/(FP+TN)
    ref.sort(columns='Edge', inplace=True)
    pred.sort(columns='Edge', inplace=True)
    if not np.array_equal(ref.Edge.values, pred.Edge.values):
        print 'Not same edges'
        return

    pred.sort(columns='W', ascending=False, inplace=True)
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