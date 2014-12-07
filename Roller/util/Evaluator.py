import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from sets import Set

import pdb

class Evaluator:

    def __init__(self, gold_standard_file, sep = '/t',interaction_label = 'regulator-target'):
        self.gs_file = gold_standard_file
        self.gs_data = pd.read_csv(gold_standard_file, sep=sep, header=None)
        self.gs_data.columns = ['regulator','target','exists']
        self.gs_data['regulator-target'] = zip(self.gs_data.regulator, self.gs_data.target)
        self.interaction_label = interaction_label
        self.gs_flat = self.gs_data[self.gs_data['exists'] > 0]['regulator-target']
        self.full_list = self.gs_data['regulator-target']

    def possible_edges(self,parents, children):
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

    def create_link_list(self,df, w):
        parent_names = df.index.values
        child_names = df.columns.values
        edges = self.possible_edges(parent_names, child_names)
        parents = edges[:, 0]
        children = edges[:, 1]
        directed_edges = df.values.flatten()
        all_edges = np.abs(directed_edges)
        ll_array = [parents, children, zip(parents, children), directed_edges, all_edges, w]
        link_list = pd.DataFrame(ll_array).transpose()
        link_list.columns = ['Parent', 'Child', 'Edge', 'Directed_Edge', 'Edge_Exists', 'W']
        #link_list.sort(columns='Edge_Exists', ascending=False, inplace=True)
        return link_list

    def calc_pr(self, pred, sort_on = 'exists'):
        # True Positive Rate (TPR) = TP/(TP+FN)
        # False Positive Rate (FPR) = FP/(FP+TN)
        ref = self.gs_data

        #ref.sort(columns=self.interaction_label, inplace=True)
        #pred.sort(columns=self.interaction_label, inplace=True)
        #if not np.array_equal(ref[self.interaction_label].values, pred[self.interaction_label].values):
            #print 'Not same edges'
            #return

        #pred.sort(columns=sort_on, ascending=False, inplace=True)

        ref_edge_list = ref[self.interaction_label].values.tolist()

        #initialize counts
        #todo: remove below, it is unnecessary
        counts = {}
        counts['tp'] = 0.0
        counts['fp'] = 0.0
        counts['fn'] = 0.0
        counts['tn'] = 0.0
        precision = []
        recall = []
        current_list = []

        #current_list is the list at a certain index.
        #at each iteration, add the edge to the current_list.
        #evaluate the current list's precision and recall value

        #it is assumed that the edges are already sorted in descending rank order

        for edge in pred['regulator-target']:
            current_list.append(edge)
            counts = self._evaluate(current_list)

            if counts['tp'] ==0 and counts['fn'] ==0:
                recall.append(0.0)
            else:
                recall.append(counts['tp']/(counts['tp']+counts['fn']))
            if counts['fp'] ==0 and counts['tp'] ==0:
                precision.append(0.0)
            else:
                precision.append(counts['tp']/(counts['tp']+counts['fp']))

            aupr = integrate.cumtrapz(precision, recall)
        return precision, recall, aupr[-1]

    def _evaluate(self, current_list):
        """ evaluate the list using sets hooray, packs it up in dict """
        gold_standard = Set(self.gs_flat)
        prediction = Set(current_list)
        full = Set(self.full_list)

        tp = len(prediction.intersection(gold_standard))
        #both in prediction and gold_standard
        fp = len(prediction.difference(gold_standard))
        #in prediction but not in gold standard
        pred_full = full.difference(prediction) # 'negatives' or in full but not pred
        gold_full = full.difference(gold_standard) # in full but not in gold

        tn = len(pred_full.intersection(gold_full))
        fn = len(pred_full.difference(gold_full))
        #compare pred - full negatives, and gold - full negatives
        #true negatives is the intersection between gold_full and pred full
        #false negatives is the number of things in pred full that are not in gold full
        results = {}
        results['tn'] = float(tn)
        results['tp'] = float(tp)
        results['fn'] = float(fn)
        results['fp'] = float(fp)
        return(results)

    def calc_roc(self,ref, pred):
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

        total_p = float(np.sum(ref.Edge_Exists.values))
        total_n = len(ref.Edge_Exists.values) - total_p

        tp = 0.0
        fp = 0.0
        tpr = []
        fpr = []

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

            tpr.append(tp/total_p)
            fpr.append(fp/total_n)

            auroc = integrate.cumtrapz(fpr, tpr)
        return tpr, fpr, auroc[-1]


if __name__ == '__main__':
    xls = pd.ExcelFile('../../goldbetter_model/adjacency_matrix.xlsx')
    df = xls.parse()

    xls2 = pd.ExcelFile('../../goldbetter_model/test_matrix.xlsx')
    df2 = xls2.parse()

    np.random.seed(8)
    weights = np.random.random(len(df2)**2)
    reference = create_link_list(df, weights)
    random = np.array([np.sum(reference.Edge_Exists.values)/256.0]*256)
    prediction = create_link_list(df2, weights)
    p, r, area = calc_pr(reference, prediction)
    plt.plot(r, p, r, random, 'r')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(['Test', 'Random'])
    print area
    plt.show()

    tpr, fpr, area = self.calc_roc(reference, prediction)
    plt.plot(fpr, tpr, fpr, fpr, 'r')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.legend(['Test', 'Random'], loc='best')
    print area
    plt.show()
