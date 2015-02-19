import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from sets import Set
import pandas as pd
import pdb

class Evaluator:

    def __init__(self, gold_standard_file, sep = '/t',interaction_label = 'regulator-target'):
        self.gs_file = gold_standard_file
        self.gs_data = pd.read_csv(gold_standard_file, sep=sep, header=None)
        self.gs_data.columns = ['regulator','target','exists']
        self.gs_data['regulator-target'] = zip(self.gs_data.regulator, self.gs_data.target)
        self.interaction_label = interaction_label
        self.gs_flat = self.gs_data[self.gs_data['exists'] > 0]['regulator-target']
        #more robust version of defining the full list
        all_regulators = self.gs_data['regulator'].unique().tolist()
        all_targets = self.gs_data['target'].unique().tolist()
        all_regulators.extend(all_targets)
        all_regulators = np.array(list(set(all_regulators)))
        
        self.full_list = tuple(map(tuple,self.possible_edges(all_regulators,
          all_regulators)))
        #remove self edges
        self.full_list = [ x for x in self.full_list if x[0] != x[1] ]
        self.full_list = pd.Series(self.full_list)

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

            aupr = integrate.cumtrapz(y=precision, x=recall)
        return precision, recall, aupr

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

    def calc_roc(self, pred):
        # True Positive Rate (TPR) = TP/(TP+FN)
        # False Positive Rate (FPR) = FP/(FP+TN)

        tp = 0.0
        fp = 0.0
        tpr = []
        fpr = []
        current_list = []
        for edge in pred['regulator-target']:
            current_list.append(edge)
            counts = self._evaluate(current_list)

            total_p = counts['tp']+ counts['fn']
            total_n = counts['fp']+ counts['tn']

            if total_n == 0:
                fpr.append(0.0)
            else:
                fpr.append(counts['fp']/total_n)
            if total_p == 0:
                tpr.append(0.0)
            else:
                tpr.append(counts['tp']/total_p)

            auroc = integrate.cumtrapz(x=fpr, y=tpr)
        return tpr, fpr, auroc

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
# Load adjacency matrices
