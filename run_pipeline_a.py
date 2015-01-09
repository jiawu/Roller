from __future__ import absolute_import
import Roller
from sklearn.preprocessing import Imputer
from Roller.util.linear_wrapper import LassoWrapper
import numpy as np
import matplotlib as mpl
from Roller.util.permutation_test import Permuter
import itertools
from Roller.util import utility_module as utility
from Roller.util.Ranker import Bootstrapper
from Roller.util.Evaluator import Evaluator
import pdb
import pickle
import scipy
import sys
import getopt

def main(argv):
    parameter_file = ''

    ### main roller parameters ###
    input_file = ''
    gene_start_column = ''
    time_label = ''
    separator = "\t"
    gene_end = None
    window_size = ''
    alpha = ''

    ### bootstrapping parameters ###
    boots = ''
    max_random = ''
    n_alphas = ''

    ### saving bootstrapped and permuted results ###
    load_file = ''
    save_file = ''

    ### AUPR calculations ###
    gold_standard_file = ''
    rank_by = ''

    try:
        opts, args = getopt.getopt(argv,"hi:w:a:",["ifile=","window_size","alpha"])
    except getopt.GetoptError:
        print 'run_pipeline_a.py -i <parameter-file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'run_pipeline_a.py -i <parameter-file>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            parameter_file = arg
        elif opt in ("-w", "--window_size"):
            window_size = arg
        elif opt in ("-a", "--alpha"):
            alpha = arg
    print 'The parameter input file is "', parameter_file
    pd = setup_pipeline(parameter_file, window_size, alpha)
    results = execute(pd)
    print(results['aupr'])

def setup_pipeline(parameter_file, window_size, alpha):
    ### read parameters from file ###
    pd = readParams(parameter_file)

    ### override file pds with command call -- useful for looping through windows ###
    if pd['window_size'] == '':
        pd['window_size'] = window_size
    if pd['alpha'] == '':
        pd['alpha'] = float(alpha)
    return(pd)

def execute(pd):

    ### part 1, create inital model ###
    initial_model, roller, pd = initialize_model(pd)
    roller.reset()

    existing_data = False

    ## part 1a, save model ###
    if pd['load_file'] != "None":
        loaded_file = pickle.load(open(pd['load_file'],"rb"))
        existing_data = True
        all_panels = loaded_file

    if existing_data == False:
        ### part 2, permutation tests ###
        permuted_model_means, permuted_model_sd = permute_model(roller, pd)
        roller.reset()
        ### part 3, bootstrapping tests ###
        stability_model = bootstrap_model(roller, pd)

        ### part 4, aggregate data ###
        all_panels = [initial_model,permuted_model_means, permuted_model_sd, stability_model]

    results_table = aggregate_model(all_panels, pd)

    ### part 5, rank edges ###
    ranked_table = rank_model(results_table, pd)

    ### part 6, calculate AUPR ###
    precision, recall, aupr = score_model(ranked_table, pd)

    ### package results into a dict ###
    final_result = { 'precision': precision, 'recall':recall, 'aupr': aupr, 'ranked_table':ranked_table, 'results_table': results_table, 'param_dict': pd}
    if existing_data == False:
        pickle.dump(final_result, open(pd["save_file"], "wb"))
    return(final_result)

def readParams(parameter_file):
    param_dict = {}
    with open(parameter_file, 'r') as p_file:
        for line in p_file:
            splitLine = line.strip().split('=')
            param_dict[str(splitLine[0])] = splitLine[1]
    #todo: test if all required parameters are specified
    #further parse parameter types
    param_dict['gene_start_column'] = int(param_dict['gene_start_column'])
    param_dict['window_size'] = int(param_dict['window_size'])
    param_dict['boots'] = int(param_dict['boots'])
    param_dict['n_alphas'] = int(param_dict['n_alphas'])
    param_dict['max_random'] = float(param_dict['max_random'])

    param_dict['alpha'] = float(param_dict['alpha'])
    if param_dict['gene_end'] == "None":
        param_dict['gene_end'] = None
    if param_dict['gold_standard_sep'] == '\\t':
        param_dict['gold_standard_sep'] = '\t'
    if param_dict['separator'] == '\\t':
        param_dict['separator'] = '\t'
    #param_dict['separator'] = unicode(param_dict['separator'], 'unicode_escape')
    return(param_dict)

def initialize_model(pd):
    file_path = pd['file_path']
    gene_start_column = pd['gene_start_column']
    time_label = pd['time_label']
    separator = pd['separator']
    gene_end = pd['gene_end']
    window_size = pd['window_size']
    alpha = pd['alpha']

    roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
    roll_me.set_window(window_size)
    total_window_number = roll_me.get_n_windows()

    imputer = Imputer(missing_values="NaN")
    mpl.rcParams.update({'font.size':8})
    mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    roll_me.remove_blank_rows()

    #todo: fix the below two lines, I don't need this to get the matrix size...
    current_window = roll_me.get_window()
    filled_matrix = imputer.fit_transform(current_window)
    n_genes = filled_matrix.shape[1]
    coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))
    gene_names=list(current_window.columns.values)

    #zscore all data

    roll_me.zscore_all_data()

    #regulator-target labels
    edge_labels = [x for x in itertools.product(gene_names,repeat=2)]
    alpha_list = []
    print("Creating initial model...")
    for nth_window in range(0,total_window_number):
        #loop gets the window, gets the coefficients for that window, then increments the window
        current_window = roll_me.get_window()

        #check if any values are completely blank
        current_window = current_window
        filled_matrix = current_window.values
        #filled_matrix = imputer.fit_transform(current_window)
        current_lasso = LassoWrapper(filled_matrix)

        alpha_range = np.linspace(0, current_lasso.get_max_alpha())
        q_list = [current_lasso.cross_validate_alpha(a1) for a1 in alpha_range]

        alpha_table = zip(alpha_range, q_list)
        alpha_list.append(alpha_table)
        if pd['override-alpha'] == "false":
            (best_alpha,Qs) = max(alpha_table, key = lambda t: t[1])
        elif pd['override-alpha'] == "true":
            best_alpha = alpha
        coeff_mat = current_lasso.get_coeffs(best_alpha)
        coeff_matrix_3d[:,:,nth_window] = coeff_mat
        #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
        roll_me.next()

    #convert coeff model to edge list
    #takes a 3D matrix and gene names list, converts it to a 3D edge list. it assumes that the 3d matrix is symmetric and square.

    #initial model conversion
    initial_model = utility.create_3D_linked_list(edge_labels, coeff_matrix_3d, 'B')
    pd['edge_labels'] = edge_labels
    pd['total_window_number'] = total_window_number
    pd['alpha'] = best_alpha
    pd['alpha_list'] = alpha_list
    return(initial_model, roll_me, pd)

def permute_model(roll_me, pd):
    print("Running permutation test...")
    #start permutation test
    permuter = Permuter()
    #give it a roller object
    alpha = pd['alpha']
    edge_labels = pd['edge_labels']

    permuter.run_permutation_test(roll_me, alpha = alpha)
    perm_means=permuter.permutation_means
    perm_sd=permuter.permutation_sd

    permuted_model_means = utility.create_3D_linked_list(edge_labels, perm_means, 'p-means')
    permuted_model_sd = utility.create_3D_linked_list(edge_labels, perm_sd, 'p-sd')

    return(permuted_model_means, permuted_model_sd)

def bootstrap_model(roll_me, pd):
    print("Running bootstrapping test...")
    booter = Bootstrapper()
    boots = pd['boots']
    max_random = pd['max_random']
    n_alphas = pd['n_alphas']
    window_size = pd['window_size']
    total_window_number = pd['total_window_number']
    edge_labels = pd['edge_labels']

    booted_alphas = booter.run_bootstrap(roll_me, window_size, boots, n_alphas, noise=max_random)
    sums = np.sum(booter.freq_matrix,axis=3)
    auc = []
    for nth_window in range(0,total_window_number):
        auc.append(booter.get_nth_window_auc(0))
    auc = np.dstack(auc)
#get 3d coeff matrix for stability scores
    stability_model = utility.create_3D_linked_list(edge_labels, auc, 'stability')
    return(stability_model)

def aggregate_model(all_panels, pd):
    total_window_number = pd['total_window_number']
    results_table = [] # aggregated results for each window. index of list is the window number

    for nth_window in range(0,total_window_number):
        #merge panels
        aggregated_window = all_panels[0][nth_window].merge(all_panels[1][nth_window], on='regulator-target').merge(all_panels[2][nth_window], on='regulator-target').merge(all_panels[3][nth_window],on='regulator-target')
        results_table.append(aggregated_window)

#z-score again, find p-value using mean and standard deviation
    results_table_pvalues = []
    for nth_window in results_table:
        #don't get self edges in the p-value calculation. this is to avoid dividing by 0.
        valid_indicies = nth_window['p-sd'] != 0

        valid_window = nth_window[valid_indicies]
        initial_B = valid_window['B']
        sd = valid_window['p-sd']
        mean = valid_window['p-means']
        valid_window['final-z-scores-perm'] = (initial_B - mean)/sd
        #ensure that z score is negative to use cdf -> pvalue
        valid_window['cdf-perm'] = (-1*abs(valid_window['final-z-scores-perm'])).apply(scipy.stats.norm.cdf)
        #calculate t-tailed pvalue
        valid_window['p-value-perm'] = (2*valid_window['cdf-perm'])
        results_table_pvalues.append(valid_window)


    return(results_table_pvalues)

def rank_model(results_table_pvalues, pd):
    #Order and rank tables
    rank_by = pd['rank_by']
    ranked_results = utility.rank_results_3D(results_table_pvalues, rank_by)
    aggr_ranks = utility.average_rank(ranked_results, rank_by+"-rank")

    #Sort tables by mean rank in ascending order
    mean_sorted_edge_list = aggr_ranks.sort(columns="mean-rank", axis = 0)
    return(mean_sorted_edge_list)

def score_model(ranked_table, pd):
    gold_standard_file = pd['gold_standard_file']
    gold_standard_sep = pd['gold_standard_sep']
    evaluator = Evaluator(gold_standard_file, sep=gold_standard_sep)
    precision, recall, aupr = evaluator.calc_pr(ranked_table)
    return(precision, recall, aupr)

if __name__ == "__main__":
    main(sys.argv[1:])

