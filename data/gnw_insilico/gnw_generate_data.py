__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import os
import subprocess
import fileinput
import pandas as pd
import numpy as np
import time
import xml.etree.ElementTree as ET
import networkx as nx
from networkx.algorithms import isomorphism
import nxpd

class GnwWrapper(object):
    """
    Object that handles interactions with GNW command line inteface
    """
    def __init__(self, jar_path, settings):
        if os.path.isfile(jar_path):
            self.jar_call = ['java', '-jar', jar_path]
        else:
            raise IOError(('File %s does not exist' % jar_path))
        self.settings_default = settings
        self.devnull = open(os.devnull, 'w')

    def make_networks(self, input_net, out_fmt, out_path, n, arg_list, settings=None, out_name=None, graphviz=False):
        # Make the top level directory if it doesn't exist
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        # Change working directory because simulating data output-path doesn't work
        temp_path = out_path+'/temp/'
        try:
            os.chdir(temp_path)
        except:
            os.mkdir(temp_path)
            os.chdir(temp_path)
        if settings is None:
            settings = self.settings_default
        if out_name is None:
            out_name = input_net.split('/')[-1].split('.')[0]

        network_list = self.get_existing_networks(out_path)
        if len(network_list) == n:
            sys.exit('The requested number of networks already exists')
        if len(network_list) > 0:
            print('%i networks already exists. Generating %i more' % (len(network_list), n-len(network_list)))

        while len(network_list) < n:
            net_num = n-len(network_list)
            net_num = 5
            # Make network
            call_string = self.extract_sub_net(input_net, out_fmt, out_path, arg_list, settings, out_name)
            network_file = out_path + out_name + '-1.xml'

            # Anonymize gene names
            self.anonymize_genes(network_file)

            # Check if same structure exists in outfolder
            # Transform
            self.transform(network_file, 3, out_path)
            transformed_file = network_file.replace('.xml', '.tsv')
            net = self.gold_to_nx(transformed_file)
            unique_net = True
            for graph in network_list:
                graph_matcher = isomorphism.DiGraphMatcher(net, graph)
                if graph_matcher.is_isomorphic():
                    unique_net = False
                    break

            if not unique_net:
                print('isomorphic')
                os.remove(network_file)
                os.remove(transformed_file)
                continue

            # Simulate data
            # Delete files in the temp folder
            for fn in os.listdir(temp_path):
                os.remove(fn)
            print('Generating data for network %i' % net_num)
            subprocess.call(self.jar_call + ['--simulate', '-c', settings, '--input-net', network_file],
                            stdout=self.devnull, stderr=self.devnull)

            # The returncode doesn't seem to work. Count the number of files. Should be 44
            # Remove the base networks
            os.remove(network_file)
            os.remove(transformed_file)
            if len(os.listdir(temp_path)) == 44:
                # Rename files
                for fn in os.listdir(temp_path):
                    new_fn = fn.replace('1.', str(net_num)+'.').replace('1_', str(net_num)+'_')
                    os.rename(temp_path+fn, out_path+new_fn)

                # Save networks and add to list
                if graphviz:
                    graph_name = out_path + out_name + '-%i.png' % net_num
                    nxpd.draw(net, filename=graph_name, format='png', show=False, layout='circo')
                network_list.append(net)
            else:
                print('not simulated properly. check settings if this persists')

        # Cleanup - remove temporary directory and make settings file
        print('Cleaning up...')
        os.removedirs(temp_path)
        with open(out_path+'settings.txt', 'w') as f:
            f.write(call_string)

    def get_existing_networks(self, path):
        network_list = [self.gold_to_nx(path+file) for file in os.listdir(path) if 'goldstandard.tsv' in file]
        return network_list

    def gold_to_nx(self, file_name):
        """
        Make a gold standard file into a networkx DiGraph
        """
        df = pd.read_csv(file_name, sep='\t', header=None)
        nodes = np.append(df.loc[:, 0].values, df.loc[:, 1].values)
        df = df[df.loc[:, 2] == 1]
        edges = zip(df.loc[:, 0], df.loc[:, 1])
        net = nx.DiGraph()
        net.add_nodes_from(nodes)
        net.add_edges_from(edges)
        return net

    def anonymize_genes(self, file_path, in_place=True):
        """

        :return:
        """
        tree = ET.parse(file_path)
        model = tree.getroot()[0]
        species_list = model[2]
        species_dict = {s.attrib['name']: ('G'+str(idx+1)) for idx, s in enumerate(species_list)
                        if s.attrib['name'] != '_void_'}
        fout = fileinput.input(file_path, inplace=in_place)
        for line in fout:
            for species, anon_id in species_dict.items():
                line = line.replace(species, anon_id)
            # Print the modified line to the output buffer
            print(line)
        fout.close()
        return

    def extract_sub_net(self, input_net, out_fmt, out_path, arg_list, settings, out_name):
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        extract_call = self.jar_call + ['--extract', '-c', settings]
        call_list = ['--input-net', input_net, '--num-subnets=1', '--output-net-format=%i' % out_fmt, '--output-path',
                     out_path, '--network-name', out_name] + arg_list
        subprocess.call(extract_call+call_list, stdout=self.devnull, stderr=self.devnull)
        return ' '.join(extract_call+call_list)

    def transform(self, in_list, out_fmt, out_path, settings=None):
        """
        Transform one network file into another type
        """
        if type(in_list) is not list:
            in_list = [in_list]

        if settings is None:
            settings = self.settings_default

        for path in in_list:
            call_list = ['--transform', '-c', settings, '--input-net', path, ('--output-net-format='+str(out_fmt)),
                         '--output-path', out_path]
            subprocess.call(self.jar_call+call_list, stdout=self.devnull, stderr=self.devnull)


if __name__ == '__main__':
    """
    READ FIRST!
    IMPORTANT NOTE - the jar file doesn't seem to use the --output-path flag appropriately. Files will be saved wherever
    this script is run
    """

    directory = '/Users/jfinkle/Documents/Northwestern/MoDyLS/Python/Roller/data/gnw_insilico/'
    jar_location = '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/'
    path_dict = {'Yeast': '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/src/ch/epfl/lis/networks/'
                          'yeast_transcriptional_network_Balaji2006.tsv',
                 'Ecoli': '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/src/ch/epfl/lis/networks/'
                          'ecoli_transcriptional_network_regulonDB_6_7.tsv'}
    network = 'Ecoli'
    n_nodes = 1000
    num_nets = 5
    scc_fraction = 0.5
    scc_num = round(n_nodes*scc_fraction)
    mode = {'Yeast': ['--scc-seed', str(scc_num)], 'Ecoli': ['--random-seed']}
    # mode = {'Yeast': ['--random-seed'], 'Ecoli': ['--random-seed']}
    base_net = path_dict[network]
    jar_file = jar_location + 'gnw-3.1.2b.jar'
    sim_settings = directory + 'settings.txt'
    data_out = directory + 'network_data/' + network + str(n_nodes) + '/'
    optional_args = mode[network]+['--greedy-selection', '--subnet-size='+str(n_nodes)]

    gnw = GnwWrapper(jar_file, sim_settings)
    gnw.make_networks(base_net, 4, data_out, num_nets, optional_args, out_name=network+str(n_nodes))
