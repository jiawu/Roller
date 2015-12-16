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

    def make_networks(self, input_net, out_fmt, out_path, n, arg_list, settings=None, out_name=None):
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

        network_list = []
        while len(network_list) < n:
            net_num = n-len(network_list)
            print net_num
            # Make network
            self.extract_sub_net(input_net, out_fmt, out_path, arg_list, settings, out_name)
            network_file = out_path + out_name + '-1.xml'

            # Anonymize gene names
            self.anonymize_genes(network_file)

            # Check if same structure exists in outfolder
            # Transform
            self.transform(network_file, 3, out_path)
            transformed_file = network_file.replace('.xml', '.tsv')
            df = pd.read_csv(transformed_file, sep='\t', header=None)
            nodes = np.append(df.loc[:, 0].values, df.loc[:, 1].values)
            df = df[df.loc[:, 2] == 1]
            edges = zip(df.loc[:, 0], df.loc[:, 1])
            net = nx.DiGraph()
            net.add_nodes_from(nodes)
            net.add_edges_from(edges)
            unique_net = True
            for graph in network_list:
                graph_matcher = isomorphism.DiGraphMatcher(net, graph)
                if graph_matcher.is_isomorphic():
                    unique_net = False
                    print 'ISOMORPHIC!'
                    break

            if not unique_net:
                os.remove(network_file)
                os.remove(transformed_file)
                continue

            # Simulate data
            # Delete files in the temp folder
            for fn in os.listdir(temp_path):
                os.remove(fn)
            subprocess.call(self.jar_call + ['--simulate', '-c', settings, '--input-net', network_file],
                            stdout=self.devnull, stderr=self.devnull)

            # The returncode doesn't seem to work. Count the number of files. Should be 44
            # Remove the base networks
            os.remove(network_file)
            os.remove(transformed_file)
            if len(os.listdir(temp_path)) == 36:
                # Rename files
                for fn in os.listdir(temp_path):
                    new_fn = fn.replace('1.', str(net_num)+'.').replace('1_', str(net_num)+'_')
                    os.rename(temp_path+fn, out_path+new_fn)
                network_list.append(net)


            # If the simulation happened properly then save the network!


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
            for species, anon_id in species_dict.iteritems():
                line = line.replace(species, anon_id)
            # Print the modified line to the output buffer
            print line
        fout.close()
        return

    def extract_sub_net(self, input_net, out_fmt, out_path, arg_list, settings, out_name):
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        extract_call = self.jar_call + ['--extract', '-c', settings]
        call_list = ['--input-net', input_net, '--num-subnets=1', '--output-net-format=%i' % out_fmt, '--output-path',
                     out_path, '--network-name', out_name] + arg_list
        subprocess.call(extract_call+call_list, stdout=self.devnull, stderr=self.devnull)
        return

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


def anonymize_genes(structure_directory):
    for net in os.listdir(structure_directory):
        if '.xml' not in net:
            continue

        current_num = net.split('_')[3].split('-')[1].split('.')[0]
        tree = ET.parse(structure_directory+net)
        model = tree.getroot()[0]
        species_list = model[2]
        species_dict = {s.attrib['name']: ('G'+str(idx+1)) for idx, s in enumerate(species_list)
                        if s.attrib['name'] != '_void_'}

        new_file = structure_directory + 'anonymous/yeast_anon_' + current_num + '.xml'
        with open(new_file, "wt") as fout:
            with open(structure_directory+net, "rt") as fin:
                for line in fin:
                    for species, anon_id in species_dict.iteritems():
                        line = line.replace(species, anon_id)
                    fout.write(line)
        fout.close()

def transform(input, ouput, settings, out_format):
    pass


if __name__ == '__main__':
    """
    READ FIRST!
    IMPORTANT NOTE - the jar file doesn't seem to use the --output-path flag appropriately. Files will be saved wherever
    this script is run
    """

    directory = '/Users/jfinkle/Documents/Northwestern/MoDyLS/Python/Roller/data/gnw_insilico/'
    jar_location = '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/'
    network = 'Ecoli'
    path_dict = {'Yeast': '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/src/ch/epfl/lis/networks/'
                          'yeast_transcriptional_network_Balaji2006.tsv',
                 'Ecoli': '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/src/ch/epfl/lis/networks/'
                          'ecoli_transcriptional_network_regulonDB_6_7.tsv'}
    base_net = path_dict[network]
    jar_file = jar_location + 'gnw-3.1.2b.jar'
    sim_settings = directory + 'settings.txt'
    structures = directory + 'network_structures/' + network + '/'
    data_out = directory + 'network_data/' + network + '/'
    optional_args = ['--random-seed', '--greedy-selection', '--subnet-size=10']
    num_nets = 10
    make_files = False

    gnw = GnwWrapper(jar_file, sim_settings)
    gnw.make_networks(base_net, 4, structures, num_nets, optional_args, out_name=network)

    sys.exit()
    # Make files
    if make_files:
        for (path, names, filenames) in os.walk(structures):
            for filename in filenames:
                if '.xml' in filename:
                    network = filename.replace('.xml', '')
                    # if not os.path.exists(output+network):
                    #     os.makedirs(output+network)
                    # Simulate
                    code = subprocess.call(['java', '-jar', jar_file, '--simulate', '-c', sim_settings, '--input-net',
                                            structures+filename])

    sys.exit()

    # Remove files
    ii=1
    jj=1
    for (path, names, filenames) in os.walk(output):
            for filename in filenames:
                if ('.py' not in filename and 'timeseries.tsv' not in filename and 'goldstandard.tsv' not in filename) \
                        or 'noise' in filename or 'proteins' in filename:
                    os.remove(filename)
                if 'tsv' in filename:
                    net_num = filename.split('_')[3].split('-')[1]
                    if 'timeseries.tsv' in filename:
                        print net_num,
                        df = pd.read_csv(filename, sep='\t')
                        print df.columns
                    if 'goldstandard.tsv' in filename:
                        df = pd.read_csv(filename, sep='\t', header=None)
                        print df.head()
                        sys.exit()