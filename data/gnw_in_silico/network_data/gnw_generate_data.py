__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import os
import subprocess


if __name__ == '__main__':
    """
    READ FIRST!
    IMPORTANT NOTE - the jar file doesn't seem to use the --output-path flag appropriately. Files will be saved wherever
    this script is run
    """

    directory = '/Users/jfinkle/Documents/Northwestern/MoDyLS/Python/Roller/data/gnw_in_silico/'
    jar_location = '/Users/jfinkle/Documents/Northwestern/MoDyLS/gnw/'
    jar_file = jar_location + 'gnw-3.1.2b.jar'
    sim_settings = jar_location + 'sandbox/settings.txt'
    structures = directory + 'network_structures/'
    output = directory + 'network_data/'

    # Get file list
    for (path, names, filenames) in os.walk(structures):
        for filename in filenames:
            if '.xml' in filename:
                network = filename.replace('.xml', '')
                # if not os.path.exists(output+network):
                #     os.makedirs(output+network)
                # Simulate
                subprocess.call(['java', '-jar', jar_file, '--simulate', '-c', sim_settings, '--input-net', structures+filename,
                                 '--output-path', output+network])
