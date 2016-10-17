#!/usr/bin/env python

from __future__ import print_function
import os
import argparse

# Function to check the existence of a file.
# Used in conjunction with argparse to check that the given parameter files exist.
def valid_file(path):
    value = str(path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value

# Initialization of the argument parser.
parser = argparse.ArgumentParser(
    description='Script to format the output of the data_analysis_v2.py script like Khai\'\s results.',
    epilog="Designed by Xavier Martinez on August 22nd, 2016")
parser.add_argument('f1', type=valid_file, help='The output of data_analysis_v2.py')
parser.add_argument('f2', type=valid_file, help='Khai\'\s results')
args = parser.parse_args()

FILE_ONE = args.f1
FILE_TWO = args.f2

f1_dict = {}

with open(FILE_ONE, 'r') as f1:
    f1.next()
    for line in f1:
        split1 = line.split()
        f1_3dna = ''
        if int(split1[7]) != 0:
            f1_3dna = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(split1[7], split1[8], split1[9], split1[10], split1[11], split1[12], split1[13])
        else:
            f1_3dna = '{}\t{}\t{}\t{}\t{}\t{}'.format(split1[7], split1[8], split1[9], split1[10], split1[11], split1[12])
        f1_dict[(split1[0][1:], split1[1][1:], split1[2][1:], split1[3])] = (split1[4], split1[5], split1[6], f1_3dna)

with open(FILE_TWO, 'r') as f2:
    with open('luteo-analysis-trial2_reformat.txt', 'w') as log_file:
        khais_format = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
        for line in f2:
            split2 = line.split()
            key = (split2[0], split2[1], split2[2], split2[3])
            if f1_dict.get(key):
                value = f1_dict.get(key)
                log_file.write(khais_format.format(key[0], key[1], key[2], key[3], value[0], value[1], value[2], value[3]))
            else:
                print('There is a missing timestep! {} {} {} {}'.format(key[0], key[1], key[2], key[3]))