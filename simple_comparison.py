#!/usr/bin/env python
'''Simple comparison script.'''

import argparse
import os


def valid_file(path):
    '''Function to check the validity of a user-provided path.'''
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(
            '\"{}\" does not exist (must be in the same directory or specify full path).' % path)
    return path

PARSER = argparse.ArgumentParser(
    description='Simple Comparison Script.',
    epilog='Designed by Xavier Martinez on April 12th, 2016')
PARSER.add_argument('txt', type=valid_file, nargs='+',
                    help='The .txt file(s) to process. Must be 2 arguments.')
ARGS = PARSER.parse_args()

FILES = ARGS.txt
NUM_FILES = len(FILES)
assert NUM_FILES == 2, 'Required arguments 2. Arguments given: {}'.format(
    NUM_FILES)
TOTAL_DIFFERENCE = 0

FILE_ONE_LINES = []
FILE_TWO_LINES = []

with open(FILES[0], mode='r') as file_one:
    with open(FILES[1], mode='r') as file_two:
        FILE_ONE_LINES = file_one.readlines()
        FILE_TWO_LINES = file_two.readlines()

for file_one_line in FILE_ONE_LINES:
    hits = 0
    line_one = file_one_line.encode('UTF-8')
    line_one = file_one_line.strip()
    for file_two_line in FILE_TWO_LINES:
        line_two = file_two_line.encode('UTF-8')
        line_two = file_two_line.strip()
        if line_one == line_two:
            hits += 1
    if hits == 0:
        print '-- Line: {} in {} not in {}'.format(line_one, FILES[0], FILES[1])
        TOTAL_DIFFERENCE += 1
print '-- Total differences: {}'.format(TOTAL_DIFFERENCE)
