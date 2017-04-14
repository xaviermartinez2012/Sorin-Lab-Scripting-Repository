#!/usr/bin/env python
'''Simple comparison script.'''

import argparse
import os
from difflib import SequenceMatcher


def valid_file(path):
    '''Function to check the validity of a user-provided path.'''
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(
            '\"{}\" does not exist (must be in the same directory or specify full path).' % path)
    return path


def comparison(file_one_name, file_two_name, file_one_lines, file_two_lines):
    '''Comparison Engine'''
    print '-- Comparing \"{}\" to \"{}\" --'.format(file_one_name, file_two_name)
    total_missing = 0
    for file_one_line in file_one_lines:
        hits = 0
        maximum_percent_match = 0
        corresponding_line = ''
        line_one = file_one_line.encode('UTF-8')
        line_one = file_one_line.strip()
        for file_two_line in file_two_lines:
            line_two = file_two_line.encode('UTF-8')
            line_two = file_two_line.strip()
            if line_one == line_two:
                hits += 1
                break
            else:
                relative_percent_match = SequenceMatcher(
                    None, line_one, line_two).ratio()
                if maximum_percent_match < relative_percent_match:
                    maximum_percent_match = relative_percent_match
                    corresponding_line = line_two
        if hits == 0:
            print '-- Line \"{}\" in {} not in {}'.format(line_one, file_one_name, file_two_name)
            print '\t\"{}\" matched \"{}\" with {} %% agreement\n'.format(line_one, corresponding_line, (maximum_percent_match * 100.0))
            total_missing += 1
    print '-- Total Missing from \"{}\" in \"{}\": {}\n'.format(file_one_name, file_two_name, total_missing)


if __name__ == "__main__":
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

    FILE_ONE_LINES = []
    FILE_TWO_LINES = []

    with open(FILES[0], mode='r') as file_one:
        with open(FILES[1], mode='r') as file_two:
            FILE_ONE_LINES = file_one.readlines()
            FILE_TWO_LINES = file_two.readlines()

    comparison(FILES[0], FILES[1], FILE_ONE_LINES, FILE_TWO_LINES)
    comparison(FILES[1], FILES[0], FILE_TWO_LINES, FILE_ONE_LINES)
