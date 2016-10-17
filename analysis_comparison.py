#!/usr/bin/env python

from __future__ import print_function
import os
import argparse
from datetime import date

time_stamp = date.today()


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
    description='Comparison script test consistency between the output of the data_analysis_v2.py script with Khai\'\s results.',
    epilog="Designed by Xavier Martinez on August 3rd, 2016")
parser.add_argument('f1', type=valid_file, help='The output of data_analysis_v2.py')
parser.add_argument('f2', type=valid_file, help='Khai\'\s results')
parser.add_argument('c', type=float, help='Cutoff value')
args = parser.parse_args()

FILE_ONE = args.f1
FILE_TWO = args.f2
CUTOFF = args.c

missing_count = 0
differences = 0
non_terminal_work_units = 0
terminal_work_units = 0

f1_dict = {}

with open(FILE_ONE, 'r') as f1:
    for line in f1:
        split1 = line.split()
        project = split1[0]
        run = split1[1]
        clone = split1[2]
        time = split1[3]
        f1_rmsd = split1[4]
        f1_gyrate = split1[5]
        f1_3dna = '{}'
        if int(split1[7]) != 0:
            f1_3dna = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(split1[7], split1[8], split1[9], split1[10], split1[11],
                                                          split1[12], split1[13])
        else:
            f1_3dna = '{}\t{}\t{}\t{}\t{}\t{}'.format(split1[7], split1[8], split1[9], split1[10], split1[11],
                                                      split1[12])
        f1_dict[(project, run, clone, time)] = (
            '{:.3f}'.format(float(f1_rmsd)), '{:.3f}'.format(float(f1_gyrate)), f1_3dna)

with open(FILE_TWO, 'r') as f2:
    with open('analysis_comparison_c-{}_{}-{}-{}.log'.format(CUTOFF, time_stamp.month, time_stamp.day, time_stamp.year),
              'w') as log_file:
        log_file.write('# Project, Run, Clone, Time, RMSD Difference, Gyrate Difference, Terminal/Non-terminal\n')
        stack = []
        for line in f2:
            split2 = line.split()
            project = split2[0]
            run = split2[1]
            clone = split2[2]
            time = split2[3]
            f2_rmsd = float(split2[4])
            f2_gyrate = float(split2[5])
            key = (project, run, clone, time)
            # if len(stack) > 0 and stack[-1][2] != clone:
            #     terminal = True
            #     while len(stack) > 1:
            #         if not stack.pop(0)[-1]:
            #             terminal = False
            #     stack.pop()
            #     if terminal:
            #         log_file.write('###########TERMINAL############\n')
            #         terminal_work_units += 1
            if f1_dict.get(key):

                gyrate_difference = abs(float(f1_dict.get(key)[1]) - f2_gyrate)
                rmsd_difference = abs(float(f1_dict.get(key)[0]) - f2_rmsd)

                if rmsd_difference > CUTOFF or gyrate_difference > CUTOFF:
                    if int(time) == 0:
                        log_file.write('{} {} {} {} {} {}'.format(project, run, clone, time, rmsd_difference, gyrate_difference))
                        differences += 1

                    elif len(stack) != 10:
                        stack.append((project, run, clone, time, rmsd_difference, gyrate_difference, True))

                    elif len(stack) == 10:
                        if stack[-1][2] == clone:
                            non_terminal = False
                            while len(stack) > 0:
                                value = stack.pop(0)
                                if value[-1]:
                                    non_terminal = True
                                    log_file.write(
                                        '{} {} {} {} {} {} NT\n'.format(value[0], value[1], value[2], value[3],
                                                                        value[4], value[5]))
                                    differences += 1
                            if non_terminal:
                                non_terminal_work_units += 1
                                log_file.write('\n')

                        else:
                            terminal = False
                            while len(stack) > 0:
                                value = stack.pop(0)
                                if value[-1]:
                                    terminal = True
                                    log_file.write(
                                        '{} {} {} {} {} {} T\n'.format(value[0], value[1], value[2], value[3],
                                                                       value[4], value[5]))
                                    differences += 1
                            if terminal:
                                terminal_work_units += 1
                                log_file.write('\n')

                        stack.append((project, run, clone, time, rmsd_difference, gyrate_difference, True))
                else:
                    if int(time) == 0:
                        pass

                    elif len(stack) != 10:
                        stack.append((project, run, clone, time, rmsd_difference, gyrate_difference, False))

                    elif len(stack) == 10:
                        if stack[-1][2] == clone:
                            non_terminal = False
                            while len(stack) > 0:
                                value = stack.pop(0)
                                if value[-1]:
                                    non_terminal = True
                                    log_file.write(
                                        '{} {} {} {} {} {} NT\n'.format(value[0], value[1], value[2], value[3],
                                                                       value[4], value[5]))
                                    differences += 1
                            if non_terminal:
                                non_terminal_work_units += 1
                                log_file.write('\n')

                        else:
                            terminal = False
                            while len(stack) > 0:
                                value = stack.pop(0)
                                if value[-1]:
                                    terminal = True
                                    log_file.write(
                                        '{} {} {} {} {} {} T\n'.format(value[0], value[1], value[2], value[3],
                                                                       value[4], value[5]))
                                    differences += 1
                            if terminal:
                                terminal_work_units += 1
                                log_file.write('\n')

                        stack.append((project, run, clone, time, rmsd_difference, gyrate_difference, False))
            else:
                log_file.write(
                    '# Time missing! Encountered with P, R, C, T: {}, {}, {}, {}\n'.format(key[0],
                                                                                           key[1],
                                                                                           key[2],
                                                                                           key[3]))
                missing_count += 1
                differences += 1
        log_file.write(
            '# Missing files: {}\n# Terminal Work Units: {}\n# Non-Terminal Work Units: {}\n# Total differences: {}'.format(
                missing_count,
                terminal_work_units,
                non_terminal_work_units,
                differences))
