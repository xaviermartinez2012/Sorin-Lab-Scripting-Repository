#!/usr/bin/env python

from __future__ import print_function
import os
import argparse


def valid_file(path):
    value = str(path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


parser = argparse.ArgumentParser(
    description="Determine whether differences between two analysis log files are from terminal work units.",
    epilog="Designed by Xavier Martinez on August 11th, 2016")
parser.add_argument('logfile', type=valid_file, help='The logfile to process.')
args = parser.parse_args()

LOG_FILE = args.logfile
triggers = {}
with open(LOG_FILE, 'r') as in_file:
    for line in in_file:
        split = line.split()[9:-3]
        split[0] = split[0][:-1]
        split[1] = split[1][:-1]
        split[2] = split[2][:-1]
        split[3] = split[3][:-1]
        in_key = (int(split[0]), int(split[1]), int(split[2]))
        if triggers.get(in_key):
            triggers.get(in_key).append(int(split[3]))
        else:
            triggers[(int(split[0]), int(split[1]), int(split[2]))] = [int(split[3])]
non_terminal_log = '{}, {}, {}, {} NO\n'
terminal_log = '{}, {}, {}, {} YES\n'
sorted_keys = sorted(triggers.keys())
non_terminal_count = 0
terminal_count = 0
with open('{}.nonterminal'.format(LOG_FILE), 'w') as non_terminal:
    with open('{}.terminal'.format(LOG_FILE), 'w') as terminal:
        for key in sorted_keys:
            timesteps = sorted(triggers.get(key))
            if len(timesteps) == 1:
                non_terminal_count += 1
                non_terminal.write(non_terminal_log.format(key[0], key[1], key[2], timesteps[0]))
            elif ((timesteps[0] - 100) % 1000) != 0:
                for time in timesteps:
                    non_terminal_count += 1
                    non_terminal.write(non_terminal_log.format(key[0], key[1], key[2], time))
            elif ((timesteps[-1] + 100) % 1000) != 0:
                for time in timesteps:
                    non_terminal_count += 1
                    non_terminal.write(non_terminal_log.format(key[0], key[1], key[2], time))
            else:
                is_terminal = True
                current = timesteps[0]
                for time in timesteps[1:]:
                    current += 100
                    if current == time:
                        pass
                    else:
                        is_terminal = False
                        break
                if is_terminal:
                    for time in timesteps:
                        terminal_count += 1
                        terminal.write(terminal_log.format(key[0], key[1], key[2], time))
                else:
                    for time in timesteps:
                        non_terminal_count += 1
                        non_terminal.write(non_terminal_log.format(key[0], key[1], key[2], time))
