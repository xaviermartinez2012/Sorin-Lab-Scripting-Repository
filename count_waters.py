#! /usr/bin/env python
from __future__ import print_function


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


import sys
import os
import argparse
import datetime as dt
import time
import numpy as np


def correct_comparison(atom):
    value = str(atom)
    if not (value == "OW" or value == "Na+"):
        raise argparse.ArgumentTypeError('AtomComparison must be HOH or Na+')
    return value


def valid_pdb(pdb_path):
    value = str(pdb_path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


parser = argparse.ArgumentParser(
    description="Obtain \"Water Number\" by calculating the distance between every water and every other non-water atom on the reference structure (.pdb). Values within the AngstromCutoff are counted.",
    epilog="Designed by Xavier Martinez on July 8th, 2016")
parser.add_argument('AngstromCutoff', type=int,
                    help='A cutoff value that is used to decide whether or not the water atom will be counted into the Water Number.')
parser.add_argument('AtomComparison', type=correct_comparison,
                    help='The atom type to compare to the RNA atoms (HOH or Na+)')
parser.add_argument('pdb', type=valid_pdb, help='The .pdb file.')
args = parser.parse_args()

ANGSTROM_CUTOFF = args.AngstromCutoff
ATOM_COMPARISON = args.AtomComparison
PDB = args.pdb

atom_comp_dict = {}
rna_dict = {}

with open(PDB) as pdb:
    for _ in xrange(4):
        next(pdb)
    for line in pdb:
        line_split = line.split()
        if len(line_split) < 3:
            pass
        else:
            if line_split[2] == ATOM_COMPARISON:
                atom_comp_dict[(int(float(line_split[1])), line_split[2])] = (
                    float(line_split[-5]), float(line_split[-4]), float(line_split[-3]))
            else:
                rna_dict[(int(float(line_split[1])), line_split[2])] = (
                float(line_split[-5]), float(line_split[-4]), float(line_split[-3]))

time_stamp = dt.datetime.now()
log_file = open("count_log_%s.txt" % time_stamp.strftime("%m-%d-%Y"), "w")
log_file.write(
    "Project: %s\nParameters: AngstromCutoff=%s AtomComparison=%s\n\n" % (PDB, str(ANGSTROM_CUTOFF), ATOM_COMPARISON))
log_file.write(
    "{}\t{}\t{}\t{}\t{}\n".format("Atom Number(1)", "Atom Type", "Measured Distance", "Atom Number(2)", "Atom Type"))

start_time = time.time()
count = 0
for atom_val in atom_comp_dict.iterkeys():
    for rna_val in rna_dict.iterkeys():
        atom = np.array(atom_comp_dict.get(atom_val))
        rna = np.array(rna_dict.get(rna_val))
        dist = np.linalg.norm(atom - rna)
        if dist <= ANGSTROM_CUTOFF:
            count += 1
            eprint("Count: %s" % str(count))
            log_file.write("{:14}\t{:9}\t{:17}\t{:14}\t{:9}\n".format(str(atom_val[0]), ATOM_COMPARISON, str(dist),
                                                                      str(rna_val[0]), rna_val[1]))

log_file.write("Total %s counted: %s" % (ATOM_COMPARISON, str(count)))
log_file.close()

eprint("Total Runtime: %s minutes." % str((time.time() - start_time) / 60.00))
