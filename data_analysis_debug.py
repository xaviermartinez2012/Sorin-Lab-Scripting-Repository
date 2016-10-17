#!/usr/bin/env python

from __future__ import print_function
import mdtraj as md
import sys
import os
import argparse
import logging
import subprocess
from time import clock

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')


def process_file(data_file):
    # <TIME, METRIC> hashtable in ps v. nm
    metric = {}

    with open(data_file, 'r') as data:
        for line in iter(data.readline, ''):
            # Skip Header Information
            if line[0] is '@' or line[0] is '#':
                continue
            # Grab significant data from line
            raw_time, raw_metric = line.split()[:2]
            time = int(float(raw_time))
            radius = 10 * float(raw_metric)
            metric[time] = radius

    return metric


def valid_file(path):
    value = str(path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


def valid_dir(path):
    value = str(path)
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


parser = argparse.ArgumentParser(
    description="Debug script for data_analysis.py. Runs on P1796_R0_C0.xtc and provides verbose output on every "
                "calculation. Designed and tested using Gromacs 3.3, X3DNA 2.2, and MDTraj 1.7.2.",
    epilog="Designed by Xavier Martinez on July 18th, 2016")
parser.add_argument('AngstromCutoff', type=float,
                    help='A cutoff value that is used to decide whether or not the water atom will be counted into the Water Number.')
parser.add_argument('pdb', type=valid_file, help='The .pdb file')
parser.add_argument('tpr', type=valid_file, help='The .tpr file.')
parser.add_argument('ndx', type=valid_file, help='The .ndx file.')
parser.add_argument('o', type=str, help='The name of the logfile.')

args = parser.parse_args()

ANGSTROM_CUTOFF = args.AngstromCutoff / 10.00
PDB = args.pdb
NDX = args.ndx
TPR = args.tpr
OUT_FILE = args.o

with open(OUT_FILE, mode='w') as logfile:
    start_time = clock()
    xtc_file = "P1796_R0_C0.xtc"
    basename = xtc_file[:-4]
    rms_file = basename + '.rms.xvg'
    gyrate_file = basename + '.gyrate.xvg'
    hbond_file = basename + '.hbond.xvg'

    rms_cmd = 'echo 1 1 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
    gyrate_cmd = 'echo 1 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
    hbond_cmd = 'echo 1 1 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'

    try:
        subprocess.check_call(
            rms_cmd.format(TPR, xtc_file, NDX, rms_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
        subprocess.check_call(
            gyrate_cmd.format(TPR, xtc_file, NDX, gyrate_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
        subprocess.check_call(
            hbond_cmd.format(TPR, xtc_file, NDX, hbond_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError:
        logging.fatal('ERROR: GROMACS not Loaded')
        sys.exit()

    rms = process_file(rms_file)
    gyrate = process_file(gyrate_file)
    hbond = process_file(hbond_file)
    print("rms")
    print(rms)
    print("gyrate")
    print(gyrate)
    print("hbond")
    print(hbond)
    os.unlink(rms_file)
    os.unlink(gyrate_file)
    os.unlink(hbond_file)

    traj = md.load(xtc_file, top=PDB)
    print("Loaded %s: %s" % (xtc_file, str(traj)))
    print("Frames within %s" % xtc_file)
    print([int(time) for time in traj.time])
    topology = traj.topology
    rna_atoms = [atom.index for atom in topology.atoms if
                 ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
    print("rna atom indices")
    print(rna_atoms)
    ow_atoms = topology.select('name O and water')
    print("ow atom indices")
    print(ow_atoms)
    neighbors = md.compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=ow_atoms)

    logfile.write('{:5} {:3} {:4} {:>6} {:<7} {:<7} {:<4}\n'.format('P', 'R', 'C', 't(ps)', 'RMSD', 'RS', 'H'))
    frames = []
    for x in xrange(traj.n_frames):
        frame = int(traj[x].time[0])
        if frame in frames:
            print("Frame: %s is not unique." % str(frame))
            pass
        else:
            frames.append(frame)
            print("Indices of ow atoms near rna")
            print(neighbors[x])
            count = len(neighbors[x])
            logdata = [
                1796,
                0,
                0,
                frame,
                rms[frame],
                gyrate[frame],
                count
            ]
            logfile.write('{:4} {:>3} {:>3} {:>6} {:0<7.6} {:0<7.6} {:<5}\n'.format(*logdata))
    print("Runtime is %s seconds." % str(round(clock() - start_time,2)))