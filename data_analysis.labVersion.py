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


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


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
    epilog="Designed by Xavier Martinez on July 18th, 2016")
parser.add_argument('AngstromCutoff', type=float,
                    help='Cutoff value used to decide if water atom will be counted into the Hydration Number.')
parser.add_argument('pdb', type=valid_file, help='The .pdb file')
parser.add_argument('tpr', type=valid_file, help='The .tpr file.')
parser.add_argument('ndx', type=valid_file, help='The .ndx file.')
parser.add_argument('wd', type=valid_dir, help='The working directory (must use FAH directory structure).')
parser.add_argument('o', type=str, help='The name of the logfile.')

args = parser.parse_args()

ANGSTROM_CUTOFF = args.AngstromCutoff / 10.00
PDB = args.pdb
NDX = args.ndx
TPR = args.tpr
WORKING_DIRECTORY = args.wd
OUT_FILE = args.o

with open(OUT_FILE, mode='w') as logfile:
    logfile.write('{:5} {:3} {:4} {:>6} {:<7} {:<7} {:<4}\n'.format('P', 'R', 'C', 't(ps)', 'RMSD', 'RS', 'H'))

    topology = md.load(PDB).topology
    rna_atoms = [atom.index for atom in topology.atoms if
                 ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
    ow_atoms = topology.select('name O and water')

    rms_cmd = 'echo 1 1 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
    gyrate_cmd = 'echo 1 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
    hbond_cmd = 'echo 1 1 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'


    check_call = subprocess.check_call
    process = process_file
    unlink = os.unlink
    load = md.load
    compute_neighbors = md.compute_neighbors
    write = logfile.write

    for project_number in sorted(os.listdir(WORKING_DIRECTORY)):
        for run_number in sorted(os.listdir('%s/%s' % (WORKING_DIRECTORY, project_number))):
            for clone_number in sorted(os.listdir('%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number))):
                for xtc in os.listdir('%s/%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number, clone_number)):
                    assert '.xtc' in xtc, "%s is not an xtc." % xtc
                    start = clock()
                    xtc_file = '%s/%s/%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number, clone_number, xtc)
                    basename = xtc[:-4]
                    rms_file = basename + '.rms.xvg'
                    gyrate_file = basename + '.gyrate.xvg'
                    hbond_file = basename + '.hbond.xvg'
                    eprint("Processing %s.xtc" % basename)
                    try:
                        check_call(
                            rms_cmd.format(TPR, xtc_file, NDX, rms_file),
                            shell=True,
                            stderr=subprocess.STDOUT
                        )
                        check_call(
                            gyrate_cmd.format(TPR, xtc_file, NDX, gyrate_file),
                            shell=True,
                            stderr=subprocess.STDOUT
                        )
                        check_call(
                            hbond_cmd.format(TPR, xtc_file, NDX, hbond_file),
                            shell=True,
                            stderr=subprocess.STDOUT
                        )
                    except subprocess.CalledProcessError:
                        logging.fatal('ERROR: GROMACS not Loaded')
                        sys.exit()
                    rms = process(rms_file)
                    gyrate = process(gyrate_file)
                    hbond = process(hbond_file)
                    unlink(rms_file)
                    unlink(gyrate_file)
                    unlink(hbond_file)

                    traj = load(xtc_file, top=PDB)
                    eprint("Loaded trajectory file: %s" % str(traj))
                    eprint("Computing neighbors....")
                    neighbors = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=ow_atoms)
                    eprint("...finished.")
                    eprint("Writing data to %s" % OUT_FILE)
                    frames = []
                    for frame_index in range(traj.n_frames):
                        frame = int(traj.time[frame_index])
                        if frame in frames:
                            pass
                        else:
                            frames.append(frame)
                            count = len(neighbors[frame_index])
                            logdata = [
                                project_number,
                                run_number,
                                clone_number,
                                frame,
                                rms[frame],
                                gyrate[frame],
                                count
                            ]
                            write('{:5} {:3} {:4} {:>6} {:0<7.6} {:0<7.6} {:<4}\n'.format(*logdata))
                    eprint("Wrote data for %s unique frames." % (len(frames) - 1))
                    eprint("Finished processing %s.xtc. Runtime: %s seconds." % (basename, str(round((clock() - start), 2))))