#!/usr/bin/env python

from __future__ import print_function
import mdtraj as md
import sys
import os
import argparse
import logging
import subprocess

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')


def valid_file(path):
    value = str(path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


def valid_outfile(name):
    value = str(name)
    if name == '3dnaout.txt':
        raise argparse.ArgumentTypeError(
            '\"{}\" is used within the program. Please choose another name for the outfile.'.format(name))
    return value


def main():
    parser = argparse.ArgumentParser(
        description='Debug script for data_analysis_3DNA.py. Runs on P1796_R0_C0.xtc and provides verbose output on '
                    'every calculation. Designed and tested using Gromacs 3.3, X3DNA 2.2, and MDTraj 1.7.2.',
        epilog='Designed by Xavier Martinez on July 20th, 2016')
    parser.add_argument('pdb', type=valid_file, help='The .pdb file')
    parser.add_argument('tpr', type=valid_file, help='The .tpr file.')
    parser.add_argument('nbp', type=valid_file, help='The native base pairs file.')
    parser.add_argument('ndx', type=valid_file, help='The .ndx file. (Protein before system)')
    parser.add_argument('o', type=valid_outfile, help='The name of the logfile.')

    args = parser.parse_args()

    PDB = args.pdb
    TPR = args.tpr
    NBP = args.nbp
    NDX = args.ndx
    OUT_FILE = args.o

    with open(OUT_FILE, mode='w') as logfile:
        xtc_file = 'P1796_R10_C102.xtc'
        trjconv_cmd = 'echo 0 | trjconv -f {} -s {} -sep -o frame.pdb -n {}'
        X3DNA_cmd = 'find_pair -z -p frame_{}.pdb pdbout_{}'
        fix_X3DNA_output_cmd = './Fix_3DNA_output.pl pdbout_{} {} > 3dnaout.txt'
        try:
            subprocess.check_call(
                trjconv_cmd.format(xtc_file, TPR, NDX),
                shell=True,
                stderr=subprocess.STDOUT
            )
        except subprocess.CalledProcessError:
            logging.fatal('ERROR: GROMACS not Loaded')
            sys.exit()

        traj = md.load(xtc_file, top=PDB)
        logfile.write('{:5} {:3} {:4} {:>6}      {:<}\n'.format('P', 'R', 'C', 't(ps)', '3DNA'))
        frames = []
        for frame_index in xrange(traj.n_frames):
            frame = int(traj.time[frame_index])
            if frame in frames:
                os.unlink('frame_{}.pdb'.format(frame_index))
                pass
            else:
                frames.append(frame)
                try:
                    subprocess.check_call(
                        X3DNA_cmd.format(frame_index, frame_index),
                        shell=True,
                        stderr=subprocess.STDOUT
                    )
                    subprocess.check_call(
                        fix_X3DNA_output_cmd.format(frame_index, NBP),
                        shell=True,
                        stderr=subprocess.STDOUT
                    )
                except subprocess.CalledProcessError:
                    logging.fatal('Fix_3DNA ERROR')
                    sys.exit()
                with open('3dnaout.txt', mode='r') as x3dna:
                    for line in x3dna:
                        logfile.write('{:5} {:3} {:4} {:>6} {:<}'.format('P1796', '0', '0', frame, line))
                
		os.unlink('frame_{}.pdb'.format(frame_index))
                try:
                    os.unlink('allpairs.ana')
                    os.unlink('allpairs.pdb')
                    os.unlink('ref_frames.dat')
                    os.unlink('mref_frames.dat')
                    os.unlink('mulbp.inp')
                    os.unlink('multiplets.pdb')
                except OSError as e:
                    logfile.write(str(e))
                    pass
                os.unlink('3dnaout.txt')
                os.unlink('pdbout_{}'.format(frame_index))


main()
