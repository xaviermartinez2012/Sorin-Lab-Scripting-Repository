#!/usr/bin/env python
# This is a test
from __future__ import print_function
import mdtraj as md
import sys
import os
import argparse
import logging
import subprocess
from time import clock

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')


# Function to print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# Function written by Dennis to parse the output of the rms, gyrate, and hbond Gromacs commands
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


# Function to check the existence of a file.
# Used in conjunction with argparse to check that the given parameter files exist.
def valid_file(path):
    value = str(path)
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


# Function to check the existence of a directory.
# Used in conjunction with argparse to check that the given parameter dataset exist.
def valid_dir(path):
    value = str(path)
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(
            '\"%s\" does not exist (must be in the same directory or specify full path).' % value)
    return value


# Function to check that the name of logfile does not conflict with the name of a temp file in use by the program.
def valid_outfile(name):
    value = str(name)
    if name == '3dnaout.txt':
        raise argparse.ArgumentTypeError(
            '\"{}\" is used within the program. Please choose another name for the outfile.'.format(name))
    return value


# Initialization of the argument parser.
parser = argparse.ArgumentParser(
    description='Generate a log containing rms, gyrate, hbond, hydation number, and X3DNA string for every frame off '
                'every .xtc file in a FAH dataset. Tested using Gromacs 3.3, X3DNA 2.2, and MDTraj 1.7.2. '
                'Relies on the Fix_3DNA_output.pl script (located on the wiki).',
    epilog="Designed by Xavier Martinez on July 20th, 2016")
parser.add_argument('AngstromCutoff', type=float,
                    help='Cutoff value used to decide if water atom will be counted into the Hydration Number.')
parser.add_argument('pdb', type=valid_file, help='The .pdb file')
parser.add_argument('tpr', type=valid_file, help='The .tpr file.')
parser.add_argument('nbp', type=valid_file, help='The native base pairs file.')
parser.add_argument('ndx', type=valid_file, help='The .ndx file. (Protein before system)')
parser.add_argument('wd', type=valid_dir, help='The working directory (must use FAH directory structure).')
parser.add_argument('o', type=str, help='The name of the logfile.')
args = parser.parse_args()

# Initialization of the variables that correspond to the arguments passed by the user.
ANGSTROM_CUTOFF = args.AngstromCutoff / 10.00
PDB = args.pdb
NDX = args.ndx
NBP = args.nbp
TPR = args.tpr
WORKING_DIRECTORY = args.wd
OUT_FILE = args.o

# Assert that the 'Fix_3DNA_output.pl' Pearl script exists within the working directory.
# This script is needed to parse output of the X3DNA program.
assert os.path.exists(
    'Fix_3DNA_output.pl'), 'Fix_3DNA_output.pl script not found in the directory and needed to continue...'

# Open both the logfile and datafile for writing.
# The log file will contain important information about the run quality of the program.
# The data file will contain purely data.
with open('{}.log'.format(OUT_FILE), mode='w') as logfile:
    # Load into memory logfile.write to improve performance
    log_write = logfile.write
    log_write(
        'Parameters:\nAngstromCutoff = {}\npdb = {}\ntpr = {}\nnbp = {}\nndx = {}\n\n'.format(
            ANGSTROM_CUTOFF, PDB, TPR, NBP, NDX))
    with open(OUT_FILE, mode='w') as datafile:
        # Load into memory function calls to improve performance
        check_call = subprocess.check_call
        process = process_file
        unlink = os.unlink
        load = md.load
        compute_neighbors = md.compute_neighbors
        data_write = datafile.write

        # Initialize the data fle with the columns of data that the program will find
        data_write(
            '{:5} {:3} {:4} {:>6} {:<7} {:<7} {:<4} {:<4}      {:<}\n'.format('P', 'R', 'C', 't(ps)', 'RMSD', 'RS', 'H',
                                                                        'Na+', 'X3DNA'))
        # Load just the topology once to select the atom indices of the rna(protein) and water.
        # This is used to calculate the wHydration number of a particular frame.
        # We do this once to save processing time since the topology does not change from frame to frame.
        topology = md.load(PDB).topology
        rna_atoms = [atom.index for atom in topology.atoms if
                     ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
        sodium_atoms = [atom.index for atom in topology.atoms if 'Na+' in atom.residue.name]
        ow_atoms = topology.select('name O and water')
        # Initialize the commands that will be used throughout the program, saving processing time.
        # Note: the trjconv command generates pdb's for every frame in a .xtc (important when using X3DNA)
        rms_cmd = 'echo 0 0 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
        gyrate_cmd = 'echo 0 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
        hbond_cmd = 'echo 0 0 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'
        trjconv_cmd = 'echo 0 | trjconv -f {} -s {} -sep -o frame.pdb -n {}'
        X3DNA_cmd = 'find_pair -z -p frame_{}.pdb pdbout_{}'
        fix_X3DNA_output_cmd = './Fix_3DNA_output.pl pdbout_{} {} > 3dnaout.txt'
        # These nested loops will walk through every xtc in every clone in every run in every project.
        for project_number in sorted(os.listdir(WORKING_DIRECTORY)):
            for run_number in sorted(os.listdir('%s/%s' % (WORKING_DIRECTORY, project_number))):
                for clone_number in sorted(os.listdir('%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number))):
                    for xtc in os.listdir(
                                    '%s/%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number, clone_number)):
                        # We assert that the file we are working with is in fact an xtc by checking its extension.
                        assert '.xtc' in xtc, "%s is not an xtc." % xtc
                        # A clock is started to measure the runtime of each iteration.
                        # This information is written to the logfile.
                        start = clock()
                        # We rebuild the directory structure of the xtc file
                        # to be able to run this program from any directory.
                        xtc_file = '%s/%s/%s/%s/%s' % (WORKING_DIRECTORY, project_number, run_number, clone_number, xtc)
                        # By 'slicing' out the last 4 characters of the xtc file, we obtain a 'basename'
                        # [Project]_[Run]_[Clone] that is used to label the rms, gyrate, and hbond temp files.
                        basename = xtc[:-4]
                        rms_file = basename + '.rms.xvg'
                        gyrate_file = basename + '.gyrate.xvg'
                        hbond_file = basename + '.hbond.xvg'
                        log_write("Processing %s.xtc\n" % basename)
                        # Call the Gromacs rms, gyrate, hbond, and trjconv commands.
                        # Exit the program if any one of the commands fails.
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
                            check_call(
                                trjconv_cmd.format(xtc_file, TPR, NDX),
                                shell=True,
                                stderr=subprocess.STDOUT
                            )
                        except subprocess.CalledProcessError as sub:
                            logging.fatal('ERROR: GROMACS not Loaded')
                            log_write(str(sub))
                            sys.exit()
                        # Parse the output of the rms, gyrate, and hbond commands.
                        # 'Unlink' (remove) these files when done.
                        # If not done, we could double or triple the size of the dataset.
                        rms = process(rms_file)
                        gyrate = process(gyrate_file)
                        hbond = process(hbond_file)
                        unlink(rms_file)
                        unlink(gyrate_file)
                        unlink(hbond_file)
                        # Using MDTraj we load the xtc and compute neighbors between the rna and waters using the cutoff
                        # This calculation gives us our hydration number
                        traj = load(xtc_file, top=PDB)
                        log_write("Loaded trajectory file: %s\n" % str(traj))
                        log_write("Computing neighbors....\n")
                        hydration_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=ow_atoms)
                        ion_density_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=sodium_atoms)
                        log_write("...finished.\n")

                        # Using the pdb's generated by trjconv we generate and parse a X3DNA string by (unique) frame.
                        # If the frame is not unique (Gromacs randomnly generates doubles of frames...) it is skipped
                        # Lastly we write into the datafile and unlink all temp files.
                        frames = []
                        append = frames.append
                        for frame_index in range(traj.n_frames):
                            frame = int(traj.time[frame_index])
                            # Check statement for repeated frames
                            if frame in frames:
                                unlink('frame_{}.pdb'.format(frame_index))
                                pass
                            else:
                                append(frame)
                                # X3DNA calcualtion and parsing (using Perl script)
                                try:
                                    check_call(
                                        X3DNA_cmd.format(frame_index, frame_index),
                                        shell=True,
                                        stderr=subprocess.STDOUT
                                    )
                                    check_call(
                                        fix_X3DNA_output_cmd.format(frame_index, NBP),
                                        shell=True,
                                        stderr=subprocess.STDOUT
                                    )
                                except subprocess.CalledProcessError as sub:
                                    logging.fatal('3DNA ERROR')
                                    log_write(str(sub))
                                    sys.exit()
                                # Obtain the hydration number by counting the number of atom indices (waters) are
                                # neighbors to the rna: d(rna, water) < angstrom_cutoff
                                hydration_number = len(hydration_values[frame_index])
                                ion_density = len(ion_density_values[frame_index])
                                # Write data to file
                                with open('3dnaout.txt', mode='r') as x3dna:
                                    for line in x3dna:
                                        data_write(
                                            '{:5} {:3} {:4} {:>6} {:0<7.6} {:0<7.6} {:<4} {:<4} {:<}'.format(project_number[1:],
                                                                                                       run_number[1:],
                                                                                                       clone_number[1:],
                                                                                                       frame,
                                                                                                       rms[frame],
                                                                                                       gyrate[frame],
                                                                                                       hydration_number,
                                                                                                       ion_density,
                                                                                                       line))
                                # Unlink all temp files
                                unlink('frame_{}.pdb'.format(frame_index))
                                # In the case that X3DNA does not find any pairs, the files in the try/catch are not
                                # generated. If not handled they cause the program to break.
                                try:
                                    unlink('allpairs.ana')
                                    unlink('allpairs.pdb')
                                    unlink('ref_frames.dat')
                                    unlink('mref_frames.dat')
                                    unlink('mulbp.inp')
                                    unlink('multiplets.pdb')
                                except OSError as e:
                                    log_write(str(e))
                                    pass
                                unlink('3dnaout.txt')
                                unlink('pdbout_{}'.format(frame_index))
                        # Indicate in the log file the the process finished and provide the runtime of the iteration.
                        log_write("Wrote data for %s unique frames.\n" % (len(frames) - 1))
                        log_write("Finished processing %s.xtc. Runtime: %s seconds.\n" % (
                            basename, str(round((clock() - start), 2))))
