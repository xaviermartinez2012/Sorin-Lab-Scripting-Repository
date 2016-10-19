#!/usr/bin/env python


from __future__ import print_function
import Queue
import threading
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


def create_queue(working_directory):
    work_queue = Queue.Queue()
    listdir = os.listdir
    for project_number in listdir(working_directory):
        for run_number in listdir('%s/%s' % (working_directory, project_number)):
            for clone_number in listdir('%s/%s/%s' % (working_directory, project_number, run_number)):
                for xtc in os.listdir(
                                '%s/%s/%s/%s' % (working_directory, project_number, run_number, clone_number)):
                    # We assert that the file we are working with is in fact an xtc by checking its extension.
                    assert '.xtc' in xtc, "%s is not an xtc." % xtc
                    # We rebuild the directory structure of the xtc file
                    # to be able to run this program from any directory.
                    xtc_file = '%s/%s/%s/%s/%s' % (working_directory, project_number, run_number, clone_number, xtc)
                    work_queue.put(xtc_file)
    return work_queue


class DataProcessor (threading.Thread):
    def __init__(self, threadID, name, q, logfile):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
        self.log = logfile

    def run(self):
        process_data(self.q, self.log)


def process_data(queue, logFile):
    while not exitFlag:
        queue_lock.acquire()
        if not queue.empty():
            data = queue.get()
            basename = data[:-4]
            rms_file = basename + '.rms.xvg'
            gyrate_file = basename + '.gyrate.xvg'
            hbond_file = basename + '.hbond.xvg'
            # Call the Gromacs rms, gyrate, hbond, and trjconv commands.
            # Exit the program if any one of the commands fails.
            try:
                check_call(
                    rms_cmd.format(TPR, data, NDX, rms_file),
                    shell=True,
                    stderr=subprocess.STDOUT
                )
                check_call(
                    gyrate_cmd.format(TPR, data, NDX, gyrate_file),
                    shell=True,
                    stderr=subprocess.STDOUT
                )
                check_call(
                    hbond_cmd.format(TPR, data, NDX, hbond_file),
                    shell=True,
                    stderr=subprocess.STDOUT
                )
                check_call(
                    trjconv_cmd.format(data, TPR, NDX),
                    shell=True,
                    stderr=subprocess.STDOUT
                )
            except subprocess.CalledProcessError:
                logging.fatal('ERROR: GROMACS not Loaded')
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

        queue_lock.release()


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
DATASET_DIRECTORY = args.wd
OUT_FILE = args.o

# Assert that the 'Fix_3DNA_output.pl' Pearl script exists within the current working directory.
# This script is needed to parse output of the X3DNA program.
assert os.path.exists(
    'Fix_3DNA_output.pl'), 'Fix_3DNA_output.pl script not found in the directory and needed to continue...'

check_call = subprocess.check_call
process = process_file
unlink = os.unlink
load = md.load
compute_neighbors = md.compute_neighbors

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

# Create a queue of xtcs for the threads to operate on.
xtc_queue = create_queue(DATASET_DIRECTORY)
queue_lock = threading.Lock()
threads = []
thread_names = ["T1", "T2"]
thread_ID = 1
exitFlag = 0

with open(OUT_FILE, mode='w') as datafile:
    for t_name in thread_names:
        thread = DataProcessor(thread_ID, t_name, xtc_queue, datafile)
        thread.start()
        threads.append(thread)
        thread_ID += 1
    while not xtc_queue.empty():
        pass
    exitFlag += 1
    for t in threads:
        t.join()
