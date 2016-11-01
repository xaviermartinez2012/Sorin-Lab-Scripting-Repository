#!/usr/bin/env python


import Queue
import threading
import mdtraj as md
import sys
import os
import argparse


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


class DataProcessor(threading.Thread):
    def __init__(self, threadID, name, queue, lock, logfile):
        threading.Thread.__init__(self)
        self.tID = threadID
        self.n = name
        self.q = queue
        self.l = lock
        self.log = logfile
        self.stop_event = threading.Event()

    def run(self):
        while not self.stop_event.isSet():
            process_data(self.q, self.l, self.log)

    def join(self, timeout=None):
        self.stop_event.isSet()
        super(DataProcessor, self).join(timeout)


def process_data(queue, lock, logfile):
    try:
        data = queue.get(True, 1)
    except Queue.Empty:
        return
    basename = (data.split("/")[-1])[:-4]
    basename_split = basename.split("_")
    project_number = basename_split[0][1:]
    run_number = basename_split[1][1:]
    clone_number = basename_split[2][1:]
    traj = load(data, top=PDB)
    hydration_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=ow_atoms)
    ion_density_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=sodium_atoms)
    with lock:
        frames = []
        append = frames.append
        for frame_index in range(traj.n_frames):
            frame = int(traj.time[frame_index])
            # Check statement for repeated frames
            if frame in frames:
                pass
            else:
                append(frame)
                # Obtain the solvent saturation and ion density by counting the number of atom indices are
                # neighbors to the rna: distance(rna, water) < angstrom_cutoff
                hydration_number = len(hydration_values[frame_index])
                ion_density = len(ion_density_values[frame_index])
                # Write data to file
                logfile.write(
                    '{:5} {:3} {:4} {:>6} {:<4} {:<4}\n'.format(project_number, run_number, clone_number, frame,
                                                                hydration_number, ion_density))
    queue.task_done()


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


# Initialization of the argument parser.
parser = argparse.ArgumentParser(
    description='Multi-threaded application to generate a log containing solvent saturation and ion density'
                'for every frame off every .xtc file in a FAH dataset. Implemented with MDTraj 1.7.2.',
    epilog="Designed by Xavier Martinez on October 31st, 2016")
parser.add_argument('AngstromCutoff', type=float,
                    help='Cutoff value used to decide if water atom will be counted into the Hydration Number.')
parser.add_argument('pdb', type=valid_file, help='The native state .pdb file')
parser.add_argument('wd', type=valid_dir, help='The working directory (must use FAH directory structure).')
parser.add_argument('o', type=str, help='The name of the logfile.')
parser.add_argument('t', type=int, default=2, help='The number of threads to run. Default is 2.')
args = parser.parse_args()

# Initialization of the variables that correspond to the arguments passed by the user.
ANGSTROM_CUTOFF = args.AngstromCutoff / 10.00
PDB = args.pdb
DATASET_DIRECTORY = args.wd
OUT_FILE = args.o
NUM_THREADS = args.t

# Initialize the commands that will be used throughout the program, saving processing time.
load = md.load
compute_neighbors = md.compute_neighbors

# Load just the topology once to select the atom indices of the rna(protein) and water.
# This is used to calculate the Hydration number and ion density of a particular frame.
# We do this once to save processing time since the topology does not change from frame to frame.
topology = md.load(PDB).topology
rna_atoms = [atom.index for atom in topology.atoms if
             ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
sodium_atoms = [atom.index for atom in topology.atoms if 'Na+' in atom.residue.name]
ow_atoms = topology.select('name O and water')

# Create a queue of xtcs for the threads to operate on.
xtc_queue = create_queue(DATASET_DIRECTORY)

# Initialization of threads and lock
threading_lock = threading.Lock()
threads = []
with open(OUT_FILE, mode='w') as datafile:
    for t in range(NUM_THREADS):
        thread = DataProcessor((t + 1), "T{}".format(t), xtc_queue, threading_lock, datafile)
        thread.setDaemon(True)
        thread.start()
        threads.append(thread)
    try:
        xtc_queue.join()
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.join()
