#!/usr/bin/env python


import Queue
import threading
import mdtraj as md
import os
import argparse
import pickle
from time import clock


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
        self.s_e = threading.Event()

    def run(self):
        while not self.s_e.isSet():
            try:
                data = self.q.get(True, 0.001)
            except Queue.Empty:
                continue
            basename = (data.split("/")[-1])[:-4]
            print("{} processing {}.xtc".format(self.n, basename))
            basename_split = basename.split("_")
            project_number = basename_split[0][1:]
            run_number = basename_split[1][1:]
            clone_number = basename_split[2][1:]
            traj = load(data, top=PDB)
            hydration_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=water_molecules)
            ion_density_values = compute_neighbors(traj, ANGSTROM_CUTOFF, rna_atoms, haystack_indices=sodium_atoms)
            with self.l:
                print('{} writing data to {}. Locking...'.format(self.n, OUT_FILE))
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
                        self.log.write(
                            '{:5} {:3} {:4} {:>6} {:<4} {:<4}\n'.format(project_number, run_number, clone_number, frame,
                                                                        hydration_number, ion_density))
                print('{} finished writing to {}. Unlocking...'.format(self.n, OUT_FILE))
            self.q.task_done()
            print('.xtc\'s in processing queue: {}'.format(self.q.qsize()))
        print('{} terminating...'.format(self.n))

    def join(self, timeout=None):
        self.s_e.set()
        super(DataProcessor, self).join(timeout)


def serialize_analysis(data_dict):
    output = open('solventSaturation.ionDensity.{}angstromCutoff.dataset.pkl'.format(ANGSTROM_CUTOFF), mode='wb')
    pickle.dump(data_dict, output, -1)


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
parser.add_argument('-t', metavar='threads', type=int, default=2, help='The number of threads to run. Default is 2.')
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
topology = load(PDB).topology
rna_atoms = [atom.index for atom in topology.atoms if
             ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
sodium_atoms = [atom.index for atom in topology.atoms if 'Na+' in atom.residue.name]
water_molecules = topology.select('name O and water')

# Create a queue of xtcs for the threads to operate on.
xtc_queue = create_queue(DATASET_DIRECTORY)
print('.xtc\'s in processing queue: {}'.format(xtc_queue.qsize()))

# Initialization of threads and lock
threading_lock = threading.Lock()
kill_switch = threading.Event()
pool = []
with open(OUT_FILE, mode='w') as datafile:
    for t in range(NUM_THREADS):
        thread = DataProcessor((t + 1), "T{}".format(t), xtc_queue, threading_lock, datafile)
        thread.setDaemon(True)
        thread.start()
        pool.append(thread)
    start = clock()
    xtc_queue.join()
    for thread in pool:
        thread.join()
    print('Total runtime: {} minutes.'.format(str(round((clock() - start)/60.00, 2))))

data_dict = {}
with open(OUT_FILE, 'r') as data:
    for line in data:
        split = line.split()
        project, run, clone, time = split[0], split[1], split[2], split[3]
        hydrationValue = split[4]
        ionDensity = split[5]
        if data_dict.get(project):
            if data_dict.get(project).get(run):
                if data_dict.get(project).get(run).get(clone):
                    data_dict.get(project).get(run).get(clone)[time] = [hydrationValue, ionDensity]
                else:
                    data_dict[project][run][clone] = {time: [hydrationValue, ionDensity]}
            else:
                data_dict[project][run] = {clone: {time: [hydrationValue, ionDensity]}}
        else:
            data_dict[project] = {run: {clone: {time: [hydrationValue, ionDensity]}}}

serialize_analysis(data_dict)
