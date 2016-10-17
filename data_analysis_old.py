#!/usr/bin/env python

from __future__ import print_function
import mdtraj as md
import numpy as np
import sys
import os
import logging
import subprocess

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


def unique_frames(traj):
    times = []
    for x in xrange(0, traj.n_frames):
        if round(traj[x].time[0]) not in times:
            times.append(round(traj[x].time[0]))
            yield x
    del times


def water_number(traj, angstrom_cutoff):
    frames = unique_frames(traj)
    topology = traj.topology
    ow_atoms = topology.select('name O and water')
    rna_atoms = [atom.index for atom in topology.atoms if
                 ('Na+' not in atom.residue.name and "HOH" not in atom.residue.name)]
    cutoff = angstrom_cutoff / 10.00
    for frame in frames:
        print("Processing Frame: %s" % (str(frame)))
        water_count = 0
        for comp_atom in ow_atoms:
            for rna in rna_atoms:
                a = traj.xyz[frame, comp_atom, :]
                b = traj.xyz[frame, rna, :]
                dist = np.linalg.norm(a - b)
                if dist <= cutoff:
                    water_count += 1
                    break
        print("Count: %s" % str(water_count))
        yield (int(round(traj.time[frame])), water_count)


with open('logfile-old.txt', mode='w') as logfile:
    xtc_file = "P1796_R0_C0.xtc"
    basename = xtc_file[:-4]
    ndx_file = "native_FAH_luteo.ndx"
    tpr_file = "frame0.tpr"
    rms_file = basename + '.rms.xvg'
    gyrate_file = basename + '.gyrate.xvg'
    hbond_file = basename + '.hbond.xvg'

    rms_cmd = 'echo 1 1 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
    gyrate_cmd = 'echo 1 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
    hbond_cmd = 'echo 1 1 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'

    try:
        subprocess.check_call(
            rms_cmd.format(tpr_file, xtc_file, ndx_file, rms_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
        subprocess.check_call(
            gyrate_cmd.format(tpr_file, xtc_file, ndx_file, gyrate_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
        subprocess.check_call(
            hbond_cmd.format(tpr_file, xtc_file, ndx_file, hbond_file),
            shell=True,
            stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError:
        logging.fatal('ERROR: GROMACS not Loaded')
        sys.exit()

    rms = process_file(rms_file)
    gyrate = process_file(gyrate_file)
    hbond = process_file(hbond_file)
    os.unlink(rms_file)
    os.unlink(gyrate_file)
    os.unlink(hbond_file)

    traj = md.load(xtc_file, top='native_FAH_luteo.pdb')

    logfile.write('{:4} {:>3} {:>3} {:>6} {:<7} {:<7} {:<5}\n'.format('P', 'R', 'C', 't(ps)', 'RMSD', 'RS', 'H'))
    for time, count in water_number(traj, 3.00):
        logdata = [
            1796,
            0,
            0,
            time,
            rms[time],
            gyrate[time],
            count
        ]
        logfile.write('{:4} {:>3} {:>3} {:>6} {:0<7.6} {:0<7.6} {:<5}\n'.format(*logdata))
