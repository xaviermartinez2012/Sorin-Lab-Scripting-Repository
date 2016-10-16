#!/usr/bin/env python3

import logging
import os
import os.path
import sys
import subprocess

try:
    from projectdata import data
except ImportError:
    print('Cannot find projectdata.py Project Descriptiors')
    sys.exit()

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


usage = 'Usage: PROJECT_DIR= OUTFILE= ./pknot-analysis.py'
project_dir = os.environ.get('PROJECT_DIR')
outfile = os.environ.get('OUTFILE')

if len(sys.argv) > 1 or not (project_dir and outfile):
    print(usage)
    sys.exit()


logging.info('PROJECT_DIR = %s', project_dir)
logging.info('OUTFILE = %s', outfile)

with open(outfile, 'w') as logfile:
    location, dirs, files = next(os.walk(project_dir))

    rms_cmd = 'echo 1 1 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
    gyrate_cmd = 'echo 1 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
    hbond_cmd = 'echo 1 1 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'

    for xtc in files:
        if xtc[-4:] != '.xtc':
            print(xtc)
            continue
        basename = xtc[:-4]
        project, run, clone = basename.split('_')
        logging.info('Procesing %s %s %s', project, run, clone)

        xtc_file = os.path.join(location, xtc)
        ndx_file = data[project]['ndx']
        tpr_file = data[project]['tpr']
        rms_file = basename + '.rms.xvg'
        gyrate_file = basename + '.gyrate.xvg'
        hbond_file = basename + '.hbond.xvg'

        try:
            subprocess.check_call(
                rms_cmd.format(tpr_file, xtc_file, ndx_file, rms_file),
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.DEVNULL
            )
            subprocess.check_call(
                gyrate_cmd.format(tpr_file, xtc_file, ndx_file, gyrate_file),
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.DEVNULL
            )
            subprocess.check_call(
                hbond_cmd.format(tpr_file, xtc_file, ndx_file, hbond_file),
                shell=True,
                stderr=subprocess.STDOUT,
                stdout=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            logging.fatal('ERROR: GROMACS not Loaded')
            sys.exit()

        rms = process_file(rms_file)
        gyrate = process_file(gyrate_file)
        hbond = process_file(hbond_file)

        logtemplate = '{:4} {:>3} {:>3} {:>6} {:0<7.6} {:0<7.6}\n'
        for time in sorted(rms):
            logdata = [
                int(project[1:]),
                int(run[1:]),
                int(clone[1:]),
                time,
                rms[time],
                gyrate[time]
            ]
            logfile.write(logtemplate.format(*logdata))
        os.unlink(rms_file)
        os.unlink(gyrate_file)
        os.unlink(hbond_file)
