#!/usr/bin/env python3

import logging
import os.path
import sys
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


with open("pknot-analysis-test.txt", 'w') as logfile:

    rms_cmd = 'echo 1 1 | g_rms -noxvgr -s {} -f {} -n {} -o {}'
    gyrate_cmd = 'echo 1 | g_gyrate -noxvgr -s {} -f {} -n {} -o {}'
    hbond_cmd = 'echo 1 1 | g_hbond -noxvgr -s {} -f {} -n {} -num {}'

    xtc_file = "P1796_R0_C0.xtc"
    basename = xtc_file[:-4]
    ndx_file = "native_FAH_luteo.ndx"
    tpr_file = "frame0.tpr"
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
            1796,
            0,
            0,
            time,
            rms[time],
            gyrate[time]
        ]
        logfile.write(logtemplate.format(*logdata))
    os.unlink(rms_file)
    os.unlink(gyrate_file)
    os.unlink(hbond_file)
