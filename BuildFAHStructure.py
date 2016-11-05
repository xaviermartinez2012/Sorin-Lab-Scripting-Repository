#! /usr/bin/env python
import os


WORKING_DIRECTORY = '/home/xavier/LUTEO'
DATASET_DIRECTORY = '/home/xavier/LUTEO/luteo-xtcs-withwater'

os.chdir(DATASET_DIRECTORY)

FAH = {}
for i in os.listdir(os.getcwd()):
    if i.endswith(".xtc"):
        line = i[:-4]
        line_split = line.split("_")
        if not FAH.has_key(line_split[0]):
            FAH[line_split[0]] = {line_split[1]: {line_split[2]: i}}
        else:
            if not FAH[line_split[0]].has_key(line_split[1]):
                FAH[line_split[0]][line_split[1]] = {line_split[2]: i}
            else:
                FAH[line_split[0]][line_split[1]][line_split[2]] = i

for p_name in FAH.iterkeys():
    os.chdir(WORKING_DIRECTORY)
    os.mkdir('%s' % p_name)

for p_name in FAH.iterkeys():  
    os.chdir('%s/%s' % (WORKING_DIRECTORY, p_name)
    for run_name in FAH.get(p_name).iterkeys():
        os.mkdir('%s' % run_name)

for p_name in FAH.iterkeys():  
    for run_name in FAH.get(p_name).iterkeys():
        os.chdir('%s/%s/%s' % (WORKING_DIRECTORY, p_name, run_name))
        for clone_name in FAH.get(p_name).get(run_name).iterkeys():
            os.mkdir('%s' % clone_name)

for p_name in FAH.iterkeys():  
    for run_name in FAH.get(p_name).iterkeys():
        for clone_name in FAH.get(p_name).get(run_name).iterkeys():
            os.rename('%s/%s' % (DATASET_DIRECTORY, FAH.get(p_name).get(run_name).get(clone_name)), '%s/%s/%s/%s/%s' % (WORKING_DIRECTORY, p_name, run_name, clone_name, FAH.get(p_name).get(run_name).get(clone_name)))

print "SUCCESS"
