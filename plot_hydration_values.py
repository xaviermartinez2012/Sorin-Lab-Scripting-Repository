import numpy as np
import matplotlib.pyplot as plt
import math
import random

# @profile
# def main():
values = {}
with open("luteo-analysis-trial2_reformat.txt", 'r') as data:
    for line in data:
        split = line.split()
        project, run, clone, time = split[0], split[1], split[2], split[3]
        hydrationValue = split[6]
        nucleotide_pairs = 0
        if split[7] != '0':
            x3dnaString = split[13]
            pairs = x3dnaString.split("_")
            pairs = pairs[:-1]
            nucleotide_pairs = len(pairs)
        if values.get(project):
            if values.get(project).get(run):
                if values.get(project).get(run).get(clone):
                    values.get(project).get(run).get(clone)[time] = (nucleotide_pairs, hydrationValue)
                else:
                    values[project][run][clone] = {time: (nucleotide_pairs, hydrationValue)}
            else:
                values[project][run] = {clone: {time: (nucleotide_pairs, hydrationValue)}}
        else:
            values[project] = {run: {clone: {time: (nucleotide_pairs, hydrationValue)}}}

random.seed(9)
maxH = 0
minH = 1000
maxPairs = 0

avg_rp = []
avg_H = []
proj = '1798'
for run in values.get(proj).keys():
    rgb = (random.random(), random.random(), random.random())
    clone_keys = values.get(proj).get(run).keys()
    clone_values = {}
    rNP = []
    rH = []
    for clone in clone_keys:
        cNP = []
        cH = []
        time_keys = values.get(proj).get(run).get(clone).keys()
        for time in time_keys:
            nucleotide_pairs, hydrationValue = values.get(proj).get(run).get(clone).get(time)
            cNP.append(nucleotide_pairs)
            cH.append(int(hydrationValue))
        clone_rp_avg = np.mean(cNP)
        clone_h_avg = np.mean(cH)
        plt.scatter(clone_rp_avg, clone_h_avg, s=10, c=rgb)
        rNP.append(clone_rp_avg)
        rH.append(clone_h_avg)
    if maxPairs < max(cNP):
        maxPairs = max(cNP)
    if maxH < max(cH):
        maxH = max(cH)
    if minH > min(cH):
        minH = min(cH)
    avg_rp.append(np.mean(rNP))
    avg_H.append(np.mean(rH))
plt.plot(np.unique(avg_rp), np.poly1d(np.polyfit(avg_rp, avg_H, 1))(np.unique(avg_rp)), c="red")
plt.xlabel("Nucleotide Pairs")
plt.ylabel("Hydration")
plt.xticks(np.arange(0, maxPairs, 1))
plt.yticks(np.arange(minH, maxH, 5))
plt.show()
# main()

