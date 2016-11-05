####################################################################################################################
# This script is to be used for checking if all folders and files exist within the Aquifex PKnot Projects (P1797/9)#
# Can also be used for projects with utilizing the same directory structure as the Aquifex PKnot Projects          #
# Ex: PROJ(N) -> RUN0 - RUN(M-1) -> CLONE0 - CLONE(P-1)                                                            #
# Where N is the project number, M is the number of RUN of directories, and P is the number of CLONE directories   #
# Written by Xavier Martinez                                                                                       #
# Version 1.0                                                                                                      #
####################################################################################################################

#!/usr/bin/env python

import os
import datetime as dt

global PROJECT_DIRECTORY
global RUN_DIRECTORIES
global CLONE_DIRECTORIES

# Obtain the working directory from the user and number of RUN & CLONE directories to check for.
correct_information = False
while not correct_information:
    PROJECT_DIRECTORY = raw_input("Enter the working directory for the project: ")
    try:
        RUN_DIRECTORIES = int(raw_input("Enter the number of RUN directories: "))
        CLONE_DIRECTORIES = int(raw_input("Enter the number of CLONE directories: "))
        if os.path.isdir(PROJECT_DIRECTORY):
            correct_information = True
            print ("\nWorking on %s. Checking for the existance of %s RUN directories and %s CLONE directories.") %(PROJECT_DIRECTORY, RUN_DIRECTORIES, CLONE_DIRECTORIES)
        else:
            print ("The working path you entered does not exist. Try again.")
    except ValueError:
        print ("You did not enter a number. Try again.")

# CD to working directory of the project and create report file
os.chdir(PROJECT_DIRECTORY)
time_stamp = dt.datetime.now()
report = open("report.%s.txt" % time_stamp.strftime("%m-%d-%Y"), "w")

# Write header/attribute information into the report (Project, RUN/CLONE directories checked, and  Date/Time of report)
report.write("Project: %s\nRUN Parameter: %s\nCLONE Parameter: %s\nDate/Time of report: %s\n\n" % (PROJECT_DIRECTORY, RUN_DIRECTORIES, CLONE_DIRECTORIES, time_stamp.strftime("%m-%d-%Y %H:%M")))

# Using nested loops, we iterate through the RUN and CLONE directories checking for their existance, as well as the existance of a .xtc file in CLONE.
## If either the directory or .xtc does not exist then write this into the report and increment the missing counter; This counter is then saved into the report as the number of files missing.
cwd = os.getcwd()
missing = 0
for x in xrange(RUN_DIRECTORIES):
    if not os.path.isdir("%s/RUN%s" % (cwd, x)):
        print "RUN%s does not exist.\n" % x
        report.write("RUN%s does not exist.\n" % x)
        missing += 1
    else:
        for y in xrange(CLONE_DIRECTORIES):
            if not os.path.isdir("%s/RUN%s/CLONE%s" % (cwd, x, y)):
                print "RUN%s/CLONE%s does not exist.\n" % (x, y)
                report.write("RUN%s/CLONE%s does not exist.\n" % (x, y))
                missing += 1
            else:
                found = False
                for files in os.listdir("%s/RUN%s/CLONE%s" % (cwd, x, y)):
                    if files.endswith(".xtc"):
                        found = True
                if not found:
                    print "No .xtc file in RUN%s/CLONE%s\n" % (x, y)
                    report.write("No .xtc file in RUN%s/CLONE%s\n" % (x, y))
                    missing += 1

report.write("\nTotal missing: %s" % missing)
report.close()
