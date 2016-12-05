#!/usr/bin/env python2.7

from __future__ import print_function

import subprocess
import time
import os

## This script will run each combination of col rate and K 100 times.
## So you'll end up with 400 (more) output directories in ./replicates.
##
## On my box the biggest model runs in a little less than 2 minutes, so
## this kind of hackishly runs 4x at a time and sleeps for 2 minutes, then
## starts the next ones.

## To get time to write stderr to an outfile you have to do something like this:
#$ (time ls) > outfile 2>&1
#$ (time ls) > ls_results 2> time_results

if __name__ == "__main__":
    NSIMS=1000000
    SIM_DIRECTORY="replicates"
    #RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    col_rates = [0.01, 0.05]
    k_vals = [1000, 5000]

    for i in xrange(100):
        print("Doing {}".format(i))
        t = str(time.time())
        for c in col_rates:
            for k in k_vals:
                sim_string = "K_"+str(k)+"-C_"+str(c) + "_" + t
                cursim_dir = os.path.join(SIM_DIRECTORY, sim_string)
                cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                        + " -c " + str(c)\
                        + " -k " + str(k)\
                        + " -o " + cursim_dir\
                        + " &"
                print(cmd)
                subprocess.call(cmd, shell=True)
        time.sleep(120)
