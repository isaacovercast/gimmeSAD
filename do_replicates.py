#!/usr/bin/env python2.7

from __future__ import print_function

import numpy as np
import subprocess
import time
import sys
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
    ## 1e6 sims normally just barely gets to equilibrium for most models
    #NSIMS=1000000
    ## Run all the way to double equilibrium
    NSIMS=0
    SIM_DIRECTORY="new_sims"
    #RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    ## Samples from log-uniform distribution
    col_rates = np.log10(np.array([0.001, 0.05]))
    col_rates = np.random.uniform(col_rates[0], col_rates[1], 3)
    col_rates = np.power(10, col_rates)
    ## Just do uniform for k values
    k_vals = np.random.random_integers(1000, 10000, 3)
    #k_vals = np.log10(np.array([1000, 10000]))
    #k_vals = np.random.uniform(k_vals[0], k_vals[1], 3)
    #k_vals = map(int, np.power(10, k_vals))

    for i in xrange(1):
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
                        + " -q "\
                        + " &"
                print(cmd)
                subprocess.call(cmd, shell=True)
        time.sleep(120)