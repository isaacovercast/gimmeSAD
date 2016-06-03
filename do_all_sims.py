#!/usr/bin/env python2.7

from __future__ import print_function

import subprocess
import os

if __name__ == "__main__":
    #NSIMS=10000
    NSIMS = 1000000000
    SIM_DIRECTORY="simout"
    RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    col_rates = [0.001, 0.01]
    k_vals = [1000, 10000, 100000]

    for c in col_rates:
        for k in k_vals:
            sim_string = "K_"+str(k)+"-C_"+str(c)
            cursim_dir = os.path.join(SIM_DIRECTORY, sim_string)
            cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                        + " -c " + str(c)\
                        + " -k " + str(k)\
                        + " -o " + cursim_dir\
                        + " -r " + str(RECORDING_INTERVAL)
            subprocess.call(cmd, shell=True)
