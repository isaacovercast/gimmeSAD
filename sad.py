#!/usr/bin/env python2.7
""" Sample object """

import os
import numpy as np

# pylint: disable=C0103
# pylint: disable=R0903


class sad(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self):
        self.immigration_probabilities = []
        self.abundances = []
        self.total_inds = 100

        self.maxabundance = 0

    def set_metacommunity(self, infile):
        if infile == "random":
            for i in range(1, self.total_inds+1):
                self.immigration_probabilities.append(i*.001)
        else:
            if os.path.isfile(infile):
                with open(infile, 'r') as inf:
                    self.abundances = [int(line.split()[0]) for line in inf]
            else:
                raise Exception("Bad metacommunity input - ".format(infile))

            self.total_inds = sum(self.abundances)
            self.immigration_probabilities = [float(self.abundances[i])/self.total_inds for i in range(len(self.abundances))]

        self.maxabunance = max(self.immigration_probabilities)

    def __str__(self):
        return "<gimmeSAD sad object {}>".format(self.name)


if __name__ == "__main__":
    data = sad()
    data.set_metacommunity("metacommunity_LS4.txt")
    print(data.total_inds)
    print(data.immigration_probabilities)
    print(data.abundances)
