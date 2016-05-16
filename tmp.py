#!/usr/bin/env python2.7
""" Sample object """

import matplotlib.pyplot as plt
import collections
import numpy as np
import itertools
import random
import uuid
import os

# pylint: disable=C0103
# pylint: disable=R0903


class sad(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self):
        self.immigration_probabilities = []
        self.abundances = []
        self.species = []

        ## Settings specific to the uniform metacommunity
        self.uniform_inds = 1000
        self.uniform_species = 1000

        ## total_inds has to be set by set_metacommunity
        self.total_inds = 0

        self.maxabundance = 0
        self.im_rate = 0.001

        self.local_community = []
        self.local_inds = 10000


    def set_metacommunity(self, infile):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there
        """
        if infile == "uniform":
            for i in range(self.uniform_inds):
                self.abundances = [self.uniform_inds] * self.uniform_species
        else:
            if os.path.isfile(infile):
                with open(infile, 'r') as inf:
                    self.abundances = [int(line.split()[0]) for line in inf]
            else:
                raise Exception("Bad metacommunity input - ".format(infile))

        self.species = [uuid.uuid4() for _ in enumerate(self.abundances)]
        self.total_inds = sum(self.abundances)
        self.immigration_probabilities = [float(self.abundances[i])/self.total_inds for i in range(len(self.abundances))]

        self.maxabunance = max(self.immigration_probabilities)


    def __str__(self):
        return "<local_community {}>".format(self.name)


    def prepopulate(self):
        ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
        init_community = np.random.multinomial(self.local_inds, self.immigration_probabilities, size=1)        
        print("Initializing local community:")
        for i, x in enumerate(init_community[0]):
            if x:
                self.local_community.append([self.species[i]] * x)
        print("N species = {}".format(len(self.local_community)))
        self.local_community = list(itertools.chain.from_iterable(self.local_community))
        print("N individuals = {}".format(len(self.local_community)))


    def step(self):
        ## If there are any members of the local community
        if self.local_community:
            ## Select the individual to die
            self.local_community.remove(random.choice(self.local_community))
        ## Check probability of an immigration event
        if np.random.random_sample() < self.im_rate:
            ## Sample from the metacommunity
            migrant_draw = np.random.multinomial(1, self.immigration_probabilities, size=1)
            print("Immigration event - {}".format(np.where(migrant_draw == 1)))
            print("Immigrant - {}".format(self.species[np.where(migrant_draw == 1)[1][0]]))
            self.local_community.append(self.species[np.where(migrant_draw == 1)[1][0]])
        else:
            ## Sample from the local community
            self.local_community.append(random.choice(self.local_community))


    def get_abundances(self, octaves=False):
        ## Make a counter for the local_community, counts the number of
        ## individuals w/in each species
        abundances = collections.Counter(self.local_community)

        ## Make a set of abundances to get all the unique values
        abundance_classes = set(abundances.values())

        ## Now for each abundance class you have to go through and
        ## count the number of species at that abundance.
        ## This is currently stupid because in python there's no
        ## straightforward way to get keys from values in a dict.
        abundance_distribution = collections.OrderedDict()
        for i in abundance_classes:
            count = 0
            for _, v in abundances.items():
                if v == i:
                    count += 1
            abundance_distribution[i] = count
        if octaves:
            dist_in_octaves = collections.OrderedDict()
            min = 1
            max = 2
            while max/2 < len(abundance_distribution):
                ## Count up all species w/in each octave
                count = 0
                ## Here `i` is the abundance class and
                ## `j` is the count for that class
                for i, j in abundance_distribution.items():
                    if (i < max) and (i >= min):
                        count += j
                dist_in_octaves[min] = count
                min = min * 2
                max = max * 2
            abundance_distribution = dist_in_octaves
        return abundance_distribution

if __name__ == "__main__":
    data = sad()
    data.set_metacommunity("uniform")
    #data.set_metacommunity("metacommunity_LS4.txt")
    data.prepopulate()
    for i in range(100000):
        if not i % 1000:
            print("Done {}".format(i))
        #print(i, len(data.local_community), len(set(data.local_community)))
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    print(abundance_distribution)
    plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    plt.show()
