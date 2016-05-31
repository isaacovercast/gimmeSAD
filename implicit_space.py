#!/usr/bin/env python2.7
""" Sample object """

import matplotlib.pyplot as plt
from ascii_graph import Pyasciigraph
import collections
import numpy as np
import itertools
import random
import os

from species import species

# pylint: disable=C0103
# pylint: disable=R0903

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 150

class implicit_space(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self, K=10000, colrate=0.001, quiet=False):
        self.quiet = quiet

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []

        ## Settings specific to the uniform metacommunity
        self.uniform_inds = 1000
        self.uniform_species = 1000

        ## Variables associated with the metacommunity (these are poorly named)
        ## total_inds has to be set by set_metacommunity
        ## This is total individuals in the metacommunity
        self.total_inds = 0
        self.immigration_probabilities = []
        self.abundances = []
        self.species = []

        self.maxabundance = 0
        self.colonization_rate = colrate

        ## Variables for tracking the local community
        self.local_community = []
        self.local_inds = K
        self.divergence_times = {}

        self.current_time = 0


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

        ## Actually using uuid is v slow
        #self.species = [uuid.uuid4() for _ in enumerate(self.abundances)]
        self.species = [x for x in enumerate(self.abundances)]
        self.total_inds = sum(self.abundances)
        self.immigration_probabilities = [float(self.abundances[i])/self.total_inds for i in range(len(self.abundances))]

        self.maxabundance = np.amax(self.immigration_probabilities)
        

    def __str__(self):
        return "<implicit_space {}>".format(self.name)


    def prepopulate(self, mode="landbridge"):
        if mode == "landbridge":
            ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
            init_community = np.random.multinomial(self.local_inds, self.immigration_probabilities, size=1)        
            print("Initializing local community:")
            for i, x in enumerate(init_community[0]):
                if x:
                    self.local_community.append([self.species[i]] * x)
            print("N species = {}".format(len(self.local_community)))
            self.local_community = list(itertools.chain.from_iterable(self.local_community))
            print("N individuals = {}".format(len(self.local_community)))

            ## Divergence times for the landbridge model are dangerously broken.
            for taxon in set(self.local_community):
                self.divergence_times[taxon] = 0
        else:
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            self.local_community = [0] * self.local_inds
            new_species = self.species[self.immigration_probabilities.index(self.maxabundance)]
            self.local_community[0] = new_species
            self.divergence_times[new_species] = 1


    def step(self, nsteps=1):
        for step in range(nsteps):
            ## If there are any members of the local community
            if self.local_community:
                ## Select the individual to die
                #victim = random.choice(self.local_community)
                #self.local_community.remove(victim)
                #print("Removing victing - {}".format(victim))
                self.local_community.remove(random.choice(self.local_community))
                
            ## Check probability of an immigration event
            if np.random.random_sample() < self.colonization_rate:
    
                ## Loop until you draw species unique in the local community
                ## The flag to tell 'while when we're done, set when you successfully 
                ## draw a non-local-doop from the metacommunity
                unique = 0
    
                ## If you set your carrying capacity too high relative to the size of your
                ## metacommunity then you'll get stuck drawing duplicates over and over
                idiot_count = 0
                while not unique:
                    ## Sample from the metacommunity
                    migrant_draw = np.random.multinomial(1, self.immigration_probabilities, size=1)
                    #print("Immigration event - {}".format(np.where(migrant_draw == 1)))
                    #print("Immigrant - {}".format(self.species[np.where(migrant_draw == 1)[1][0]]))
                    new_species = self.species[np.where(migrant_draw == 1)[1][0]]
    
                    ##TODO: Should set a flag to guard whether or not to allow multiple colonizations
                    if new_species in self.local_community:
                        #print("multiple colonization events are forbidden, for now")
                        unique = 0
    
                        if idiot_count > MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY:
                            msg = """Metacommunity is exhausted w/ respect to local
                            community. Either expand the size of the metacommunity,
                            decrease the carrying capacity, or switch on multiple
                            migration (unimplemented)."""
                            sys.exit()
                        idiot_count +=1
                    else:
                        self.local_community.append(new_species)
                        self.divergence_times[new_species] = self.current_time
                        unique = 1
            else:
                ## Sample from the local community, including empty demes
                ## Sample all available from local community (community grows slow in volcanic model)
                self.local_community.append(random.choice(self.local_community))
    
                ## Sample only from available extant species (early pops grow quickly in the volcanic model)
                ## If you do this, the original colonizer just overwhelms everything else
                ## This is more similar to the Rosindell and Harmon model, in which they simply
                ## prepopulate the island entirely with one species. This is effectively the same
                #self.local_community.append(random.choice([x for x in self.local_community if x]))
    
            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False):
        ## Make a counter for the local_community, counts the number of
        ## individuals w/in each species
        abundances = collections.Counter(self.local_community)

        ## If we were doing mode=volcanic then there may be some remaining
        ## space in our carrying capacity that is unoccupied (indicated by
        ## zeros in the ubundances.keys(), if there are no
        try:
            abundances.pop(0)
        except KeyError:
            pass

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


    def simulate_seqs(self):
        ## Setting colonization_time as a scaling factor rather than as a raw tdiv
        ## The old way of doing this is `self.current_time - tdiv`
        #self.species_objects = [species(UUID=UUID, colonization_time=1/float(tdiv), abundance=self.local_community.count(UUID),\
        self.species_objects = [species(UUID=UUID, colonization_time=self.current_time - tdiv, abundance=self.local_community.count(UUID),\
                            meta_abundance=self.abundances[self.species.index(UUID)]) for UUID, tdiv in self.divergence_times.items() if self.local_community.count(UUID)]
        for s in self.species_objects:
            s.simulate_seqs()
            s.get_sumstats()


    def get_species(self):
        return(self.species_objects)


if __name__ == "__main__":
    data = implicit_space()
    data.set_metacommunity("uniform")
    #data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(10000):
        if not i % 1000:
            print("Done {}".format(i))
        #print(i, len(data.local_community), len(set(data.local_community)))
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    print("Species abundance distribution:\n{}".format(abundance_distribution))
    print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()
    print("Species:\n{}".format(data.get_species()))
