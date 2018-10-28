#!/usr/bin/env python2.7
""" Sample object """

import matplotlib.pyplot as plt
from ascii_graph import Pyasciigraph
from scipy.stats import logser
import collections
import numpy as np
import itertools
import random
import glob
import sys
import os

from species import species

# pylint: disable=C0103
# pylint: disable=R0903

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500

class implicit_BI(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self, K=5000, colrate=0.01, allow_multiple_colonization=True, exponential=False, quiet=False):
        self.quiet = quiet

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []
        self.exponential = exponential

        ## Settings specific to the uniform metacommunity
        ## This is individuals per species
        self.uniform_inds = 10000000
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
        self.migrants = {}
        self.abundances_through_time = {}
        self.allow_multiple_colonization = True

        self.extinctions = 0
        self.colonizations = 0
        self.current_time = 0
        
        ## Vector for tracking lifetimes of excinct species
        ## We can plot this and average it to figure out how long
        ## species hang around in the local community
        self.extinction_times = []

        ## The invasive species identity
        self.invasive = -1
        ## Track how many invasives differentially survived
        self.survived_invasives = 0
        self.invasion_time = 0

        ## empirical data
        self.empirical_sample_sizes = []


    def set_empirical(self, empirical_dir):
        for f in glob.glob(empirical_dir + "/spider-fasta/*.fasta"):
            with open(f) as infile:
                ct = len(infile.readlines())
                self.empirical_sample_sizes.append(ct/2)
        #print("empircal sizes: {}".format(self.empirical_sample_sizes))        


    def set_metacommunity(self, infile):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there
        """
        if infile == "logser":
            ## Parameter of the logseries distribution
            p = .98
            self.abundances = logser.rvs(p, size=100)
        elif infile == "uniform":
            #for i in range(self.uniform_inds):
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
        return "<implicit_BI {}>".format(self.name)


    ## Every recording period log abundances per species in the local community
    def _log(self):
        for sp in set(self.local_community):
            idx = sp[0]
            if sp in self.abundances_through_time:
                self.abundances_through_time[idx].append(self.local_community.count(sp))
            else:
                self.abundances_through_time[idx] = [self.local_community.count(sp)]


    def prepopulate(self, mode="landbridge"):
        if mode == "landbridge":
            ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
            init_community = np.random.multinomial(self.local_inds, self.immigration_probabilities, size=1)        
            print("Initializing local community:")
            for i, x in enumerate(init_community[0]):
                if x:
                    self.local_community.extend([self.species[i][0]] * x)
                    ## Set founder flag
            self.local_community = [tuple([x, True]) for x in self.local_community]
            print("N species = {}".format(len(set([x[0] for x in self.local_community]))))
            #self.local_community = list(itertools.chain.from_iterable(self.local_community))
            print("N individuals = {}".format(len(self.local_community)))

            ## All species diverge simultaneously upon creation of the island.
            for taxon in self.local_community:
                self.divergence_times[taxon[0]] = 1
                self.migrants[taxon[0]] = 0
                self.abundances_through_time[taxon[0]] = [self.local_community.count(taxon)]
        else:
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            ## This is the old way that acts weird
            #self.local_community = [0] * self.local_inds
            #self.local_community = [((x,0), True) for x in self.local_community]
            new_species = (self.species[self.immigration_probabilities.index(self.maxabundance)][0], True)
            self.local_community.append(new_species)
            for i in range(1,self.local_inds):
                self.local_community.append((None,True))
            self.divergence_times[new_species[0]] = 1
            self.migrants[new_species[0]] = 0
            self.abundances_through_time[new_species[0]] = [1]


    def step(self, nsteps=1, time=0, invasion_time=100000, invasiveness=0.1):
        for step in range(nsteps):
            ## If there are any members of the local community
            if self.local_community:
                ## Select the individual to die
                victim = random.choice(self.local_community)
                ## If no invasive hasn't invaded then just do the normal sampling
                if self.invasive == -1:
                    self.local_community.remove(victim)
                else:
                    ## If invasiveness is less than the random value remove the invasive individual
                    ## else choose a new individual
                    if victim == self.invasive and np.random.rand() < invasiveness:
                        self.survived_invasives += 1
                        victim = random.choice(self.local_community)
                    self.local_community.remove(victim)
                ## Record local extinction events
                if not victim in self.local_community:
                    ## This was supposed to not record "extinctions" of empty deme space
                    ## but it fucks up the calculation of extinction rate
                    ##if not victim[0] == None:
                    if True:
                        self.extinctions += 1
                        try:
                            self.extinction_times.append(self.current_time - self.divergence_times[victim])
                        except:
                            ## The empty deme will make this freak
                            pass
                    ## If the invasive prematurely goes extinct just pick a new one
                    if victim == self.invasive:
                        print("invasive went extinct")
                        self.invasive = -1
                    ## blank all the recording for the extinct species
                    self.divergence_times[victim] = 0
                    self.migrants[victim[0]] = 0
                    self.abundances_through_time[victim[0]] = []

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
                    if new_species[0] in [x[0] for x in self.local_community]:
                        #print("multiple colonization events are forbidden, for now")
                        if self.allow_multiple_colonization:
                            self.migrants[new_species[0]] += 1
                        else:
                            ## This is the "every migrant is a new species block of code. FUCK!
                            #new_species = (self.current_time, new_species[1])
                            #self.species.append(new_species)
                            #self.divergence_times[(new_species[0], False)] = self.current_time
                            #self.migrants[new_species[0]] = 0
                            #self.abundances_through_time[new_species[0]] = [1]

                            self.colonizations += 1
                            self.migrants[new_species[0]] += 1
                        self.local_community.append((new_species[0], False))
                        unique = 1
    
                        if idiot_count > MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY:
                            msg = """\nMetacommunity is exhausted w/ respect to local
                            community. Either expand the size of the metacommunity,
                            decrease the carrying capacity, or switch on multiple
                            migration (unimplemented)."""
                            sys.exit(msg)
                        idiot_count +=1
                    else:
                        ## Only set the invasive species once at the time of next migration post invasion time
                        if self.invasive == -1 and time >= invasion_time and not invasion_time < 0:
                            self.invasive = (new_species[0], False)
                            print("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                            self.invasion_time = self.current_time
                        self.local_community.append((new_species[0], False))
                        self.divergence_times[(new_species[0], False)] = self.current_time
                        self.migrants[new_species[0]] = 0
                        self.abundances_through_time[new_species[0]] = [1]
                        self.colonizations += 1
                        unique = 1
            else:
                ## Sample from the local community, including empty demes
                ## Sample all available from local community (community grows slow in volcanic model)
                ## Also, lots of early turnover
                ## This all was only true before i implemented the full rosindell/harmon model,
                ## There are no more empty demes in the current config
                self.local_community.append(random.choice(self.local_community))
    
                ## Sample only from available extant species (early pops grow quickly in the volcanic model)
                ## If you do this, the original colonizer just overwhelms everything else
                ## This is more similar to the Rosindell and Harmon model, in which they simply
                ## prepopulate the island entirely with one species. This is effectively the same
                #self.local_community.append(random.choice([x for x in self.local_community if not x[0] == 0]))
    
            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False):
        ## Make a counter for the local_community, counts the number of
        ## individuals w/in each species
        abundances = collections.Counter([x[0] for x in self.local_community])

        ## If we were doing mode=volcanic then there may be some remaining
        ## space in our carrying capacity that is unoccupied (indicated by
        ## None in the abundances.keys()
        try:
            abundances.pop(None)
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

    ## How strong is the bottleneck? Strength should be interpreted as percent of local
    ## community to retain
    def bottleneck(self, strength=1):
        reduction = int(round(self.local_inds * strength))
        self.local_community = self.local_community[:reduction]

        ## First remove the extinct species from the species list
        pre = len(self.species_objects)
        self.species_objects = [s for s in self.species_objects if s.uuid in self.local_community]
        ## Update the extinction counter
        self.extinctions += (pre - len(self.species_objects))

        sp = self.species_objects
        ## Update abundances per species that survived the bottleneck
        for i, s in enumerate(sp):
            if s.uuid in self.local_community:
                abund = self.local_community.count(s.uuid)
                s.update_abundance(abund)
                self.species_objects[i] = s

    def simulate_seqs(self):
        self._log()
        self.species_objects = []
        ## Setting colonization_time as a scaling factor rather than as a raw tdiv
        ## The old way of doing this is `self.current_time - tdiv`
        #self.species_objects = [species(UUID=UUID, colonization_time=1/float(tdiv), abundance=self.local_community.count(UUID),\
        for UUID, sample_count in zip(self.divergence_times.keys(), self.empirical_sample_sizes):
            tdiv = self.divergence_times[UUID]
            #print(self.local_community)
            if UUID in self.local_community:
                meta_abundance = -1
                for x, y in self.species:
                    if UUID[0] == x:
                        meta_abundance = y
                tdiv = self.current_time - tdiv
                migration_rate = self.migrants[UUID[0]]/float(tdiv)

                ## Scale time to K-steps per generation
                tdiv = np.ceil(tdiv/float(self.local_inds))
                #meta_abundance = [x[1] for x in self.abundances if x[0] == UUID[0]]
                #meta_abundance = self.abundances[self.species.index(UUID[0])]
                abundance = self.local_community.count(UUID)
                #print(self.local_community)
                self.species_objects.append(species(UUID=UUID, colonization_time=tdiv,\
                                        exponential=self.exponential, abundance=abundance,\
                                        meta_abundance=meta_abundance, migration_rate=migration_rate,\
                                        abundances_through_time=self.abundances_through_time[UUID[0]],\
                                        sample_size=sample_count))

        for s in self.species_objects:
            try:
                s.simulate_seqs()
                s.get_sumstats()
                ## For debugging invasives
                #if s.abundance > 1000:
                #    print("\n{}".format(s))
            except:
                #print("  Sample size = {}".format(s.island_sample_size))
                pass

    def set_species(self, species_objects):
        self.species_objects = species_objects

    def get_species(self):
        return(self.species_objects)


if __name__ == "__main__":
    data = implicit_BI()
    #data.set_metacommunity("uniform")
    data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(100000):
        if not i % 10000:
            print("Done {}".format(i))
            #print(i, len(data.local_community), len(set(data.local_community)))
            #print(data.local_community)
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    print("Species abundance distribution:\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()
    print("Species:\n{}".format(data.get_species()))
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))
