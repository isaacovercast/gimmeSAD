#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
import numpy as np
import collections      # For Counter
import msprime
import names
# pylint: disable=C0103
# pylint: disable=R0903


class species(object):

    def __init__(self, UUID = "", abundance = 1, meta_abundance = 1, colonization_time = 0):
        self.name = names.names().get_name()
        self.uuid = UUID
        self.abundance = abundance
        self.meta_abundance = meta_abundance
        self.colonization_time = colonization_time
        self.mutation_rate = .00000001
        self.Ne = self.colonization_time / 4 * abundance
        self.sequence_length = 500
        self.tree_sequence = []
        self.island_sample_size = 10
        self.meta_sample_size = 10

        ## Stats
        self.pi = 0
        self.S = 0
        self.S_island = 0
        self.S_meta = 0

    def __str__(self):
        return "<species {}>".format(self.name)

    def __repr__(self):
        return self.__str__()

    def simulate_seqs(self):
        island_pop = msprime.PopulationConfiguration(sample_size=self.island_sample_size, initial_size=1)
        meta_pop = msprime.PopulationConfiguration(sample_size=self.meta_sample_size, initial_size=float(self.meta_abundance)/self.abundance)
        ## TODO: Source and dest are reversed in current version of msprime, so if you update
        ## make sure this bug hasn't been fixed, if it has switch source and dest
        split_event = msprime.MassMigrationEvent(time=0.75, source=0, destination=1, proportion=1)
        self.tree_sequence = msprime.simulate(sample_size=20, length=self.sequence_length, Ne=self.Ne, mutation_rate=self.mutation_rate, \
                                population_configurations=[island_pop, meta_pop],\
                                demographic_events=[split_event])

        #self.tree_sequence = msprime.simulate(sample_size=10, Ne=self.Ne, length=self.sequence_length, mutation_rate=self.mutation_rate)

    def get_sumstats(self):
        self.pi = self.tree_sequence.get_pairwise_diversity()
        self.S = len(next(self.tree_sequence.haplotypes()))
#        total = 0.0
#        for i in range(1, n+1):
#            total += 1.0 / i
        all_haps = self.tree_sequence.haplotypes()
        ## Get population specific haplotypes
        island_haps = [next(all_haps) for _ in range(self.island_sample_size)]
        meta_haps = [next(all_haps) for _ in range(self.meta_sample_size)]

        ## Calculate S for each population
        ## Turn to an np array and transpose it, then test if all elements
        ## are the same by making a set from each row and checking if len > 1
        ihaps_t = np.transpose(np.array([list(x) for x in island_haps]))
        mhaps_t = np.transpose(np.array([list(x) for x in meta_haps]))

        ## Counter makes a dict, so just get the counts for 2, which indicates 
        ## sites segregating in the pop
        self.S_island = collections.Counter([len(set(ihaps_t[x])) for x in range(len(ihaps_t))])[2]
        self.S_meta = collections.Counter([len(set(mhaps_t[x])) for x in range(len(mhaps_t))])[2]
        

if __name__ == "__main__":
    from tabulate import tabulate
    import implicit_space
    from ascii_graph import Pyasciigraph

    data = implicit_space.implicit_space()
    #data.set_metacommunity("uniform")
    data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(50000):
        if not i % 1000:
            print("Done {}".format(i))
            #print(i, len(data.local_community), len(set(data.local_community)))
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    #print("\n\nSpecies abundance distribution (fq class, count):\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()

    ## Cool ascii bar graph
    ## The bar grapher buddy doesn't like ints for description so you have to transform it
    abundance_distribution = [(str(k), v) for k, v in abundance_distribution.items()]
    graph = Pyasciigraph(graphsymbol="|")
    print("\n")
    for line in  graph.graph('Simulated Species Abundance Distribution', abundance_distribution):
        print(line)
    print("###############################################################################\n")

    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()
    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Colonization time", "Local Abundance", "Metapopulation Abundance", "pi", "S", "S_island", "S_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, s.meta_abundance, s.pi, s.S, s.S_island, s.S_meta] for s in sp]
    print(tabulate(acc, headers, floatfmt=".4f"))
