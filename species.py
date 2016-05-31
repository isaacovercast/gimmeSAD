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
        self.abundance = abundance * 1000
        self.meta_abundance = meta_abundance
        self.colonization_time = np.log(colonization_time)
        self.mutation_rate = .00001
        self.Ne = self.colonization_time / 4 * abundance
        self.sequence_length = 800
        self.tree_sequence = []
        self.island_sample_size = 10
        self.meta_sample_size = 10
        self.split_time = 0.1

        ## Stats
        self.pi = 0
        self.S = 0
        self.S_island = 0
        self.S_meta = 0
        self.pi_net = 0
        self.dxy = 0

    def __str__(self):
        return "<species {}>".format(self.name)

    def __repr__(self):
        return self.__str__()

    def simulate_seqs(self):
        ## This is hackish
        if self.Ne <= 0:
            self.Ne = 100
        island_pop = msprime.PopulationConfiguration(sample_size=self.island_sample_size, initial_size=1)
        meta_pop = msprime.PopulationConfiguration(sample_size=self.meta_sample_size, initial_size=float(self.meta_abundance)/self.abundance)
        ## TODO: Source and dest are reversed in current version of msprime, so if you update
        ## make sure this bug hasn't been fixed, if it has switch source and dest
        split_event = msprime.MassMigrationEvent(time=self.split_time, source=0, destination=1, proportion=1)
        self.tree_sequence = msprime.simulate(sample_size=20, length=self.sequence_length, Ne=self.Ne, mutation_rate=self.mutation_rate, \
                                population_configurations=[island_pop, meta_pop],\
                                demographic_events=[split_event])

        #self.tree_sequence = msprime.simulate(sample_size=10, Ne=self.Ne, length=self.sequence_length, mutation_rate=self.mutation_rate)

    def get_sumstats(self):
        self.pi = self.tree_sequence.get_pairwise_diversity() / self.sequence_length
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
        ## In the transposed array rows are bases and columns are individuals
        ihaps_t = np.transpose(np.array([map(int, list(x)) for x in island_haps]))
        mhaps_t = np.transpose(np.array([map(int, list(x)) for x in meta_haps]))

        ## Counter makes a dict, so just get the counts for 2, which indicates 
        ## sites segregating in the pop
        ## S will not always == S_island + S_meta. If a site is fixed in one pop and not
        ## present in the other then S will be less than the total. If a site is segragating
        ## in both pops then S will be greater than the total.
        ## There's probably a smarter way to do this....
        self.S_island = collections.Counter([len(set(ihaps_t[x])) for x in range(len(ihaps_t))])[2]
        self.S_meta = collections.Counter([len(set(mhaps_t[x])) for x in range(len(mhaps_t))])[2]

        ## Pass in the transposed arrays, since we already have them
        self.pi_island = get_pi(ihaps_t) / self.sequence_length
        self.pi_meta = get_pi(mhaps_t) / self.sequence_length

        ## get pairwise differences between populations while ignoring differences
        ## within populations (Dxy)
        self.dxy = get_dxy(ihaps_t, mhaps_t) / self.sequence_length

        ## pi_net
        self.pi_net = self.dxy - (self.pi_island + self.pi_meta)/2


def get_pi(haplotypes):
    ## If no seg sites in a pop then haplotypes will be 0 length
    if haplotypes.size == 0:
        return 0
    n = len(haplotypes[0])
    n_comparisons = float(n) * (n - 1) / 2

    pi = 0
    for hap in haplotypes:
        k = np.count_nonzero(hap)
        pi += float(k) * (n - k) / n_comparisons
    return(pi)


def get_dxy(ihaps_t, mhaps_t):
    ## If no seg sites in a pop then haplotypes will be 0 length
    if ihaps_t.size == 0 or mhaps_t.size == 0:
        return 0

    ## Number of comparisons is == to n_island * n_metacommunity`
    ## n_metacommunity
    n_island = ihaps_t.shape[1]
    n_meta = mhaps_t.shape[1]
    n_comparisons = n_island * n_meta

    dxy = 0
    len_seq = ihaps_t.shape[0]

    ## ibases and mbases are now a list of all bases at a particular
    ## site within each population
    ## There's probably a more elegant way to do this but I was
    ## gunning for readibility. Probably failed.
    for ibases, mbases in zip(ihaps_t, mhaps_t):
        nonzeros_island = np.count_nonzero(ibases)
        nonzeros_meta = np.count_nonzero(mbases)
        zeros_island = n_island - nonzeros_island
        zeros_meta = n_meta - nonzeros_meta

        dxy += (nonzeros_island * zeros_meta \
                + zeros_island * nonzeros_meta) / float(n_comparisons)
    return dxy

if __name__ == "__main__":
    from tabulate import tabulate
    import implicit_space

    data = implicit_space.implicit_space()
    data.set_metacommunity("metacommunity_LS4.txt")
    data.prepopulate(mode="volcanic")
    for i in range(50000):
        if not i % 10000:
            print("Done {}".format(i))
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    implicit_space.plot_abundances_ascii(abundance_distribution)

    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()
    tabulate_sumstats(data)
