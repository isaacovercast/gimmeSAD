#!/usr/bin/env python2.7

import implicit_space
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from ascii_graph import Pyasciigraph

# pylint: disable=C0103
# pylint: disable=R0903


class gimmeSAD(object):

    def __init__(self):
        pass        

    def __str__(self):
        return "<gimmeSAD {}>".format(self.name)


def plot_abundances_ascii(abundance_distribution):
    """
    Plot the abundance distribution as returned by the get_abundances
    method. The abundance dist in this case is a dict of abundance
    classes and counts
    """
    ## Cool ascii bar graph
    ## The bar grapher buddy doesn't like ints for description so you have to transform it
    abundance_distribution = [(str(k), v) for k, v in abundance_distribution.items()]
    graph = Pyasciigraph(graphsymbol="|")
    print("\n")
    for line in  graph.graph('Simulated Species Abundance Distribution', abundance_distribution):
        print(line)
    print("###############################################################################\n")


def tabulate_sumstats(data):
    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()
    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Col time", "Local Abund", "Metapop Abund", "pi", "pi_net",  "S", "S_island", "pi_island", "S_meta", "pi_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, s.meta_abundance, s.pi, s.pi_net, s.S, s.S_island, s.pi_island, s.S_meta, s.pi_meta] for s in sp]
    print(tabulate(acc, headers, floatfmt=".4f"))

    
def heatmap_pi_dxy(data):
    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.pi, x.pi_island) for x in sp])
    max_pi = max([x[0] for x in pis])
    max_pi_island = max([x[1] for x in pis])
    print max_pi, max_pi_island

    ## Make the bins
    pi_bins = np.linspace(0, np.log(max_pi), 10)
    pi_island_bins = np.linspace(0, np.log(max_pi_island), 10)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for pi, pi_island in pis:
        count_pi = 0
        count_pi_island = 0
        while not np.log(pi) <= pi_bins[count_pi]:
            count_pi += 1
        while not np.log(pi_island) <= pi_island_bins[count_pi_island]:
            count_pi_island += 1
        ## increment the heatmap point this corresponds to
        heat[count_pi][count_pi_island] += 1
    plt.pcolor(heat,cmap=plt.cm.Reds)
    plt.axis([0, np.log(max_pi), 0, np.log(max_pi_island)])
    
    plt.show()
    plt.close()


if __name__ == "__main__":

    data = implicit_space.implicit_space()
    #data.set_metacommunity("uniform")
    data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(50000):
        if not i % 10000:
            print("Done {}".format(i))
            #print(i, len(data.local_community), len(set(data.local_community)))
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    #print("\n\nSpecies abundance distribution (fq class, count):\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()

#    ## Cool ascii bar graph
#    ## The bar grapher buddy doesn't like ints for description so you have to transform it
#    abundance_distribution = [(str(k), v) for k, v in abundance_distribution.items()]
#    graph = Pyasciigraph(graphsymbol="|")
#    print("\n")
#    for line in  graph.graph('Simulated Species Abundance Distribution', abundance_distribution):
#        print(line)
#    print("###############################################################################\n")
    implicit_space.plot_abundances_ascii(abundance_distribution)

    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()
    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Col time", "Local Abund", "Metapop Abund", "pi", "pi_net",  "S", "S_island", "pi_island", "S_meta", "pi_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, s.meta_abundance, s.pi, s.pi_net, s.S, s.S_island, s.pi_island, s.S_meta, s.pi_meta] for s in sp]
    ##TODO: Sort by colonization time?
    print(tabulate(acc, headers, floatfmt=".4f"))
