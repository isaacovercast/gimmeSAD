#!/usr/bin/env python2.7

from __future__ import print_function
from ascii_graph import Pyasciigraph
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
import argparse
import sys
import os

import implicit_space

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


def plot_abundances_gui(abundance_distribution):
    plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    plt.show()


def tabulate_sumstats(data):
    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()
    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Col time", "Local Abund", "Metapop Abund", "pi", "pi_net",  "S", "S_island", "pi_island", "S_meta", "pi_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, s.meta_abundance, s.pi, s.pi_net, s.S, s.S_island, s.pi_island, s.S_meta, s.pi_meta] for s in sp]
    return tabulate(acc, headers, floatfmt=".4f")

    
def heatmap_pi_dxy(data):
    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.pi, x.pi_island) for x in sp])
    max_pi = max([x[0] for x in pis])
    max_pi_island = max([x[1] for x in pis])
    print(max_pi, max_pi_island)

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




def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 

    gimmeSAD -n 100000                  ## Run the simulation for 100000 steps

    gimmeSAD -m metacommunity_LS4.txt   ## Source metacommunity from a file of ints
                                        ## the other option is to pass "uniform"
                                        ## which will generate a random metacommunity
                                        ## with all species at uniform frequency

    gimmeSAD -p volcanic                ## Set the mode for prepopulating (or not)
                                        ## the local community.

    gimmeSAD -c 0.001                   ## Set colonization rate

    gimmeSAD -o                         ## Output abundances in octaves
    gimmeSAD -q                         ## Don't be so chatty

   """)

    ## add model arguments 
    parser.add_argument('-n', metavar='nsims', dest="nsims", type=int,
        default=50000,
        help="Number of demographic events to simulate")

    parser.add_argument('-k', metavar='K', dest="K", type=int,
        default=10000,
        help="Carrying capacity of the island (max # individuals in the local community")

    parser.add_argument('-m', metavar='meta', dest="meta", type=str,
        default="metacommunity_LS4.txt",
        help="Source metacommunity from file or generate uniform")

    parser.add_argument('-c', metavar='colrate', dest="colrate", type=float,
        default=0.001,
        help="Set colonization rate")

    parser.add_argument('-p', metavar='mode', dest="mode", type=str,
        default="volcanic",
        help="Select mode for prepopulating the island (volcanic/landbridge)")

    parser.add_argument('-o', dest="octaves", action='store_true',
        help="Count species abundances in octaves")

    ## Add standard quiet/force/version args
    parser.add_argument('-q', "--quiet", action='store_true',
        help="do not print to stderror or stdout.")
    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")
    parser.add_argument('-v', '--version', action='version',
        version="0.0.0.1")

    ## parse args
    args = parser.parse_args()

    ## Check args
    if not args.meta == "uniform":
        if not os.path.exists(args.meta):
            parser.print_help()
            sys.exit("Metacommunity file doesn't exist - {}".format(args.meta))

    if args.mode not in ["volcanic", "landbridge"]:
        parser.print_help()
        sys.exit("-p option must be either volcanic or landbridge. You said - {}".\
            format(args.mode))

    return args


def progressbar(nsims, finished, msg=""):
    """ prints a progress bar """
    progress = 100*(finished / float(nsims))
    hashes = '#'*int(progress/5.)
    nohash = ' '*int(20-len(hashes))
    print("\r  [{}] {:>3}% {} ".format(hashes+nohash, int(progress), msg), end="")
    sys.stdout.flush()


if __name__ == "__main__":

    ## Parse command line arguments
    args = parse_command_line()

    ## Implicit space is the only model we have at this time
    data = implicit_space.implicit_space(K=args.K, colrate=args.colrate, quiet=args.quiet)

    ## Set model parameters
    data.set_metacommunity(args.meta)
    data.prepopulate(mode=args.mode)

    for i in range(args.nsims):
        if not i % 10000:
            progressbar(args.nsims, i)
            #print(i, len(data.local_community), len(set(data.local_community)))
        data.step()
    progressbar(100, 100, "done!")

    abundance_distribution = data.get_abundances(octaves=args.octaves)

    plot_abundances_ascii(abundance_distribution)

    sp = data.get_species()
    for s in sp:
        s.simulate_seqs()
        s.get_sumstats()

    ## This is dumb right now cuz it resimulates seq
    print(tabulate_sumstats(data))
