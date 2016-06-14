#!/usr/bin/env python2.7

from __future__ import print_function
from ascii_graph import Pyasciigraph
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
import subprocess
import collections
import argparse
import datetime
import shutil
import time
import sys
import os

import implicit_BI
import implicit_CI

# pylint: disable=C0103
# pylint: disable=R0903

## $ python -m cProfile -s cumtime lwn2pocket.py

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
    msg = "\n\n"
    for line in  graph.graph('Simulated Species Abundance Distribution', abundance_distribution):
        msg = msg + line + "\n"
    msg = msg + "###############################################################################\n"
    return msg


def plot_abundances_gui(abundance_distribution):
    plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    plt.show()


def tabulate_sumstats(data):
    sp = data.get_species()
    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Col time", "Loc Abund", "Meta Abund", "pi", "pi_net", "Dxy",  "S", "S_island", "pi_island", "S_meta", "pi_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, s.meta_abundance, s.pi, s.pi_net, s.dxy, s.S, s.S_island, s.pi_island, s.S_meta, s.pi_meta] for s in sp]
    return tabulate(acc, headers, floatfmt=".4f")


## This actually is doing pi x dxy, but some of the variable
## names are goofy cuz i developed it for pi x pi_w_island 
def heatmap_pi_dxy(data, write="", title=""):
    """ The write flag specifies whether or not to write the image
    to a file. You must pass in a file name. If this file exists
    it gets overwritten. Normally write should be a full path.
    No need to pass in the file extension."""
    sp = data.get_species()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.dxy, x.pi_island) for x in sp])
    max_pi = max([x[0] for x in pis])
    max_pi_island = max([x[1] for x in pis])

    ## Make the bins
    pi_bins = np.linspace(0, max_pi, 10)
    pi_island_bins = np.linspace(0, max_pi_island, 10)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for pi, pi_island in pis:
        count_pi = 0
        count_pi_island = 0
        while not pi <= pi_bins[count_pi]:
            count_pi += 1
        while not pi_island <= pi_island_bins[count_pi_island]:
            count_pi_island += 1
        ## increment the heatmap point this corresponds to
        heat[count_pi][count_pi_island] += 1
    plt.pcolor(heat,cmap=plt.cm.Reds)
    plt.xlabel('Dxy')
    plt.ylabel('Pi_w Island')
    plt.colorbar()
    plt.xticks(np.arange(len(pi_bins)), ["{0:.4f}".format(x) for x in pi_bins], rotation='vertical')
    plt.yticks(np.arange(len(pi_bins)), ["{0:.4f}".format(x) for x in pi_island_bins])
    
    ## If writing to a file, don't bother displaying it, plus it hangs the program
    if write:
        plt.title(title)
        plt.savefig(write+".png")
    else:
        plt.show()
    plt.close()


def normalized_pi_dxy_heatmaps(outdir, sp_through_time, equilibria):
    """ Normalize x and y axes for the heatmaps. Only take into account extant species.
    Inputs are the output directory to write to and an ordered dict 
    of the species at every recording duration timepoint """

    print("Generating pi x dxy heatmap animation")
    ## Make the output directory for heatmaps inside the top level output directory
    heat_out = os.path.join(outdir, "normalized_heatmaps")
    if not os.path.exists(heat_out):
        os.mkdir(heat_out)

    equilibria = equilibria.values()
    print(equilibria)

    ## GET MAX pi and dxy
    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.uuid[0] for x in sp_through_time.values()[-1]]
    ## find the max_pi and max_dxy for all extant species through all timepoints
    max_pi_island = 0
    max_dxy = 0
    ## For each recorded timeslice
    my_dxys = []
    my_pi_islands = []
    for sp_list in sp_through_time.values():
        ## Get 
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
        ## pis will be empty if this timeslice includes no extant species
        if pis.any():
            my_max_dxy = max([x[0] for x in pis])
            my_max_pi_island = max([x[1] for x in pis])
            if max_dxy < my_max_dxy:
                max_dxy = my_max_dxy
            if max_pi_island < my_max_pi_island:
                max_pi_island = my_max_pi_island
            my_dxys.append(my_max_dxy)
            my_pi_islands.append(my_max_pi_island)

    print("Got\tmax_dxy - {}\t max_pi_island - {}".format(max_dxy, max_pi_island))
    max_dxy = np.average(my_dxys)
    max_pi_island = np.average(my_pi_islands)
    print("Got\tmax_dxy - {}\t max_pi_island - {}".format(max_dxy, max_pi_island))
    ## Make the heatmaps, one for each timeslice
    ## TODO: Should we include all species or just the ones that live to the end?
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the heatmap
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))
    print(nslices, file_index)
    for i, sp_list in enumerate(sp_through_time.values()):
        title = "Time_"+str(file_index+i)
        write = os.path.join(heat_out, title)
        #print("Doing", title)
        
        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list])

        ## Empty heatmap we'll write into
        heat = np.zeros((10,10), dtype=np.int)

        ## Make the bins
        dxy_bins = np.linspace(0, max_dxy, 20)
        pi_island_bins = np.linspace(0, max_pi_island, 20)

        ## Now you have the bins each value belongs in, but you need to 
        ## go through and populate the heat matrix
        for dxy, pi_island in pis:
            count_dxy = 0
            count_pi_island = 0
            try:
                while not dxy <= dxy_bins[count_dxy]:
                    count_dxy += 1
                while not pi_island <= pi_island_bins[count_pi_island]:
                    count_pi_island += 1
                ## increment the heatmap point this corresponds to
                heat[count_dxy][count_pi_island] += 1
            except Exception as inst:
                ## Got a value bigger than our current max pi/dxy. ignore.
                pass
        plt.pcolor(heat,cmap=plt.cm.Reds)
        plt.xlabel('Dxy')
        plt.ylabel('Pi_w Island')
        tix = np.linspace(0, 30, 6)
        plt.colorbar(ticks=tix)
        plt.xticks(np.arange(len(dxy_bins)), ["{0:.4f}".format(x) for x in dxy_bins], rotation='vertical')
        plt.yticks(np.arange(len(pi_island_bins)), ["{0:.4f}".format(x) for x in pi_island_bins])

        #plt.title(title)
        print("Doing eq - {}".format(equilibria[i]))
        plt.title(title + "\n%equilibrium = {}".format(equilibria[i]))
        plt.savefig(write+".png")
        plt.close()

    ## Do the imagemagick conversion, if possible
    ## `convert -delay 100 outdir/* anim.gif`
    cmd = "convert -delay 100 "\
            + heat_out + "/*.png "\
            + outdir + "/pi_dxy_anim.gif"
    try:
        subprocess.Popen(cmd.split())
    except Exception as inst:
        print("Trouble creating pi x dxy animated gif - {}".format(inst))


def heatmap_pi_dxy_ascii(data, labels=False):
    """ This is kind of a toy. It logs to a file in the output directory
    at recording_period interval, but the axes aren't normalized, so it's
    more semi-informational. """
    sp = data.get_species()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.dxy, x.pi_island) for x in sp])
    ## Set a reasonable default
    max_pi = max_pi_island = 0.1
    if pis.any():
        max_pi = max([x[0] for x in pis])
        max_pi_island = max([x[1] for x in pis])
        #print(max_pi, max_pi_island)

    ## Make the bins
    pi_bins = np.linspace(0, max_pi, 10)
    pi_island_bins = np.linspace(0, max_pi_island, 10)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for pi, pi_island in pis:
        count_pi = 0
        count_pi_island = 0
        while not pi <= pi_bins[count_pi]:
            count_pi += 1
        while not pi_island <= pi_island_bins[count_pi_island]:
            count_pi_island += 1
        ## increment the heatmap point this corresponds to
        heat[count_pi][count_pi_island] += 1

    ret = ""
    ## ascii format the data, and make a weak attempt to convey some information
    if labels:
        ret += "\nDxy\n"

    for i, row in enumerate(heat):
        if labels:
            ret += "{0:.4f}".format(pi_bins[i])
        ret += str(row) + "\n"

    if labels:
        ret += "\t" + str(["{0:.4f}".format(x) for x in pi_island_bins])
        ret += "\n\t\tpi_w island"

    return ret


def write_sizechanges(outdir, yoyo):
    """ Write out sizechanges through time for all species extant
    at the end of the simulation. Takes in an ordered dict of collections """

    out = open(os.path.join(outdir, "pop_yoyos.txt"), 'w')
    keys = yoyo.keys()
    keys.reverse()
    ## Get a list of all extant species, the only ones we care about pop size change in.
    extant = yoyo[keys[0]].keys()
    sizes = {}
    for x in extant:
        sizes[x] = []
    ## For each time slice going back in time, record the current pop size of each 
    ## species of interest
    for k in keys:
       for sp, siz in yoyo[k].items():
            if sp in extant:
                sizes[sp].append(siz)
    for k, v in sizes.items():
        out.write("{} - {}\n".format(k, v))
    

def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 

    gimmeSAD -n 100000                  ## Run the simulation for 100000 steps

    gimmeSAD -o outdir                  ## Directory where tons of junk will get written
    gimmeSAD -r recording_period        ## Num steps to sample and write out after

    gimmeSAD -m metacommunity_LS4.txt   ## Source metacommunity from a file of ints
                                        ## the other option is to pass "uniform"
                                        ## which will generate a random metacommunity
                                        ## with all species at uniform frequency

    gimmeSAD -p volcanic                ## Set the mode for prepopulating (or not)
                                        ## the local community.

    gimmeSAD -i 4                       ## Set # of colonizing inds per colonization event
                                        ## if unset basic immigration is assumed (1 individual)
    gimmeSAD -c 0.001                   ## Set colonization rate

    gimmeSAD -a                         ## Output abundances in octaves
    gimmeSAD -q                         ## Don't be so chatty

   """)

    ## add model arguments 
    parser.add_argument('-n', metavar='nsims', dest="nsims", type=int,
        default=50000,
        help="Number of demographic events to simulate")

    parser.add_argument('-k', metavar='K', dest="K", type=int,
        default=10000,
        help="Carrying capacity of the island (max # individuals in the local community")

    parser.add_argument('-i', metavar='colonizers', dest="colonizers", type=int,
        default=0,
        help="Switch mode to clustered colonization and set the # of colonizers per event")

    parser.add_argument('-r', metavar='recording_period', dest="recording_period", type=int,
        default=0,
        help="Length of timeslice between samples for logging")

    parser.add_argument('-m', metavar='meta', dest="meta", type=str,
        default="metacommunity_LS4.txt",
        help="Source metacommunity from file or generate uniform")

    parser.add_argument('-c', metavar='colrate', dest="colrate", type=float,
        default=0.003,
        help="Set colonization rate")

    parser.add_argument('-o', metavar='outdir', dest="outdir", type=str,
        default="output",
        help="Output directory for log/data files and pngs")

    parser.add_argument('-p', metavar='mode', dest="mode", type=str,
        default="volcanic",
        help="Select mode for prepopulating the island (volcanic/landbridge)")

    parser.add_argument('-a', dest="octaves", action='store_true',
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

    if args.colonizers:
        ## Implicit space and clustered immigration
        data = implicit_CI.implicit_CI(K=args.K, colrate=args.colrate, mig_clust_size=args.colonizers, quiet=args.quiet)
    else:
        ## Implicit space, one colonizer per event
        data = implicit_BI.implicit_BI(K=args.K, colrate=args.colrate, quiet=args.quiet)

    ## Set model parameters
    data.set_metacommunity(args.meta)
    data.prepopulate(mode=args.mode)

    ## Setup output directory and files
    try:
        if os.path.exists(args.outdir) and args.force:
            shutil.rmtree(args.outdir)
            os.mkdir(args.outdir)
        elif not os.path.exists(args.outdir):
            os.mkdir(args.outdir)
        else:
            sys.exit("Output directory exists - {}\nUse the force flag -f to overwrite".format(args.outdir))
        out = open(os.path.join(args.outdir, "pi_x_dxy.log"), "w")
        yoyofile = open(os.path.join(args.outdir, "sizechange_through_time.log"), 'w')
    except Exception as inst:
        sys.exit("problem opening output for writing - {}\n{}".format(args.outdir, inst))

    ## If you don't specify a recording period, then just record
    ## 10% of the time
    if not args.recording_period:
        if args.nsims < 1:
            args.recording_period = 100000
        else:
            args.recording_period = args.nsims/10

    ## Start the main loop
    start = time.time()
    elapsed = 0
    reached_equilib = False

    ## Track pop size changes through time
    ## Also track species sumstats through time so we can normalize the
    ## heatmap plots
    yoyo = collections.OrderedDict()
    sp_through_time = collections.OrderedDict()
    equilibria = collections.OrderedDict()

    ## if args.nsims == -1 just run until double equilibrium
    ## or 10e7 steps (effectively forever)
    if args.nsims < 1:
        args.nsims = 1000000000
    #for i in range(1, args.nsims):
    i = 0
    while True:
        if i == args.nsims:
            break
        i += 1
        data.step()

        ## Print the progress bar every once in a while
        ## Set a flag for equilibrium. If you've reached it, flip all the
        ## founder flags back to True and keep running til next equilibrium
        ## then stop
        percent_equil = 0
        if not i % 10000:
            ## Record abundance of each species through time
            ## Make a counter for the local_community, counts the number of
            ## individuals w/in each species
            abundances = collections.Counter([x[0] for x in data.local_community])
            yoyofile.write(str(abundances)+"\n")
            yoyo[i] = abundances

            ## Test for equilibrium
            founder_flags = [x[1] for x in data.local_community]
            percent_equil = float(founder_flags.count(False))/len(founder_flags)            
            if (not any(founder_flags)) and reached_equilib:
                print("\nReached second equilibrium")
                break
            elif not any(founder_flags):
                print("\nReached first equilibrium, reset founder flags")
                data.local_community = [(x[0], True) for x in data.local_community]
                reached_equilib = True

            ## Update progress bar
            secs = time.time()-start
            elapsed = datetime.timedelta(seconds=secs)
            ## Calculate the remaining time
            rate = float(i)/secs
            remaining_secs = (args.nsims - i) / rate
            progressbar(args.nsims, i, " | %equilib - {} | elapsed - {} | remaining - {}".format(percent_equil, elapsed, datetime.timedelta(seconds=int(remaining_secs))))

        ## Recording data every once in a while
        if not i % args.recording_period and not i == 0: 
            ## Every once in a while write out useful info
            data.simulate_seqs()
            sp_through_time[i] = data.get_species()
            equilibria[i] = percent_equil
            out.write("step {}\n".format(i))
            out.write(heatmap_pi_dxy_ascii(data, labels=False)+"\n")
            ## Unnormalized axes version of heatmap. Don't use.
            #heatmap_pi_dxy(data, os.path.join(args.outdir, "plot-"+str(i)+".png"), title="time = " + str(i))

    progressbar(100, 100, "  |  {}  steps completed  |  Total runtime   {}".format(i, elapsed))

    if not reached_equilib:
        founder_flags = [x[1] for x in data.local_community]
        percent_equil = float(founder_flags.count(False))/len(founder_flags)
        print("How close to equilibrium? {}".format(percent_equil))
    
    ## When finished simulate the final set of sequences
    data.simulate_seqs()
    sp_through_time[i] = data.get_species()
    equilibria[i] = percent_equil
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))

    ## Print out some informative business
    ## Get all results and write out final sumstats
    abundance_distribution = data.get_abundances(octaves=args.octaves)
    print(plot_abundances_ascii(abundance_distribution))
    if not args.quiet:
        print(tabulate_sumstats(data))

    ## Write out to log files
    write_sizechanges(args.outdir, yoyo)
    with open(os.path.join(args.outdir, "gimmeSAD.out"), "w") as stats:
        stats.write("Parameters - {}\n".format(args))
        stats.write("Raw abundance dist - {}\n".format(data.get_abundances(octaves=False)))
        stats.write("Abundance in octaves - {}\n".format(data.get_abundances(octaves=True)))
        stats.write(plot_abundances_ascii(abundance_distribution))
        stats.write("\n")
        stats.write(tabulate_sumstats(data))

    ## Make the normalized pi_x_dxy heatmaps
    normalized_pi_dxy_heatmaps(args.outdir, sp_through_time, equilibria)

