## gimmeSAD - Joint modelling of abundance and genetic diversity. 

gimmeSAD is an integrated model of population genetics and community ecology. It allows one to model abundance in an ecological community under a neutral process while simultaneously modelling neutral population genetic variation.

## Quick install
After installing conda for python2.7, you can create a new conda env and install gimmeSAD inside it:
```
conda create -n py27-gimmeSAD python=2.7
conda activate py27-gimmeSAD
conda install -c iovercast gimmesad
```

## Basic usage
In the basic usage of `gimmeSAD` the most fundamental arguments are the size of the local community (`-k`, which is measured in number of individuals), and the colonization rate (`-c`, which is the probability per timestep of a colonization event).  

At any time the `-h` flag can be passed to the `gimmeSAD` command to display copious usage information, as well as a few example command line uses. 

```
usage: gimmeSAD [-h] [-a] [-c colrate] [-C colonizers] [-k K] [-m meta]
                [-n nsims] [-o outdir] [--empirical_dir empiricaldir]
                [--model model] [-p mode] [-r recording_period]
                [-i invasion_time] [-I invasiveness] [-e] [-b bottleneck] 
                [--do_plots] [--curves] [--plot_models] [-q] [-v] [-f]
                [--version]               
optional arguments:                                                                                                                                                       
  -h, --help            show this help message and exit
  -a                    Print species abundances in octaves
  -c colrate            Set colonization rate
  -C colonizers         Switch mode to clustered colonization and set the # of
                        colonizers per event
  -k K                  Carrying capacity of the island (max # individuals in
                        the local community
  -m meta               Source metacommunity from file or generate uniform
  -n nsims              Number of demographic events to simulate
  -o outdir             Output directory for log/data files and pngs
  --empirical_dir empiricaldir             
                        Directory with the empirical data
  --model model         Data model for writing output files.
  -p mode               Select mode for prepopulating the island
                        (volcanic/landbridge)
  -r recording_period   Length of timeslice between samples for logging
  -i invasion_time      Timestep of invasion of invasive species
  -I invasiveness       Invasiveness of the invasive species
  -e                    Do exponential growth, otherwise constant
  -b bottleneck         Strength of the bottleneck
  --do_plots            Generate a bunch of plots and animations. Default is
                        not to do plots.
  --curves              Plot rank abundance as curves rather than points
  --plot_models         Add expectations under different statistical models to
                        the rank abundance plot
  -q, --quiet           do not print to stderror or stdout.
  -v, --verbose         print lots of debug information
  -f, --force           force overwrite of existing data
  --version             show program's version number and exit

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
    gimmeSAD -v                         ## Be very chatty
```

## License
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)  
gimmeSAD is licensed under under the Creative Commons Attribution-Share Alike license.

## Citing gimmeSAD
The manuscript describing the simulation method and an approximate Bayesian computation framework for parameter estimation may be cited as:

`Overcast, I., Emerson, B., Hickerson, M. "An integrated model of population genetics and community ecology". (In Review).`

