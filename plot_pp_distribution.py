#########################################
#
# plot_pp_distribution.py
#
# Victor Hanson-Smith
# University of California, San Francisco
# victor.hanson-smith@ucsf.edu
# victorhansonsmith@gmail.com
#
#
# This script will print summary statistics about the posterior probability distribution for a reconstructed
# ancestral sequence.
#
# USAGE (from the command line):
# %> python plot_pp_distribution.py X.dat
#
# . . . where X.dat is your query ancestral node.
#
################################################################

import math, os, sys, re

def floatToString(f):
    return "%.3f"%f

def floatToString_long(f):
    return f.__str__()

# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )



# returns a bin number for this P value
def binForProb(p):
    return int(p / 0.05)

# return the P value of the floor of this bin
def probForBin(b):
    x = float(b*5) / float(100)
    if x == 1.00:
        return x
    return x + 0.025

# INPUT: a *.dat file, with a distribution of ancestral states for all sites.
# OUTPUT: prints to screen the frequencies of PPs (binned into 5% bins).
def Plot_pp_proportions(datpath):
    fin = open(datpath, "r")
    lines = fin.readlines()
    fin.close()
    
    ppbins = {}
    for i in range(0,21):
        ppbins[i] = []
    
    count_total_sites = 0
    
    allpps = []
    
    for l in lines:
        l = l.strip()
        tokens = l.split()
        state = tokens[1]
        if state != "-":
            pp = float(tokens[2])
            allpps.append(pp)
            ppbins[ binForProb(pp) ].append( pp )
            count_total_sites += 1

    print "\n=========================================="
    print "\n. Posterior probability (PP) summary"
    print "  for ancestor " + datpath
    print "\n\t", count_total_sites, " total sites." 


    sum90 = 0.0
    for i in range(18,21):
        sum90 += ppbins[ i ].__len__()
    
    sum80 = 0.0
    for i in range(16,21):
        sum80 += ppbins[ i ].__len__()
        
    print "\t", int(sum90), "(%.3f%%)"%(sum90 * 1.0 / count_total_sites * 100), "sites have PP >= 0.90 "
    print "\t", int(sum80), "(%.3f%%)"%(sum80 * 1.0 / count_total_sites * 100), "sites have PP >= 0.80 "

    print "\n\tmean PP = \t\t\t%.3f"%mean(allpps)
    print "\tstandard deviation = \t\t%.3f"%sd(allpps)
    print "\tvariance = \t\t\t%.3f"%var(allpps)
    print "\tstandard error of the mean = \t%.3f"%stderr(allpps)


    points = {}
    for i in range(0, 21):
        points[ i ] = float(ppbins[ i ].__len__())/ float(count_total_sites)
    
    print "\n. PP Distribution:\n"
    print "\tPP range\tproportion of sites"
    print "\t------------\t-------------------"
    for i in range(0, 21):
        if i < 20:
            print "\t[%.2f"%(probForBin(i) - 0.025) + ", %.2f)"%(probForBin(i) + 0.025), "\t%.3f"%points[i]
        else:
            print "\t[1.00]\t\t", "%.3f"%points[i]

Plot_pp_proportions(sys.argv[1])