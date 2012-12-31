################################################
#
# 
#
#

import os, sys, pickle, re

inpath = sys.argv[1]
outkern = sys.argv[2]
R = sys.argv[3]

def floatToString(f):
    return "%.3f"%f

def floatToString_long(f):
    return f.__str__()
    #return "%.8f"%f

# set is an array of floats
def calculateAverage(set):
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

def calculateStandardDeviation(set):
    avg = calculateAverage(set)
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

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
# OUTPUT: a PDF plot will be written to "outpath", plotting PP bins versus
# the proportion of sites in the ancestor with that PP.
def Plot_pp_proportions(datpath, outpath):
    fin = open(datpath, "r")
    lines = fin.readlines()
    fin.close()
    
    ppbins = {}
    for i in range(0,21):
        ppbins[i] = []
    
    count_total_sites = 0
    
    for l in lines:
        l = l.strip()
        tokens = l.split()
        state = tokens[1]
        if state != "-":
            pp = float(tokens[2])
            ppbins[ binForProb(pp) ].append( pp )
            count_total_sites += 1

    print "Total number of sites = ", count_total_sites

    sum90 = 0.0
    for i in range(18,21):
        sum90 += ppbins[ i ].__len__()
    
    sum80 = 0.0
    for i in range(16,21):
        sum80 += ppbins[ i ].__len__()
        
    print "sites with PP >= 0.90: ", sum90, sum90 * 1.0 / count_total_sites
    print "sites with PP >= 0.80: ", sum80, sum80 * 1.0 / count_total_sites

    points = {}
    for i in range(0, 21):
        points[ i ] = float(ppbins[ i ].__len__())/ float(count_total_sites)
    
    # create a CRAN script:
    f = open(outpath + ".cran", "w")

    # add points
    xvals = points.keys()
    xvals.sort()
    x = "x <- c("
    for k in xvals:
        k = probForBin(k)
        x += floatToString(k) + ","
    x = re.sub(",$", "", x)
    x += ")"
    f.write( x + "\n")
    y = "y <- c("
    for k in xvals:
        y += floatToString(points[k])  + ","
    y = re.sub(",$", "", y)
    y += ")"
    f.write( y + "\n")

    f.write("pdf(\"" + outpath +".pdf" + "\", width=6, height=4) \n")
    f.write("barplot(y, names.arg=x, horiz=FALSE, xlab=\"ASR posterior probability\", ylab=\"proportion of sites\", ylim=c(0.0, 1.0) )\n")  
    f.write("dev.off()\n")
    f.close()
    
    os.system(R + " --save < " + outpath + ".cran > " + outpath + ".cran_log") 

Plot_pp_proportions(inpath, outkern)