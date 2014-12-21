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

#
# Identify restriction sites.
#
rsites = []
if sys.argv[2:].__len__() > 0 and sys.argv[2:].__len__() < 3:
    if sys.argv[2].isdigit() and sys.argv[3].isdigit():
        limstart = int(sys.argv[2])
        limstop = int(sys.argv[3])
        for i in range(limstart, limstop+1):
            rsites.append(i)
            #print i
else:
    for i in sys.argv[2:]:
        #print i
        rsites.append( int(i) )

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
    mlseq = ""
    fin = open(datpath, "r")
    lines = fin.readlines()
    fin.close()
    
    count_total_sites = 0
    
    ppbins = {}
    for i in range(0,21):
        ppbins[i] = []
        
    allpps = []
    
    for l in lines:
        if l.__contains__("#") or l.__len__() < 2:
            continue
        #print l
        l = l.strip()
        tokens = l.split()
        site = int(tokens[0])
#        if limstart != False:
#            if site < limstart:
#                continue
#        if limstop != False:
#            if site > limstop:
#                continue
        if rsites.__len__() > 0 and site not in rsites:
            continue
        
#        if limstart != False:
#            if site < limstart:
#                continue
#        if limstop != False:
#            if site > limstop:
#                continue
        state = tokens[1]
        mlseq += state
        if state != "-":
            #print site, state
            pp = float(tokens[2])
            allpps.append(pp)
            ppbins[ binForProb(pp) ].append( pp )
            count_total_sites += 1
    print "\n=========================================="
    print ". ML sequence =", mlseq
    print "\n"

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

    sum70 = 0.0
    for i in range(14,21):
        sum70 += ppbins[ i ].__len__()
        
    print "\t", int(sum90), "(%.3f%%)"%(sum90 * 1.0 / count_total_sites * 100), "sites have PP >= 0.90 "
    print "\t", int(sum80), "(%.3f%%)"%(sum80 * 1.0 / count_total_sites * 100), "sites have PP >= 0.80 "
    print "\t", int(sum70), "(%.3f%%)"%(sum70 * 1.0 / count_total_sites * 100), "sites have PP >= 0.70 "

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
    return [ppbins,count_total_sites,allpps]

def R_barplot(data):
    xlab = "PP"
    ylab = "P"
    filekeyword = "anc"
    
    pointsets = data.keys()
    pointsets.sort()
    
    finalset = pointsets[ pointsets.__len__()-1 ]
    tablepath = "barplot.table." + filekeyword + ".txt"
    fout = open(tablepath, "w")
    for p in pointsets:
        if p != finalset:
            fout.write(p.__str__() + "\t")
        else:
            fout.write(p.__str__() )
    fout.write("\n")
    for p in pointsets:
        if p != finalset:
            fout.write( data[p].__str__() + "\t")
        else:
            fout.write( data[p].__str__() )            
    fout.write("\n")
    fout.close()
    
            
    pdfpath = "barplot." + filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=6, height=3);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"
    
    #print pointsets


    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += (p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, col=\"black\", ylim=range(0,1.0), names.arg=pointsets);\n"
    cranpath = "barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    #os.system("r --no-save < " + cranpath)


def R_lineplot(allpps):    
    cranpath = "anc.seqplot.cran"
    cranout = open(cranpath, "w")
    x = "x <- c("
    y = "y <- c("

    miny = 0.0
    maxy = 1.0
    minx = 1
    maxx = allpps.__len__()

    print allpps
    print minx
    print maxx

    for s in range(minx, maxx+1):
        x += s.__str__() + ","
        y += allpps[s-1].__str__() + ","
    x = re.sub(",$", "", x)
    x += ")"
    cranout.write( x + "\n")
    y = re.sub(",$", "", y)
    y += ")"
    cranout.write( y + "\n")
        

    cranout.write("pdf('" + cranpath + ".pdf', width=6, height=3);\n")
    cranout.write("plot(x,y,xlab='sequence sites',ylab='P of ML state', type='l', lwd=2);\n")
    cranout.write("dev.off();\n")
    cranout.close()
    #os.system("r --no-save < " + cranpath)
    return cranpath

[ppbins,count,allpps] = Plot_pp_proportions(sys.argv[1])
pp_p = {}
for pp in ppbins:
    pp_p[pp] = float(ppbins[pp].__len__())/count
R_barplot(pp_p)
R_lineplot(allpps)