##############################################################################
#
# Node class:
# a class representing a single phylogenetic node
#
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
import os
from sequenceSite import Site

class Node:
    def __init__(self, id):
        self.nodeNumber = id
        self.sites = {} # key = site number, value = a Site object

    def addSite(self, siteNumber):
        self.sites[siteNumber] = Site(siteNumber)
        
        
    def writeSequenceDetail(self, distributionPath):
        self.writeSequenceDetail_withHeader(distributionPath, None)
        
    # 'header' is a special message to append to the top of the file.
    def writeSequenceDetail_withHeader(self, distributionPath, header):
        if os.path.exists(distributionPath):
            os.system("rm " + distributionPath)
        fd = open(distributionPath, "w")
        if header != None:
            fd.write(header + "\n")
        for s in self.sites.keys():
            #f.write(self.sites[s].toString() + "\n")
            x = self.sites[s].toString()
            fd.write(x + "\n")
        fd.close()
     
    #
    # This reads a file, as written by the engine.spewTree method
    #
    def readSequenceDetail(self, filepath):
        f = open(filepath, "r")
        lines = f.readlines()
        for l in lines:
            l = l.strip()
            site = Site(None)
            site.fromString(l)
            self.addSite( site.siteNumber )
            self.sites[ site.siteNumber ] = site     
        f.close()