##############################################################################
#
# DataWharehouse class: 
# a class for storing and accessing trees and ancestral sequences.
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
import math
from tree import *

class DataWharehouse:
    def __init__(self, argParser):
        self.argParser = argParser    
        self.trees = {} # key = tree number, value = a Tree object
    
    #
    # add a newick tree to our collection.
    # 'treeID' = the tree number (an integer).  This ID will be used in the future for other data structures.
    #    
    def addTree(self, treeID):
        self.trees[treeID] = Tree(treeID)
    
    #    
    # at the end of this method, the post. prob for each tree X will be stored to trees[X].postprob
    def calculateTreePosteriorProbabilities(self):
        # 0. get the non-excluded trees:
        useTheseTrees = {}
        for t in self.trees:
            if self.trees[t].exclude == False:
                useTheseTrees[t] = self.trees[t]
            else:
                self.trees[t].postprob = 0.0
        
        #print "useTheseTrees lnl = " + self.trees[1].lnL.__str__()
        
        # 1. find the max lnL:
        maxLnl = None
        for t in useTheseTrees:
            if maxLnl == None:
                maxLnl = self.trees[t].lnL
            elif maxLnl < self.trees[t].lnL:
                maxLnl = self.trees[t].lnL

                
        # 2. scale all lnLs by the max lnL and then convert to normal likelihoods (non log):
        scaledLnls = {} # key = tree ID, value = scaled lnL
        for t in useTheseTrees:
            #print "dataWharehouse 50:", t
            #print "dataWharehouse 51:", self.trees[t].lnL
            #print "dataWharehouse 52:", maxLnl
            scaledLnls[t] = math.exp( self.trees[t].lnL - maxLnl )
        
        # 3. sum the scaled lnLs:
        sum = 0
        for t in useTheseTrees:
            sum += scaledLnls[t]
        
        # 4. calculate post probs:
        for t in useTheseTrees:
            self.trees[t].postprob = scaledLnls[t] / sum
        
    #
    # returns a Node object, by reading the NodeX file from disk.
    #    
    def getNode(self, treeid, nodeid):        
        n = Node(nodeid)
        n.readSequenceDetail( self.argParser.getArg("--outputdir") + "/tree" + treeid.__str__() + "/node" + nodeid.__str__() + ".dat" )
        return n
    
    #
    # Returns the tree number of the tree with the highest likelihood
    #
    def getMlTreeNumber(self):
        maxL = None # highest found log likelihood
        maxT = None # the tree associated with maxL
        for t in self.trees:
            if maxL == None:
                maxL = self.trees[t].lnL
                maxT = self.trees[t].treeNumber
            elif maxL < self.trees[t].lnL:
                maxL = self.trees[t].lnL
                maxT = self.trees[t].treeNumber
        return maxT        
    
    #
    # Returns a Node object, which contains the ancestral sequence and its distribution
    #
    # ingroup = a list of ingroup taxa IDs
    def calculateMapAncestralSequence(self, ingroup, outgroup):
        treeAncs = {} # key = tree ID, value = node name of analog ancestor
        for t in self.trees:
            n = self.trees[t].getLCA(ingroup, outgroup)
            if n == False: # does this tree have an anlogous node?
                print "Warning: Excluding tree " + t.__str__() + ", for LSA of ingroup " + ingroup.__str__()
                self.trees[t].exclude = True
            else:
                #print "tree: " + t.__str__() + ", ingroup: " + ingroup.__str__() + ", lca: " + n.Name
                treeAncs[t] = n.Name
                #print "treeAncs[" + t.__str__() + "] = " + n.Name
                            
        self.calculateTreePosteriorProbabilities()
            
        topoprior = self.argParser.getOptionalArg("topoprior")
        if topoprior == False:
            topoprior = "pp"
                
        mapAnc = Node(-1)
        for t in treeAncs:
            n = self.getNode(t, treeAncs[t])
            for s in n.sites:
                if mapAnc.sites.__contains__(s) == False:
                    mapAnc.addSite(s)
                for c in n.sites[s].stateProbs:
                    if c == "-":
                        continue
                    else:
                        if mapAnc.sites[s].stateProbs.__contains__(c) == False:
                            mapAnc.sites[s].addState(c, 0.0)
                        if topoprior == "pp":
                            mapAnc.sites[s].stateProbs[c] += float(self.trees[t].postprob) * float(n.sites[s].stateProbs[c])   
                        elif topoprior == "flat":
                            mapAnc.sites[s].stateProbs[c] += float(1 / self.trees.keys().__len__()) * float(n.sites[s].stateProbs[c]) 
        return mapAnc
             
        