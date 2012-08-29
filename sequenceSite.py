##############################################################################
#
# Site class:
# a class representing a single sequence site and its probability distribution
# of states
#
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
class Site:
    def __init__(self, id):
        self.siteNumber = id
        self.stateProbs = {} # key = a character, value = a posterior probability value

    def addState(self, character, postProb):
        self.stateProbs[character] = postProb
        
    def toString(self):
        """Returns a single line which can be pasted into the nodeX.dat ancestral posterior distribution files."""
        x = []  
        stateProbsCopy = self.stateProbs.copy() # we copy so as to be non-destructive to the original data structure
#        if stateProbsCopy.__contains__("-"):
#            line = self.siteNumber.__str__() + "  - (NA)"
#            return line 
        while stateProbsCopy.__len__() > 0:
            maxp = 0.0
            maxc = ""
            # find max prob:
            for s in stateProbsCopy.keys():
                if float(stateProbsCopy[s]) > maxp:
                    maxp = stateProbsCopy[s]
                    maxc = s
            if maxp == 0.00:
                stateProbsCopy.clear()
            else:
                x.append([maxc,maxp])
                stateProbsCopy.__delitem__(maxc)
        line = self.siteNumber.__str__() + " "
        for i in x:
            c = i[0]
            p = i[1]
            if p == 0.00:
                continue
            line += " "
            line += c.__str__()
            if p > 100.0:
                line += " (NA)"
            elif p > 1.0 and p < 2.0:
                line += " 1.000"
            else:
                line += " %.3f" % p   
        return line

    
    #
    # useful for postorder and preorder traversal in engine.
    #
    def find_ML_state(self):
        stateProbsCopy = self.stateProbs.copy() # we copy so as to be non-destructive to the original data structure     
        maxp = 0.0
        maxc = ""
        # find max prob:
        for s in stateProbsCopy.keys():
            if float(stateProbsCopy[s]) > maxp:
                maxp = stateProbsCopy[s]
                maxc = s
        return maxc

    
    def fromString(self, line):
        tokens = line.split()
        self.siteNumber = int(tokens[0])
        i = 1
        while( i+1 < tokens.__len__() ):
            if False == (tokens[i+1] == "(NA)"):
               tokens[i+1] = float(tokens[i+1])     
            self.addState( tokens[i].__str__(), tokens[i+1] )
            i = i+2
        
        