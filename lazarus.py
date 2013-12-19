#!/usr/bin/python

######################################################
#
# Lazarus
# is an ancestral reconstruction tool,
# which parses the output from PAML,
# parsimoniously places gaps at ancestral sites,
# and integrates reconstructions from a distribution of trees.
#
# 
# Victor Hanson-Smith
# victorhs@cs.uoregon.edu
#
#
# Required command-line arguments:
# --alignment <path to FASTA alignment>
# --tree <path to newick tree file>
# --model <path to .dat substitution matrix>
# --outputdir <path to desired output directory>
#
# Optional command-line arguments:
# --branch_lengths <value> // where value = 'estimate' or 'fixed'
# --asrv <value> // where value = the number of gamma categories, or '0' for no among site rate variation]
# --alpha <value> // the alpha value for the gamma model
# --fix_asrv <True/False> // if True, then use the fixed alpha value.  If False, then find the Ml estimate of the value.
# --codeml // to run amino acid analysis
# --baseml // to run nucleotide analysis
# --verbose
# --gapcorrect // parsimoniously place gaps at ancestral sites, according to Fitch's parsimony.
# --outgroup <outgroup>  // where outgroup is formatted as [A,B,C] (use brackets, and no spaces),
#                        // and A, B, and C are names of extant taxa.
# --getanc <value> // where value is an ancestral node number on the ML tree.
#                  // alternatively, you can specify --getanc --ingroup [A,B,C], where "[A,B,C]" is the bracket-bound list of ingroup taxa with no spaces.
#
# --topoprior <value> // where the default value is "pp" (to weight ancestral reconstructions by the PP of their tree
#                     // alternatively, you can use "flat" to weight ancestral reconstructions equally
#
# Depricated command-line arguments:
# --postpaml
######################################################
from includeMe import *
from util import *

# Load a previously serialized engine object...
def loadState(path):
    print "--> Loading your Lazarus workspace. . ."
    if os.path.exists( path ) == False:
        print "\nI cannot find any saved workspaces.\nI'm expecting to find a file named " + path + "\n"
        exit(1)
    f = open(path, "r")
    x = pickle.load(f)
    print ". . .this workspace contains " + x.data.trees.__len__().__str__() + " trees."
    return x

#######################################################
#
# "main" starts here:
#
#######################################################
try:
    if sys.argv.__contains__("--outputdir") == False:
        sys.argv.append("--outputdir")
        sys.argv.append( os.getcwd() )    
    
    argParser = ArgParser(sys.argv)
    engine = None
        
    x = argParser.doesContainArg("--codeml")
    if x != False:
        engine = Engine(argParser)
        print "--> Importing your files. . ."
        engine.importFiles()
        print "--> Starting PAML jobs. . ."    
        engine.startPamlJobs()
        print "--> Parsing the results from PAML. . ."    
        engine.parsePamlResults()
        print "--> Saving your Lazarus workspace. . ."
        engine.saveState()
        if argParser.getOptionalToggle("--cleanup"):
            engine.cleanupPamlResults()
        
    x = argParser.doesContainArg("--baseml")
    if x != False:
        engine = Engine(argParser)
        print "--> Importing your files. . ."
        engine.importFiles(useAminoAcids=False)
        print "--> Starting PAML jobs. . ."    
        engine.startPamlJobs()
        print "--> Parsing the results from PAML. . ."    
        engine.parsePamlResults()
        print "--> Saving your Lazarus workspace. . ."
        engine.saveState()
        if argParser.getOptionalToggle("--cleanup"):
            engine.cleanupPamlResults()
        
    # skip PAML, jump to the post-PAML analysis:
    x = argParser.getOptionalArg("--postpaml")
    if x != False:
        engine = Engine(argParser)
        print "--> Importing your files. . ."
        engine.importFiles()
        print "--> Parsing PAML results. . ."    
        engine.parsePamlResults()
        print "--> Saving your Lazarus workspace. . ."
        engine.saveState()
            
    x = argParser.getOptionalArg("--getanc")
    if x != False:
        if engine == None:
            engine = loadState(argParser.getArg("--outputdir") + "/lazarus_workspace.save")
            engine.argParser = argParser
        engine.getAncestralSequence()          
        if argParser.getOptionalToggle("--cleanup"):
            engine.cleanupPamlResults()
    
except AssertionError:
    print "\nYikes! Something went wrong.\nLazarus is stopping."
    print "Scroll up to find error messages.\n"
    exit(1)

print "\n\nLazarus is done.\n"
