##############################################################################
#
# Lazarus Unit Tests
#
# This test file is deprecated!
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
from includeMe import *

class testEngine:
    
    def __init__(self):
        pass
        pass
    
    def test_importFiles(self):
        print "\n--> test_importFiles:"
        
        #
        # test 0: we attempt to import non-existent files:
        #
        sys.stdout.write("\n. test 0: this test is intended to throw an (ERROR); please ignore it . .\n")
        sys.argv.append("--alignment")
        sys.argv.append("doesnotexist.phy")
        sys.argv.append("--tree")
        sys.argv.append("doesnotexist.tre")
        sys.argv.append("--model")
        sys.argv.append("doesnotexist.dat")
        
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        foundError = False
        try:
            self.engine.importFiles() 
        except AssertionError:
            print ". . .SUCCESS."
            foundError = True
        if foundError == False:
            print "failed!"      
            print "This test was intended to throw an (ERROR), but that didn't happen."
            raise AssertionError
                                
        #
        # test 1a: our tree contains a taxon which is not in our alignment
        #
        sys.stdout.write("\n. test 1a: this test is intended to throw an (ERROR); please ignore it . .\n")
        
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA,taxonB,taxonC,taxonE);\n")
        ft.close()
        
        ffasta = open(self.outDir + "/testalign.fasta", "w")
        ffasta.write(">taxonA\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonB\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonC\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonD\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write("")
        ffasta.close()

        sys.argv[sys.argv.index("doesnotexist.phy")] = self.outDir + "/testalign.fasta"
        sys.argv[sys.argv.index("doesnotexist.tre")] = self.outDir + "/testtrees.tree"
        sys.argv[sys.argv.index("doesnotexist.dat")] = "./paml/dat/jones.dat"      

        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        foundError = False
        try:
            self.engine.importFiles() 
        except AssertionError:
            print ". . .SUCCESS."
            foundError = True
        if foundError == False:
            print "failed!"      
            print "This test was intended to throw an (ERROR), but that didn't happen."
            raise AssertionError
        
        #
        # test 1b: our alignment contains a taxon which is not in the trees.
        #
        sys.stdout.write("\n. test 1b: this test is intended to throw an (ERROR); please ignore it . .\n")
        
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        ffasta = open(self.outDir + "/testalign.fasta", "w")
        ffasta.write(">taxonA\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonB\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonC\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonD\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonE\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write("")
        ffasta.close()      
        
        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        foundError = False
        try:
            self.engine.importFiles() 
        except AssertionError:
            print ". . .SUCCESS."
            foundError = True
        if foundError == False:
            print "failed!"      
            print "This test was intended to throw an (ERROR), but that didn't happen."
            raise AssertionError

        #
        # test 2: our alignment contains sequences with different lengths
        #
        sys.stdout.write("\n. test 2: this test is intended to throw an (ERROR); please ignore it . .\n")
        
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        ffasta = open(self.outDir + "/testalign.fasta", "w")
        ffasta.write(">taxonA\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonB\n")
        ffasta.write("AKKKFVK\n") # this is the short sequence
        ffasta.write(">taxonC\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonD\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write("")
        ffasta.close()     

        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        foundError = False
        try:
            self.engine.importFiles() 
        except AssertionError:
            print ". . .SUCCESS."
            foundError = True
        if foundError == False:
            print "failed!"      
            print "This test was intended to throw an (ERROR), but that didn't happen."
            raise AssertionError

        #
        # test 3a: we try to import correctly-formed input
        #
        sys.stdout.write("\n. test 3a. . .")

        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        ffasta = open(self.outDir + "/testalign.fasta", "w")
        ffasta.write(">taxonA\n")
        ffasta.write("ARKKFV-Y\n")
        ffasta.write(">taxonB\n")
        ffasta.write("AKKRFVKY\n")
        ffasta.write(">taxonC\n")
        ffasta.write("AKKLFVKY\n")
        ffasta.write(">taxonD\n")
        ffasta.write("AK-KFVDY\n")
        ffasta.write("")
        ffasta.close()       

        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        
        sys.argv[sys.argv.index("./paml/dat/jones.dat")] = "./paml/dat/jones.dat"
        self.engine = Engine(self.argParser)
        
        self.engine.importFiles() 

        expectedTrees = [1,2,3]
        for t in expectedTrees:
            pj = self.engine.pamljobs[t]
            if os.path.exists(pj.treePath) == False:
                print "failed!"      
                print "The following file does not exist: " + pj.treePath
                raise AssertionError
            if os.path.exists(pj.alignmentPath) == False:
               print "failed!"      
               print "The following file does not exist: " + pj.alignmentPath
               raise AssertionError
            if os.path.exists(pj.modelPath) == False:
                print "failed!"
                print "The following files does not exist: " + pj.modelPath
                raise AssertionError
        print "SUCCESS."
        
        #
        # test 3b: we'll import files with relative paths, instead of absolute paths
        #
        sys.stdout.write("\n. test 3b:. . .")
        
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        ffasta = open(self.outDir + "/testalign.fasta", "w")
        ffasta.write(">taxonA\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonB\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonC\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write(">taxonD\n")
        ffasta.write("AKKKFVKY\n")
        ffasta.write("")
        ffasta.close()

    
        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        self.engine.importFiles() 

        expectedTrees = [1,2,3]
        for t in expectedTrees:
            pj = self.engine.pamljobs[t]
            if os.path.exists(pj.treePath) == False:
                print "failed!"      
                print "The following file does not exist: " + pj.treePath
                raise AssertionError
            if os.path.exists(pj.alignmentPath) == False:
               print "failed!"      
               print "The following file does not exist: " + pj.alignmentPath
               raise AssertionError
            if os.path.exists(pj.modelPath) == False:
                print "failed!"
                print "The following files does not exist: " + pj.modelPath
                raise AssertionError
        print "SUCCESS."
        
        #
        # test 4: same as test #3, but with different file formats
        #
        sys.stdout.write("\n. test 4. . .")

        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        fphy = open(self.outDir + "/testalign.phy", "w")
        fphy.write("4 10\n")
        fphy.write("taxonA  ARKKFV-Y\n")
        fphy.write("taxonB  AKKRFVKY\n")
        fphy.write("taxonC  AKKLFVKY\n")
        fphy.write("taxonD  AK-KFVDY\n")
        fphy.write("")
        fphy.close()     

        # (we need to reset the engine, to wipe away the results from test #1)
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        self.engine.importFiles() 

        expectedTrees = [1,2,3]
        for t in expectedTrees:
            pj = self.engine.pamljobs[t]
            if os.path.exists(pj.treePath) == False:
                print "failed!"      
                print "The following file does not exist: " + pj.treePath
                raise AssertionError
            if os.path.exists(pj.alignmentPath) == False:
               print "failed!"      
               print "The following file does not exist: " + pj.alignmentPath
               raise AssertionError
            if os.path.exists(pj.modelPath) == False:
                print "failed!"
                print "The following files does not exist: " + pj.modelPath
                raise AssertionError
        print "SUCCESS."
        
        #
        # test 5: nucleotide data
        #
        sys.argv[sys.argv.index("./paml/dat/jones.dat")] = "JC69"
        sys.stdout.write("\n. test 5. . .")
        print sys.argv
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        fphy = open(self.outDir + "/testalign.phy", "w")
        fphy.write("4 8\n")
        fphy.write("taxonA  CTGATCTC\n")
        fphy.write("taxonB  ATCGATCG\n")
        fphy.write("taxonC  ATCGAGTA\n")
        fphy.write("taxonD  TGCATTAT\n")
        fphy.write("")
        fphy.close()     
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        self.engine.importFiles(useAminoAcids=False) 
        
        expectedTrees = [1,2,3]
        for t in expectedTrees:
            pj = self.engine.pamljobs[t]
            if os.path.exists(pj.treePath) == False:
                print "failed!"      
                print "The following file does not exist: " + pj.treePath
                raise AssertionError
            if os.path.exists(pj.alignmentPath) == False:
               print "failed!"      
               print "The following file does not exist: " + pj.alignmentPath
               raise AssertionError
        print "SUCCESS."
        sys.argv[sys.argv.index("JC69")] = "./paml/dat/jones.dat"

   
    def test_startPamlJobs(self):
        #
        # using data from test #4 in test_importFiles. . .
        #
        print "\n--> test_startPamlJobs with AA:"
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        fphy = open(self.outDir + "/testalign.phy", "w")
        fphy.write("4 10\n")
        fphy.write("taxonA  ARKKFV-Y\n")
        fphy.write("taxonB  AKKRFVKY\n")
        fphy.write("taxonC  AKKLFVKY\n")
        fphy.write("taxonD  AK-KFVDY\n")
        fphy.write("")
        fphy.close()     
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        self.engine.importFiles() 
        
        sys.stdout.write("\n. test 1. . .")
        self.engine.startPamlJobs()    
        for t in self.engine.pamljobs:
            if os.path.exists(self.engine.pamljobs[t].controlPath) == False:
                print "failed!"      
                print "I could not find a PAML control file here: " + self.engine.pamljobs[t].controlPath
                raise AssertionError
        print "SUCCESS." 
        
        sys.stdout.write("\n. test 2. . .")
        os.system("mv " + self.engine.pamljobs[1].executionDirectory + "/rst " + self.engine.pamljobs[1].executionDirectory + "/rst.backup")
        flag = self.engine.pamljobs[1].verifyResults()
        if flag != True:
            print ". . .SUCCESS."
        else:
                print "failed!"      
                print "The method engine.pamljobs[...].verifyResults should have reported an error because the file named 'rst' does not exist."
                raise AssertionError                       
        os.system("mv " + self.engine.pamljobs[1].executionDirectory + "/rst.backup " + self.engine.pamljobs[1].executionDirectory + "/rst ")
   
    
    def test_splitRatesAndRst(self):
        print "\n--> test_splitRatesAndRst with nucleotides:"
        # does it work with nucleotide data?
        os.system("rm -rf " + self.outDir + "/*")
        sys.stdout.write("\n. test 1. . .")
        sys.argv[sys.argv.index("./paml/dat/jones.dat")] = "JC69"
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        
        fphy = open(self.outDir + "/testalign.fasta", "w")
        #fphy.write("4 8\n")
        fphy.write(">taxonA\nCTGATCTC\n")
        fphy.write(">taxonB\nATCGATCG\n")
        fphy.write(">taxonC\nATCGAGTA\n")
        fphy.write(">taxonD\nTGCATTAT\n")
        fphy.write("")
        fphy.close()    
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        self.engine.importFiles(useAminoAcids=False)
        self.engine.startPamlJobs()
        self.engine.parsePamlResults()

        sys.stdout.write("\n. test 1a. . .")
        for t in self.engine.data.trees:           
            if os.path.exists(self.outDir + "/tree" + t.__str__() ) == False:
                print "failed!"      
                print "The following directory does not exist: " + self.outDir + "/tree" + t.__str__()
                raise AssertionError
        print "SUCCESS."
        
        sys.stdout.write("\n. test 1b. . .")
        expectedTreeNodes = { 1:["5", "6", "7"], 2:["5", "6", "7"], 3:["5"] }
        for t in expectedTreeNodes:
            for n in expectedTreeNodes[t]:
                if os.path.exists(self.outDir + "/tree" + t.__str__() + "/node" + n.__str__() + ".dat" ) == False:
                    print "failed!"      
                    print "The following file does not exist: " + self.outDir + "/tree" + t.__str__() + "/node" + n.__str__() + ".dat" 
                    raise AssertionError
        print "SUCCESS."
        sys.argv[sys.argv.index("JC69")] = "./paml/dat/jones.dat"
        

        # amino acid data
        ft = open(self.outDir + "/testtrees.tree", "w")
        ft.write("((taxonA:0.1,taxonB:0.1):0.1,(taxonC:0.1,taxonD:0.1):0.1);\n")
        ft.write("((taxonA:0.1,taxonC:0.1):0.1,(taxonB:0.1,taxonD:0.1):0.1);\n")
        ft.write("(taxonA:0.1,taxonB:0.1,taxonC:0.1,taxonD:0.1);\n")
        ft.close()
        fphy = open(self.outDir + "/testalign.phy", "w")
        fphy.write("4 10\n")
        fphy.write("taxonA  ARKKFV-Y\n")
        fphy.write("taxonB  AKKRFVKY\n")
        fphy.write("taxonC  AKKLFVKY\n")
        fphy.write("taxonD  AK-KFVDY\n")
        fphy.write("")
        fphy.close()     
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        self.engine.importFiles() 
        sys.stdout.write("\n. test 2. . .")
        self.engine.startPamlJobs()
        print "\n--> test_splitRatesAndRst:"
        self.engine.logger.verbosityLevel = 9
        print "trees: " + self.engine.data.trees.__str__()
        self.engine.parsePamlResults()
        
        # test 1: did Lazarus write 'treeX' directories for every tree?
        sys.stdout.write("\n. test 2a. . .")
        for t in self.engine.data.trees:
            if os.path.exists(self.outDir + "/tree" + t.__str__() ) == False:
                print "failed!"      
                print "The following directory does not exist: " + self.outDir + "/tree" + t.__str__()
                raise AssertionError
        print "SUCCESS."
        

        # test 2: did Lazarus pick the correct internal nodes?
        sys.stdout.write("\n. test 2b. . .")
        expectedTreeNodes = { 1:["5", "6", "7"], 2:["5", "6", "7"], 3:["5"] }
        for t in expectedTreeNodes:
            for n in expectedTreeNodes[t]:
                if os.path.exists(self.outDir + "/tree" + t.__str__() + "/node" + n.__str__() + ".dat" ) == False:
                    print "failed!"      
                    print "The following file does not exist: " + self.outDir + "/tree" + t.__str__() + "/node" + n.__str__() + ".dat" 
                    raise AssertionError
        print "SUCCESS."

        

        
    def test_saveState(self):
        print "\n--> test_saveState"
        self.engine.saveState()
        
        sys.stdout.write("\n. test 1. . .")
        if os.path.exists(self.argParser.getArg("--outputdir") + "/lazarus_workspace.save") == False:
            print "failed!"
            print "I cannot find your saved state at: " + self.argParser.getArg("--outputdir") + "/lazarus_workspace.save"
            raise AssertionError
        print "SUCCESS."
        
    # this test requires you to first run test_calculateTreePosteriorProbabilities
    def test_writeTreeSummary(self):
        print "\n--> test_writeTreeSummary"
        
        sys.stdout.write("\n. test 1. . .")
        if False == os.path.exists(self.argParser.getArg("--outputdir") + "/tree_summary.txt" ):
            print "failed!"
            print "I did not find the file named " + self.argParser.getArg("--outputdir") + "/tree_summary.txt" 
            raise AssertionError
        print "SUCCESS."
        
        sys.stdout.write("\n. test 2. . .")
        f = open(self.argParser.getArg("--outputdir") + "/tree_summary.txt" , "r")
        lines = f.readlines()
        if self.engine.data.trees.__len__() != lines.__len__() - 1:
            print "failed!"
            print "The file named " + self.argParser.getArg("--outputdir") + "/tree_summary.txt" + " contains too many trees."
            print "I'm expecting to find " + lines.__len__().__str__() + " trees."
            raise AssertionError 
        f.close()
        print "SUCCESS."
        
    def test_calculateTreePosteriorProbabilities(self):
        print "\n--> test_calculateTreePosteriorProbabilities"
        sys.stdout.write("\n. test 1. . .")
        self.engine.data.calculateTreePosteriorProbabilities()
        sum = 0.0
        for t in self.engine.data.trees:
            sum += self.engine.data.trees[t].postprob
        if sum > 1.0:
            print "failed!"
            print "The sum of posterior probability values for trees is " + sum.__str__() + " (it's greater than 1.0)"
            raise AssertionError
        if sum < 0.99:
            print "failed!"
            print "The sum of posterior probability values for trees is " + sum.__str__() + " (it's less than 1.0)"
            raise AssertionError
        print "SUCCESS."
    
    
    def test_calculateMlAncestralSequence(self):
        print "\n--> test_calculateMlAncestralSequence"
        
        sys.stdout.write("\n. test 1. . .")
        self.argParser.setArg("--outgroup","[taxonD]")
        self.argParser.setArg("--ingroup","[taxonA,taxonB,taxonC]")
        self.engine.getAncestralSequence() 
        #self.engine.calculateMlAncestralSequences(["taxonA", "taxonB", "taxonC"], ["taxonD"])
        if os.path.exists(self.argParser.getArg("--outputdir") + "/ancestor.out.txt") == False:
            print "failed!"
            print "I cannot find the ML ancestral sequence at: " + self.argParser.getArg("--outputdir") + "/ancestor.out.txt"
            raise AssertionError
        print "SUCCESS."

        sys.stdout.write("\n. test 2. . .")
        self.engine.calculateMlAncestralSequences(["taxonA", "taxonC"], ["taxonD"])
        if os.path.exists(self.argParser.getArg("--outputdir") + "/ancestor.out.txt") == False:
            print "failed!"
            print "I cannot find the ML ancestral sequence at: " + self.argParser.getArg("--outputdir") + "/ancestor.out.txt"
            raise AssertionError
        print "SUCCESS."
         
    # this test requires that your first call test_calculateTreePosteriorProbabilities()
    def test_calculateMapAncestralSequence(self):
        print "\n--> test_calculateMapAncestralSequence"
        
        sys.stdout.write("\n. test 1. . .")
        self.engine.calculateMapAncestralSequence(["taxonA", "taxonB", "taxonC"], ["taxonD"])

        if os.path.exists(self.argParser.getArg("--outputdir") + "/ancestor-map.dat") == False:
            print "failed!"
            print "I cannot find the MAP ancestral sequence at: " + self.argParser.getArg("--outputdir") + "/ancestor-map.dat"
            raise AssertionError
        print "SUCCESS."
        
        f = open(self.argParser.getArg("--outputdir") + "/ancestor-map.dat", "r")
        mapAnc = Node(-1)
        lines = f.readlines()
        for l in lines:
            # skip any header information, look for data:
            if l == re.match("^\d+\s", l):
                l = l.strip()
                tokens = l.split()
                mapAnc.addSite( int(tokens[0]) )
                for i in range(1, tokens.__len__() - 1):
                    mapAnc.sites[ int(tokens[0]) ].addState( tokens[i], float(tokens[i+1]) ) 
        f.close()

        sys.stdout.write("\n. test 2. . .")
        self.engine.calculateMapAncestralSequence(["taxonA", "taxonC"], ["taxonD"])

        if os.path.exists(self.argParser.getArg("--outputdir") + "/ancestor-map.dat") == False:
            print "failed!"
            print "I cannot find the MAP ancestral sequence at: " + self.argParser.getArg("--outputdir") + "/ancestor-map.dat"
            raise AssertionError
        print "SUCCESS."
        
        f = open(self.argParser.getArg("--outputdir") + "/ancestor-map.dat", "r")
        mapAnc = Node(-1)
        lines = f.readlines()
        for l in lines:
            # skip any header information, look for data:
            if l == re.match("^\d+\s", l):
                l = l.strip()
                tokens = l.split()
                mapAnc.addSite( int(tokens[0]) )
                for i in range(1, tokens.__len__() - 1):
                    mapAnc.sites[ int(tokens[0]) ].addState( tokens[i], float(tokens[i+1]) ) 
        f.close()

        

        sys.stdout.write("\n. test 3. . .")
        for s in mapAnc.sites:
            sum = 0.0
            for c in mapAnc.sites[s].stateProbs:
                sum += mapAnc.sites[s].stateProbs[c]
            if sum > 1.0:
                print "failed!"
                print "The sum of posterior probability values for site " + s.__str__() + " is " + sum.__str__() + " (it's greater than 1.0)"
                raise AssertionError
            if sum < 0.99:
                print "failed!"
                print "The sum of posterior probability values for site " + s.__str__() + " is " + sum.__str__() + " (it's less than 1.0)"
                raise AssertionError
        print "SUCCESS."
        
        
#        sys.stdout.write("\n. test 4. . .")
#        
#        # This test assumes that all three trees in self.engine have post-prob = 0.3333
#        
#        f = open(self.argParser.getArg("--outputdir") + "/tree1/node6.dat", "w")
#        f.write("1  L 0.9  I 0.1\n")
#        f.write("2  V 0.7  A 0.3\n")
#        f.close()
#        
#        f = open(self.argParser.getArg("--outputdir") + "/tree2/node7.dat", "w")
#        f.write("1  L 0.6  I 0.4\n")
#        f.write("2  V 0.1  A 0.9\n")
#        f.close()
#
#        f = open(self.argParser.getArg("--outputdir") + "/tree3/node5.dat", "w")
#        f.write("1  L 0.6  I 0.4\n")
#        f.write("2  V 0.1  A 0.9\n")
#        f.close()
#                        
#        ingroup = ["taxonA","taxonB"]
#        outgroup = ["taxonD"]  
#        mapAnc = self.engine.data.calculateMapAncestralSequence(ingroup, outgroup)
#    
#        if mapAnc.sites[1].stateProbs["L"] != (0.9 / 3.0 ) + (0.6 / 3.0) + (0.6 / 3.0):
#            print "failed!"
#            print "The empirical Bayesian posterior probability for state L at site 1 is: " + float(mapAnc.sites[1].stateProbs["L"]).__str__() + ", but it should be 0.7"
#            raise AssertionError 
#        else:
#            print "SUCCESS."
    
    #
    # main(). . .
    #
    
    def test_engine(self):
        self.outDir = os.getcwd() + "/lazarus_test"
        if os.path.exists(self.outDir) == False:
            os.system("mkdir " + self.outDir)
        sys.argv.append("--outputdir")
        sys.argv.append(self.outDir) 
        self.argParser = ArgParser(sys.argv)
        self.engine = Engine(self.argParser)
        
        self.test_importFiles()
        self.test_startPamlJobs()
        self.test_splitRatesAndRst()
        self.test_saveState()

        self.test_calculateTreePosteriorProbabilities()
        self.test_writeTreeSummary()
        self.test_calculateMlAncestralSequence()
        self.test_calculateMapAncestralSequence()
        
        #os.system("rm -r " + self.outDir)
    