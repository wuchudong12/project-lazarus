##############################################################################
#
# Engine class: 
# a class for initializing and managing the PAML execution, and then
# parsing the results.
#
# Author: Victor Hanson-Smith
# Contact: victorhs@cs.uoregon.edu
#
###############################################################################
from logger import *
from dataWharehouse import *
from argParser import *
from pamlJob import *
from treelib import *
from util import *
import os, pickle, re

class Engine:
    def __init__(self, argParser):
        self.argParser = argParser
        self.logger = Logger(argParser.getOptionalArg("--verbose"), argParser.getArg("--outputdir") )
        self.data = DataWharehouse(argParser)
        self.pamljobs = {} # key = tree number, value = CodemlJob object or Baseml object
    
    def saveState(self):
        outputDirectory = self.argParser.getArg("--outputdir")
        if os.path.exists(outputDirectory) == False:
            os.system("mkdir " + outputDirectory)
        
        fasrL = open(outputDirectory + "/lazarus_workspace.save", "w")
        pickle.dump(self, fasrL)
        fasrL.close()
        
    def floatToString(self, f):
        """Rounds to three decimal places."""
        return "%.3f"%f
            
      
    def importFiles(self, useAminoAcids=True):
        """Import the alignment, tree(s), and model specified by the
        command-line parameters --alignment, --tree, and --model.
        If any of the files are problematic, an error message will be reported
        and the process will exit.  Upon successful import, this function will
        write PAML-friendly files named reformatted_alignment.phy and 
        reformatted_tree.tre."""
        
        alignmentpath = self.argParser.getArg("--alignment")
        treepath = self.argParser.getArg("--tree")
        if useAminoAcids:
            modelpath = self.argParser.getArg("--model")
        else:
            modelName = self.argParser.getArg("--model")
        
        
        al = self.getAlignmentObject()
                
                
        # load the tree(s):
        if os.path.exists(treepath) == False:
            self.logger.throwError("I cannot find your tree file at: " + treepath.__str__() )
            
        f = open(treepath, "r")
        self.logger.printMessage("I found your tree file: " +  treepath, 1)
        lines = f.readlines()
        self.logger.printMessage("...it appears to contain " + lines.__len__().__str__() + " trees.", 1)
        f.close()
        self.logger.updateStatus("Loading your tree(s): Success")        

        
        # Verify that the tree(s) and the alignment contain the same taxa:
        countTrees = 0
        for line in lines:
            if line.__len__() > 3: # only look at non-empty lines...
                try:
                    thisTree = LoadTree(treestring = line.strip())
                    
                    treeNames = []
                    for t in thisTree.tips():
                        treeNames.append(t.Name)
                    
                    countTrees += 1
                    # are all the tips of the tree in our alignment?
                    for tip in treeNames:
                        if al.getSeqNames().__contains__(tip.__str__()) == False:
                            for n in al.getSeqNames():
                                print "'", n, "'"
                            for n in treeNames:
                                print ".", n, "."
                            self.logger.throwError("Your tree contains a taxon named " + tip.__str__() + ", which I cannot find in your alignment.")
                            
                        
                    # the reverse question: are all the alignment taxa in this tree?
                    for taxon in al.getSeqNames():
                        if treeNames.__contains__(taxon) == False:
                            self.logger.throwError("Your alignment contains a taxon named " + taxon + ", which I cannot find in tree # " + countTrees.__str__() + " of your tree file." )  
                except parse.record.FileFormatError:
                    self.logger.throwError("I had some trouble reading your tree file at this line:\n" + line)  

        # check the model:
        if useAminoAcids:
            r = relative_path_to_abs_path(modelpath)
            #print modelpath
            #print r
            if r == False:
                self.logger.throwError("I cannot find your evolutionary model." )
            elif os.path.exists(r) == False:
                self.logger.throwError("I cannot find your model file at: " + r.__str__() )
            else:
                self.logger.printMessage("I found your evolutionary model at: " + r.__str__(), 1)
                modelpath = r
        else:
            known_models = ["JC69", "K80", "F81", "F84", "HKY85", "T92", "TN93", "REV", "UNREST", "REVu", "UNRESTu"]
            if False == known_models.__contains__(modelName):
                self.logger.throwError("BaseML doesn't know about the nucleotide model " + modelName + ".  Known models are: " + known_models.__str__())
            self.logger.updateStatus("Loading your evolutionary model: Success")

        # Before we write any output, verify the output directory
        outputDirectory = self.argParser.getArg("--outputdir")
        if os.path.exists(outputDirectory) == False:
            os.system("mkdir " + outputDirectory)

        # Write the alignment to outputDirectory + "/reformatted_alignment.phy"
        tokens = al.__str__().split("\n")
        cleantokens = []
        for t in tokens:
            if False == t.__contains__(">"):
                t = re.sub("X", "?", t)
                t = re.sub("B", "?", t)
                t = re.sub("J", "?", t)
                t = re.sub("Z", "?", t)
            cleantokens.append(t)
        ag = open(outputDirectory + "/reformatted_alignment.phy", "w")
        ag.write(al.getSeqNames().__len__().__str__() + "  " + al.__len__().__str__() + "\n")
        for t in cleantokens:
            #print "engine 132", t
            if t.__contains__(">"):
                t = re.sub(">", "", t)
                ag.write("\n" + t + "   ")
            else:
                ag.write(t)
        ag.close()
             
        # Write the tree to thisExecutionDir + "/reformatted_tree.tre"   
        countTrees = 0
        for line in lines:
            if line.__len__() > 3: # ignore empty lines:
                countTrees += 1
                
                # create one subdirectory per tree
                thisTreeDir = outputDirectory + "/tree" + countTrees.__str__()
                if os.path.exists( thisTreeDir ) == False:
                    os.system("mkdir " + thisTreeDir )
                
                # build a subdirectory in which PAML will operate
                thisExecutionDir = thisTreeDir + "/pamlWorkspace"
                if os.path.exists( thisExecutionDir ) == False:
                    os.system("mkdir " + thisExecutionDir )
                
                # print this tree into the PAMl workspace
                ts = line.strip()
                ts = re.sub("\)\d+\.*\d+\:", "):", ts)
                ts = re.sub("\[", "", ts)
                ts = re.sub("\]", "", ts)
                thisTree = LoadTree(treestring = ts)
                tg = open(thisExecutionDir + "/reformatted_tree.tre", "w")
                tg.write( thisTree.tips().__len__().__str__() + "  1" + "\n" )
                if False == line.__contains__(";"):
                    line = line + ";"
                tg.write(line + "\n")
                tg.close()
                
                # build a soft link to the alignment (located in the output directory):
                #if os.path.exists(thisExecutionDir + "/reformatted_alignment.phy"):
                #    os.system("rm " + thisExecutionDir + "/reformatted_alignment.phy")
                # for debugging:
                #os.system("ls -alh " + outputDirectory)
                #os.system("ls -alh " + thisExecutionDir)
                #print "\n\n engine 173 \n\n"
                #os.system("pwd")
                #os.system("ls -alh")
                os.system("cp " + outputDirectory + "/reformatted_alignment.phy " + thisExecutionDir + "/reformatted_alignment.phy")
                #os.system("cp " + outputDirectory + "/reformatted_alignment.phy " + thisExecutionDir + "/")#reformatted_alignment.phy")
                #os.system("ln -s " + outputDirectory + "/reformatted_alignment.phy " + thisExecutionDir + "/reformatted_alignment.phy")
                
                # build a soft link to the model, in each execution directory:
                if useAminoAcids:
                    # remove any previous soft links:
                    if os.path.exists(thisExecutionDir + "/model.dat"):
                        os.system("rm " + thisExecutionDir + "/model.dat")
                    #print "\n. making softlink for model here: " + thisExecutionDir + "/model.dat"
                    command = "ln -sf " + modelpath + " " + thisExecutionDir + "/model.dat"
                    #print command
                    os.system(command)
        
                self.data.addTree(countTrees)
                if useAminoAcids:
                    self.buildCodemlJob(countTrees, thisExecutionDir, "reformatted_tree.tre", "reformatted_alignment.phy", "model.dat")
                else:
                    self.buildBasemlJob(countTrees, thisExecutionDir, "reformatted_tree.tre", "reformatted_alignment.phy", modelName)

    #
    # This is a helper method, used by the method self.importFiles
    # This method builds a PAML job, which performs ASR on a single tree
    #
    # INPUT PARAMETERS:
    # treenumber = the tree ID that this PAML job is operating on.  This treenumber corresponds to the tree IDs in self.data.trees
    # executionDir = the PAML workspace for this job.
    #   (we assume here that executionDir is a valid directory)
    # treeFileName, alignFileName, modelFileName = pathless filenames for these necessary files.
    #
    # At the end of this method, self.pamljobs contains a valid CodemlJob object.
    #
    def buildCodemlJob(self, treenumber, executionDir, treeFileName, alignFileName, modelFileName):     
        """This is a helper method, used by the method self.importFiles."""
        self.pamljobs[treenumber] = CodemlJob()
        self.pamljobs[treenumber].treeNumber = treenumber
        self.pamljobs[treenumber].outputDirectory = executionDir
        self.pamljobs[treenumber].treePath = executionDir + "/" + treeFileName
        self.pamljobs[treenumber].alignmentPath = executionDir + "/" + alignFileName
        self.pamljobs[treenumber].modelPath = executionDir + "/" + modelFileName
        self.pamljobs[treenumber].executionDirectory = executionDir
        if self.argParser.doesContainArg("--branch_lengths"):
           x = self.argParser.getArg("--branch_lengths")
           if x == "estimate":
               self.pamljobs[treenumber].estimate_branch_lengths = True
           elif x == "fixed":
               self.pamljobs[treenumber].estimate_branch_lengths = False
        if self.argParser.doesContainArg("--asrv"):
            x = self.argParser.getArg("--asrv")
            if x == "0":
                self.pamljobs[treenumber].among_site_rate_variation = False
                self.fix_asrv = False
            else:
                self.pamljobs[treenumber].among_site_rate_variation = True
                self.pamljobs[treenumber].ncat_gamma = float(x)
        else:
            self.pamljobs[treenumber].among_site_rate_variation = False
            self.fix_asrv = False          
        if self.argParser.doesContainArg("--cleandata"):
            x = self.argParser.getArg("--cleandata")
            self.pamljobs[treenumber].cleandata = int(x)
        if self.argParser.doesContainArg("--fix_asrv"):
            self.pamljobs[treenumber].fix_asrv = True
        if self.argParser.doesContainArg("--alpha"):
            x = self.argParser.getArg("--alpha")
            self.pamljobs[treenumber].alpha = float(x)

    def buildBasemlJob(self, treenumber, executionDir, treeFileName, alignFileName, modelName):     
        """This is a helper method, used by the method self.importFiles."""
        self.pamljobs[treenumber] = BasemlJob()
        self.pamljobs[treenumber].treeNumber = treenumber
        self.pamljobs[treenumber].outputDirectory = executionDir
        self.pamljobs[treenumber].treePath = executionDir + "/" + treeFileName
        self.pamljobs[treenumber].alignmentPath = executionDir + "/" + alignFileName
        self.pamljobs[treenumber].modelName = modelName
        self.pamljobs[treenumber].executionDirectory = executionDir
        if self.argParser.doesContainArg("--branch_lengths"):
           x = self.argParser.getArg("--branch_lengths")
           if x == "estimate":
               self.pamljobs[treenumber].estimate_branch_lengths = True
           elif x == "fixed":
               self.pamljobs[treenumber].estimate_branch_lengths = False
        if self.argParser.doesContainArg("--asrv"):
            x = self.argParser.getArg("--asrv")
            if x == "0":
                self.pamljobs[treenumber].among_site_rate_variation = False
                self.fix_asrv = False
            else:
                self.pamljobs[treenumber].among_site_rate_variation = True
                self.pamljobs[treenumber].ncat_gamma = float(x)
        else:
            self.pamljobs[treenumber].among_site_rate_variation = False
            self.fix_asrv = False 
        if self.argParser.doesContainArg("--cleandata"):
            x = self.argParser.getArg("--cleandata")
            self.pamljobs[treenumber].cleandata = int(x)
        if self.argParser.doesContainArg("--fix_asrv"):
            self.pamljobs[treenumber].fix_asrv = True
        if self.argParser.doesContainArg("--alpha"):
            x = self.argParser.getArg("--alpha")
            self.pamljobs[treenumber].alpha = float(x)
    
    def startPamlJobs(self):
        """Starts all jobs in self.pamljobs"""
        # Here is an obvious opportunity for parallelism!
        # The jobs can be distributed to other nodes, etc.
        for j in self.pamljobs:
            self.logger.updateStatus("Reconstructing ancestral states on tree " + j.__str__() + ". . ." )
            if self.argParser.doesContainArg("--skippaml"):   
                pass
            else:
                self.pamljobs[j].startJob()
            flag = self.pamljobs[j].verifyResults()
            if flag != True:
                self.logger.throwError(flag)

    def cleanupPamlResults(self):
        for job in self.pamljobs:
            self.logger.updateStatus("Cleaning-up the PAML output for tree " + job.__str__() + ". . .")
            path1 = self.pamljobs[job].outputDirectory + "/out.paml"
            path2 = self.pamljobs[job].outputDirectory + "/rst"
            path3 = self.pamljobs[job].outputDirectory + "/rates"
            os.system("rm -f " + path1)
            os.system("rm -f " + path2)
            os.system("rm -f " + path3)

    #
    # Read the files named 'rst' and 'out.paml', which are assumed to exist in the output directory.
    #
    def parsePamlResults(self):
        """For each PamlJob 't' in self.pamljobs, read the files named 'rst' and 'out.paml', which are assumed to exist in self.pamljobs[t].outputDirectory."""
        for t in self.pamljobs:
            self.logger.updateStatus("Parsing the PAML output for tree " + t.__str__() + ". . .")
            self.parsePamlFiles(self.pamljobs[t].outputDirectory + "/out.paml", self.pamljobs[t].outputDirectory + "/rst", t)        
        self.writeTreeSummary()
            

    #
    # This method parses the output files from PAML ('rst' and 'out.paml'), and builds
    # a directory hierarchy containing small text files with information about ancestral
    # states on the the tree numbered 'treenumber'.
    #
    # This method must be called once for each tree.
    # (in our current architecture, we invoke a new PAML session for the analysis of each tree.
    # Therefore, with N trees, we have N rst files and N out.paml files).
    #
    # This method makes many assumptions about the format of the PAML files.
    # If you upgrade PAML, this method might break.
    #
    # 'rstPath' is the filepath to the rst file
    # 'outPath' is the filepath to the out.paml file
    # 'treenumber' is the number of this tree.  After the method is finished, the contents
    #    of the tree and its nodes will be stored in self.engine.data.trees[treenumber]
    #
    def parsePamlFiles(self, outPath, rstPath, treenumber):
        """This is a helper function for the method self.parsePamlResults."""
        self.logger.printMessage("", 1)
        
        if os.path.exists( rstPath ) == False:
            self.logger.throwError("Your rst file does not exist at: " + rstPath)
        
        if os.path.exists( outPath ) == False:
            self.logger.throwError("Your out.paml file does not exist at: " + outPath)
        
        if self.argParser.getArg("--outputdir") == False:
            self.logger.throwError("You did not specify the required command-line argument: --outputdir")    
        else:
            outputDirectory = self.argParser.getArg("--outputdir")
            if os.path.exists( outputDirectory ) == False:
                self.logger.printMessage("mkdir " + outputDirectory, 3)
                os.system("mkdir " + outputDirectory)
        
        # start by parsing out.paml      
        f = open(outPath, "r")
        found_string = False
        found_lnl = False
        while(True):
            try:
                if found_string and found_lnl:
                    break
                line = f.next()
                if None != re.match("^tree length =", line):
                    line = f.next()
                    line = f.next()
                    line = f.next()
                    line = f.next()
                    if False == self.data.trees.__contains__(treenumber):
                        self.data.addTree(treenumber)
                    self.data.trees[treenumber].newickString = line.strip()
                    #print "FOUND NEWICK tree " + treenumber.__str__() + " = " + line.strip()
                    found_string = True
                if None != re.match("^lnL", line): # does this line contain a log likelihood value?
                    x = re.sub("lnL(.*)\:", "", line)
                    tokens = x.split()
                    likelihood = float( tokens[0] )
                    print "Log likelihood of tree " + treenumber.__str__() + " = " + likelihood.__str__()
                    if False == self.data.trees.__contains__(treenumber):
                        self.data.addTree(treenumber)
                    self.data.trees[treenumber].lnL = likelihood
                    found_lnl = True
            except StopIteration:
                break
            except IndexError:
                self.logger.printMessage("An IndexError occurred while parsing the PAML output.  See engine.py near line 329.")
                print "An error occurred in Lazarus.  See engine.py, near line 330."
                print "line = " + line
                exit(1)        
        f.close()
        
        # next parse rst:
        f = open(rstPath, "r")
        currentTreeNumber = treenumber     
        currentNodeNumber = None
        while(True):
            try:
                line = f.next()
                if None != re.match("^TREE", line):
                    tokens = line.split()
                    self.logger.printMessage("Collecting ancestors from tree " + currentTreeNumber.__str__() + ". . .", 2)                
                
                if None != re.match("^Prob of best state at each node", line):
                    # THESE TWO LINES ARE CRITICAL!
                    # Here, we spew out a subdirectory containing text files about the tree and ancestral nodes
                    self.spewTree(currentTreeNumber)
                    self.data.trees[currentTreeNumber].nodes = {}
                    
                    currentTreeNumber = None
                
                if None != re.match("^tree with node labels for Rod Page's TreeView", line) and currentTreeNumber != None: # the internal node numbering scheme is unfortunately embedded in the TreeView string
                    line = f.next()
                    line = line.strip()
                    # strip away the node numbering provided by TreeView labels:
                    # These labels add "XX_" before each taxon name, where XX is TreeView node number.
                    # Here, I look for three specific cases of "XX_", to protect against incorrectly removing
                    # "XX_" from taxa names.
                    # 1. "(XX_"
                    line = re.sub("\(\d+\_", "(", line)
                    # 2. " XX_"
                    line = re.sub("\s\d+\_", " ", line)
                    # 3. ",XX_"
                    line = re.sub("\,\d+\_", ", ", line)
                    #print "\n. 423: Importing cladogram", line
                    self.data.trees[currentTreeNumber].treeviewString = line
                    
                if None != re.match("^Prob distribution at node", line) and currentTreeNumber != None:
                    tokens = line.split()
                    currentNodeNumber = tokens[4].replace(",", "")
                    # we found a node!
                    self.data.trees[currentTreeNumber].addNode(currentNodeNumber)
                    
                if None != re.match("^\s+\d+\s+\d+\s+\S+\:\s", line) and currentTreeNumber != None:
                    tokens = line.split()
                    siteNumber = int(tokens[0])
                    # we found a site!                
                    self.data.trees[currentTreeNumber].nodes[currentNodeNumber].addSite(siteNumber)
                    
                    for i in range(3, tokens.__len__() ):
                        subtokens = tokens[i].split("(")
                        state = subtokens[0]
                        postprob = float(subtokens[1].replace(")", ""))

                        # we found a state!
                        self.data.trees[currentTreeNumber].nodes[currentNodeNumber].sites[siteNumber].addState(state, postprob)
                        
                        #print "tree: " + currentTreeNumber.__str__() + ", node: " + currentNodeNumber.__str__() + ", site: " + siteNumber.__str__() + ", state: " + state.__str__() + ", prob: " + postprob.__str__()
            
            except StopIteration:
                # This line is important, otherwise we don't print the Nth tree.
                f.close()
                break
        f.close()
    
    #
    # This method is a helper function for 'splitRatesAndRst'
    #
    # This creates a tree subdirectory in the outputDirectory.
    # Inside this subdirectory, the tree and the nodes are written
    # into text files.
    #
    # 't' is the ID number of a tree, not the tree object itself.
    #         
    def spewTree(self, t):
        # ensure the output directory exists:
        outputDirectory = self.argParser.getArg("--outputdir")
        if os.path.exists( outputDirectory ) == False:
            os.system("mkdir " + outputDirectory)
        # ensure the tree subdirectory exists:
        treeDirectory = outputDirectory + "/tree" + t.__str__()
        if os.path.exists( treeDirectory ) == False:
            os.system("mkdir " + treeDirectory)
        
        self.logger.printMessage("Building output for tree " + t.__str__() + ". . .", 1)
        
        # write treeX.newick. . .    
        f = open(treeDirectory + "/tree" + t.__str__() + ".txt", "w")
        f.write("TREE #" + t.__str__() + " with branch lengths:\n")
        f.write( self.data.trees[t].newickString + "\n")
        f.write("TREE #" + t.__str__() + " as a cladogram:\n")
        f.write( self.data.trees[t].treeviewString + "\n")
        f.close()
        
        # should we adjust the nodes for indel placement?
        if self.argParser.doesContainArg("--gapcorrect"):
            print "Placing ancestral gaps...."
            self.place_ancestral_gaps(t)
    
        
        # write nodeX.dat (one for each node). . .
        # and write nodeX.motif
        for n in self.data.trees[t].nodes:
            self.data.trees[t].nodes[n].writeSequenceDetail(treeDirectory + "/node" + n.__str__() + ".dat")
            
            #alph = "aa"
            #if self.argParser.doesContainArg("--baseml"):
            #    alph = "nt"
            #self.data.trees[t].nodes[n].writeSequenceMotif(treeDirectory + "/node" + n.__str__() + ".motif", alphabet = alph)
    
    
    # pre assumptions: 
    # 1. --outgroup is given at the command-line
    # 2. the specified outgroup is valid.
    # 3. the tree file contains only a single tree
    def place_ancestral_gaps(self, t):
        """ Parsimoniously places ancestral gaps in the tree indexed by self.data.trees[t], according to the alignment (specified by the --alignment parameter)"""                
        outputDirectory = self.argParser.getArg("--outputdir")
        if os.path.exists( outputDirectory ) == False:
            os.system("mkdir " + outputDirectory)
        # get the outgroup
        outgroup = self.argParser.getArg("--outgroup")        
        if outgroup != False:
            # convert the string y into a list of taxa names:
            if outgroup.__contains__("[") == False:
                self.logger.throwError("Your list of outgroup taxa does not appear to be formatted correctly.\nI think you forgot to include the opening bracket '['.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            if outgroup.__contains__("]") == False:
                self.logger.throwError("Your list of outgroup taxa does not appear to be formatted correctly.\nI think you forgot to include the closing bracket ']'.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            ys = re.sub("[\[\]]", "", outgroup)
            outgroup = ys.split(",")
            self.checkTaxa(outgroup, [t])
        else:
            self.logger.throwError("The argument --gapcorrect requires the argument --outgroup")
        
#         # load the tree as an MIT Tree object (see treelib.py):       
#         f = open(outputDirectory + "/laztree.tre", "w")
#         f.write(self.data.trees[t].treeviewString + "\n")
#         f.close()
#         tree = MITTree()
#         tree.readNewick( outputDirectory + "/laztree.tre")
#         outgroup_nodes = []
#         for o in outgroup:
#             outgroup_nodes.append( tree.nodes[ o ] )
# 
#         print "\n. 533:", self.data.trees[t].treeviewString
#         # re-root the tree according to the outgroup
#         # This introduces a *new* node with an unknown name.
#         outgroup_root = lca(outgroup_nodes)         # using treelib.py
#         r = reroot(tree, outgroup_root.name, onBranch=True)         # using treelib.py    
        
        # This is the problem line:
        #self.data.trees[t].treeviewString = r.getOnelineNewick()
        
        #
        # Tangent, to reroot the tree using dendropy
        # 
        
        
        #
        # December 2014: new way:
        #
#         import dendropy
#         dtree = dendropy.Tree.get_from_string(self.data.trees[t].treeviewString,"newick")
#         dtree.update_splits(delete_outdegree_one=False)
#         
#         """Root the tree, temporarily, at a terminal node."""
#         dtree.reroot_at_midpoint(update_splits=True, delete_outdegree_one=True)
#         
#         """And now re-root at the outgroup mrca"""
#         mrca = dtree.mrca(taxon_labels=outgroup)
#         candidate_edges = []
#         for edge in dtree.postorder_edge_iter():
#             if edge.tail_node == mrca or edge.head_node == mrca:
#                 candidate_edges.append( edge )           
#         dtree.reroot_at_edge(mrca.edge, update_splits=False, delete_outdegree_one=True)


        #
        # old way:
        #
        import dendropy
        #print "546:", self.data.trees[t].treeviewString
        dtree = dendropy.Tree.get_from_string(self.data.trees[t].treeviewString,"newick")
        dtree.update_splits(delete_outdegree_one=True)
        #print "572:", dtree.as_string('newick')
        
        #print "573:", outgroup
        #outgroup = ['DEHA2A04378g', 'PGUG00211.1','CTRG04502.3']
#         not_outgroup = []
#         for n in dtree.nodes():
#             print n.taxon
#             if n.taxon not in outgroup and n.taxon != None:
#                 not_outgroup.append( n.taxon )
#         not_out_mrca = dtree.mrca(taxon_labels = not_outgroup)

        #dtree.reroot_at_node(new_root_node=not_out_mrca)
        
        mrca = dtree.mrca(taxon_labels=outgroup)
        
        """New for March 2015:"""
        if outgroup.__len__() == 1:
            print "590", outgroup
            """Some old versions of Dendropy don't allow rerooting at a leaf,
                so we need to ensure that we're actually dealing with an internal node."""
            for ln in dtree.leaf_nodes():
                print "594 leaf label:", ln.name, ln
                if ln.label == outgroup[0]:
                    mrca = ln.parent_node
                    print "595:", mrca.label
                    break
        
        #dtree.update_splits(delete_outdegree_one=False)
        #dtree.reroot_at_midpoint(update_splits=False, delete_outdegree_one=False)

        #print "579:", dtree.__str__()
        
        dtree.reroot_at_node(mrca, update_splits=True)
        #dtree.reroot_at_edge(mrca.edge, update_splits=True)
        
        
        #print "577", dtree.__str__()
        
        dtree_newick = dtree.as_string('newick')
        #print "551:", dtree_newick
        self.data.trees[t].treeviewString = dtree_newick
        
        #
        # end tangent
        # 
        
        #print "\n. 547:", self.data.trees[t].treeviewString
        # load the rerooted tree as a PyCogent Tree object:
        self.data.trees[t].temporary_tree_object = LoadTree(treestring = self.data.trees[t].treeviewString)
    
        #
        # continue here
        #

        # recursively place gaps
        al = self.getAlignmentObject()
        for s in range(0, al.__len__() ):
            print "\n. Placing gaps for column ", s+1, " of ", al.__len__()
            self.postorder_traversal(t, self.data.trees[t].temporary_tree_object.Name, al[s], s + 1) # s ranges from [0,...], but column numbers range from [1,...]
            self.preorder_traversal(t, self.data.trees[t].temporary_tree_object.Name, al[s], s + 1)
        
        # clear our temporary variables
        self.data.trees[t].temporary_tree_object = None
        self.data.trees[t].temporary_states = {} 
        # and remove the temporary tree file:
        os.system("rm " + outputDirectory + "/laztree.tre")
            
    # this recursion is horribly inefficient, but it works for now.
    #
    # tree_number = the index in self.data.trees
    # mit_tree = an initialized treelib object
    # node_name = a node name or node number
    # column = such as al[i], where al = LoadSeqs(...), and i in range [0, al.__len__()]
    # column_number = a valid index in the self.data.trees[...].nodes[...].sites[ column_number ]
    #
    # Post Condition:
    # self.data.trees[t] will be modified, such that gaps are placed at some sites.
    def postorder_traversal(self, tree_number, node_name, column, column_number):
        """a helper function for places_ancestral_gaps"""

        n = None

        # 1. is this node the root of the entire tree? 
        # this case is important because the re-rooting process likely added a node whose
        # name is unknown.
        if False == self.data.trees[tree_number].temporary_tree_object.getNodeNames().__contains__( node_name ):
            n = self.data.trees[tree_number].temporary_tree_object
        # 2. is this node a leaf?
        elif self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName(node_name).isTip():
            if column.getGappedSeq(node_name).__str__() == "-":
                return 0 # no state at this site
            else:
                return 1 # yes, there exists a state at this site
        # 3. we're dealing with an internal node...
        else:
            if n == None:
                n = self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName( node_name )
            union = None
            for child in n.Children:    
                x = self.postorder_traversal(tree_number, child.Name, column, column_number)
                #print child.Name, x
                self.data.trees[tree_number].temporary_states[child.Name] = x
                if union == None:
                    union = x
                elif union == 0.5: # the union is ambigious, assign the dominant state
                    union = x
                elif x != union and x != 0.5: # some of my children have indels, some do not:
                    return 0.5
        return union # all my children have the same indel status, return this status up to my parent
    
    def preorder_traversal(self, tree_number, node_name, column, column_number):
        """a helper function for places_ancestral_gaps"""
        
        #print "602:", self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName(node_name)
        #print "603:", self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName(node_name).isTip()
        
        #print "preorder_traversal, node_name=", node_name, " column_number=", column_number
        n = None
        # 1. is this node the root of the entire tree? 
        # this case is important because the re-rooting process likely added a node whose
        # name is unknown.
        if False == self.data.trees[tree_number].temporary_tree_object.getNodeNames().__contains__( node_name ):
            n = self.data.trees[tree_number].temporary_tree_object
        # 2. is this node a leaf?
        elif self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName(node_name).isTip():
            return
        # 3. else, we're dealing with an internal node...
        else:
            if n == None:
                n = self.data.trees[tree_number].temporary_tree_object.getNodeMatchingName( node_name )
            if n.Parent != None and self.data.trees[tree_number].temporary_states.__contains__(n.Parent.Name):
                if self.data.trees[tree_number].temporary_states[n.Parent.Name] == self.data.trees[tree_number].temporary_states[n.Name]:
                    if self.data.trees[tree_number].temporary_states[node_name] == 0:
                        #print "adding gap to node:", node_name, ", column:", column_number
                        self.data.trees[tree_number].nodes[node_name].sites[column_number].stateProbs["-"] = 100.0
                    if self.data.trees[tree_number].temporary_states[node_name] == 0.5:
                        mlstate = self.data.trees[tree_number].nodes[node_name.__str__()].sites[column_number].find_ML_state()
                        #print "mlstatepp = self.data.trees[",tree_number,"].nodes[",node_name.__str__(),"].sites[",column_number,"].stateProbs[",mlstate,"]"
                        mlstatepp = self.data.trees[tree_number].nodes[node_name.__str__()].sites[column_number].stateProbs[mlstate]
                        #print "623:", node_name
                        self.data.trees[tree_number].nodes[node_name.__str__()].sites[column_number].stateProbs[mlstate.swapcase()] = mlstatepp
                elif self.data.trees[tree_number].temporary_states[n.Parent.Name] == 0.5 and self.data.trees[tree_number].temporary_states[n.Name] == 0:
                    #print "626:", node_name
                    self.data.trees[tree_number].nodes[node_name.__str__()].sites[column_number].stateProbs["-"] = 100.0
                elif self.data.trees[tree_number].temporary_states[n.Parent.Name] == 0 and self.data.trees[tree_number].temporary_states[n.Name] == 0.5:
                    #print "629:", node_name, column_number
                    self.data.trees[tree_number].nodes[node_name.__str__()].sites[column_number].stateProbs["-"] = 100.0
                    self.data.trees[tree_number].temporary_states[n.Name] == 0
        
        for child in n.Children:
            self.preorder_traversal(tree_number, child.Name, column, column_number)
    
    
    def writeTreeSummary(self):
        """Write a file named tree_summary.txt into the output directory.  This file
        contains a list of trees in your analysis, their repsective log likelihood, and
        their Bayesian posterior probability."""
        outputDirectory = self.argParser.getArg("--outputdir")
        if os.path.exists( outputDirectory ) == False:
            os.system("mkdir " + outputDirectory)
            
        self.data.calculateTreePosteriorProbabilities()
                
        treeLikelihoods = {} # key = likelihood, value = array of trees with that value 
        for t in self.data.trees:
            l = self.data.trees[t].lnL
            if treeLikelihoods.__contains__( l ):
                treeLikelihoods[l].append(t)
            else:
                treeLikelihoods[l] = []
                treeLikelihoods[l].append(t)
        f = open(outputDirectory + "/tree_summary.txt", "w")
        f.write(" tree |   lnL  | Bayesian post. prob\n")            
        z = treeLikelihoods.keys()
        z.sort()
        z.reverse()
        for lvalue in z:
            for tree in treeLikelihoods[lvalue]:
                f.write(tree.__str__() + "  " + self.floatToString(lvalue) + "  " + self.floatToString(self.data.trees[ tree ].postprob) + "\n")                     
        f.close()       
    
    def getAlignmentObject(self):
        """Fetch an alignment object by parsing the --alignment parameter and then calling cogent.LoadSeqs"""
        alignmentpath = self.argParser.getArg("--alignment")
        if False == os.path.exists(alignmentpath):
            self.logger.throwError("I cannot find your alignment file at " + alignmentpath.__str__() )
        try:
            al = LoadSeqs(alignmentpath, aligned = True)
            self.logger.printMessage("I found your alignment file: " +  alignmentpath, 1)
            self.logger.printMessage("...it appears to contain " + al.getNumSeqs().__str__() + " taxon.", 1)
            print al.getSeqNames()
            self.logger.updateStatus("Loading your alignment: Success")
            return al
        except parse.record.FileFormatError:
            self.logger.throwError("Your alignment file " + alignmentpath + " does not appear to be a supported format.")
            exit(1)
        except core.alignment.DataError:
            self.logger.throwError("The sequences in your alignment (" + alignmentpath + ") are not all the same length!\nThis probably occurred because your sequences are not correctly aligned.\n")
            exit(1)

    def getAncestralSequence(self):
        """This method retrieves the ML and MAP ancestral sequences for the ingroup/outgroup definition in the command-line arguments.
        Returns a two-item array: [mlsequence, mapsequence]"""
        
        ingroup = None 
        outgroup = []

        #
        # --outgroup is required.
        #
        y = self.argParser.getArg("--outgroup")
        if y != False:
            # convert the string y into a list of taxa names:
            if y.__contains__("[") == False:
                self.logger.throwError("Your list of outgroup taxa does not appear to be formatted correctly.\nI think you forgot to include the opening bracket '['.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            if y.__contains__("]") == False:
                self.logger.throwError("Your list of outgroup taxa does not appear to be formatted correctly.\nI think you forgot to include the closing bracket ']'.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            ys = re.sub("[\[\]]", "", y)
            outgroup = ys.split(",")
            self.checkTaxa(outgroup, self.data.trees.keys())
        
        #
        # Two ways to specify the desired ancestor: --ingroup, or --ancnode
        #
        y = self.argParser.getOptionalArg("--ingroup")
        if y != False:
            # convert the string y into a list of taxa names:
            if y.__contains__("[") == False:
                self.logger.throwError("Your list of ingroup taxa does not appear to be formatted correctly.\nI think you forgot to include the opening bracket '['.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            if y.__contains__("]") == False:
                self.logger.throwError("Your list of ingroup taxa does not appear to be formatted correctly.\nI think you forgot to include the closing bracket ']'.\nMake sure it looks like this: [taxaName,taxaName,taxaName]\n")
            ys = re.sub("[\[\]]", "", y)
            ingroup = ys.split(",")
            self.checkTaxa(ingroup, self.data.trees.keys())
        y = self.argParser.getOptionalArg("--ancnode")
        if y != False:
            nt = self.data.getMlTreeNumber()        
            ingroup = self.data.trees[nt].getDescendantsOfNode(y)
            
        if ingroup == None or outgroup == None:
            self.logger.throwError("You need to specify the '--ancnode' argument, or the '--ingroup' and '--outgroup' arguments")
        
        trees_sequences = self.calculateMlAncestralSequences(ingroup, outgroup)
        map_sequence = self.calculateMapAncestralSequence(ingroup, outgroup)

        ml_node_name = self.data.trees[ self.data.getMlTreeNumber() ].getLCA(ingroup, outgroup).Name
        
        print "\n. Oh, you want node #", ml_node_name.__str__()
        
        if os.path.exists(self.argParser.getArg("--outputdir") + "/ancestor.out.txt"):
            os.system("rm " + self.argParser.getArg("--outputdir") + "/ancestor.out.txt")
        fout = open(self.argParser.getArg("--outputdir") + "/ancestor.out.txt", "w")
        fout.write("===========================================================================\n")
        fout.write("This file was auto-generated by Lazarus on " + time.localtime().__str__()  + "\n")
        fout.write("\n")
        fout.write("Below are the most-recent-common-shared ancestral sequences\nfor the ingroup on each tree.\n")
        fout.write("Ingroup: \n" + ingroup.__str__() + "\n")
        fout.write("The outgroup is used to root the tree.\n")
        fout.write("Outgroup:\n" + outgroup.__str__() + "\n")
        fout.write("On the ML tree (#" + self.data.getMlTreeNumber().__str__() + "), this is node #" + ml_node_name.__str__() + "\n")
        fout.write("===========================================================================\n")
        for t in trees_sequences.keys():
            fout.write("from TREE " + t.__str__() + ", lnL = " + self.data.trees[t].lnL.__str__() + "\n")
            fout.write(trees_sequences[t] + "\n")
            fout.write("\n")
        if trees_sequences.keys().__len__() > 1:
            fout.write("Averaged Over Trees (weighing each tree by its posterior probability):\n")
            fout.write(map_sequence + "\n")
            fout.close()
        
        print "The ML ancestral sequence for every tree was written to " + self.argParser.getArg("--outputdir") + "/ancestor.out.txt"
        print "The ML posterior ancestral state vector was written to " + self.argParser.getArg("--outputdir") + "/ancestor-ml.dat"
        print "The MAP (over trees) posterior ancestral state vector was written to" + self.argParser.getArg("--outputdir") + "/ancestor-map.dat"
        
    # for each taxa named 't' in the array group: does 't' exist?    
    def checkTaxa(self, group, use_these_trees_numbers):
        for t in group:
            for tr in use_these_trees_numbers:
                if self.data.trees[tr].doesExtantTaxaExist(t) == False:
                    self.logger.throwError("I cannot find the taxon named " + t + " in your trees.\nAre you sure you spelled it correctly?\n")
        
    #
    # This is a helper method for 'getAncestralSequence'  
    # Upon completion, the ML ancestor will be a soft-link the output directory
    #
    # This method returns an array of strings, corresponding the ML ancestral sequence on each tree.
    #
    def calculateMlAncestralSequences(self, ingroup, outgroup):        
        mltree_number = self.data.getMlTreeNumber()
        
        trees_sequences = {}
        for t in self.data.trees.keys():
            n = self.data.trees[t].getLCA(ingroup, outgroup)
            
            if n == False:
                self.logger.throwError("The ancestral node, given your ingroup and outgroup, is ambiguous on tree " + t.__str__() + "\nYour ingroup=" + ingroup.__str__() + ".\nYour outgroup="+outgroup.__str__())
                trees_sequences[t] = "The ancestral node is ambiguous on tree " + t.__str__()
                continue
            
            node = n.Name
            
            if False == os.path.exists(self.argParser.getArg("--outputdir") + "/tree" + t.__str__() + "/node" + node.__str__() + ".dat"):
                self.logger.throwError("The file for ancestral node " + node.__str__() + " on tree " + t.__str__() + " does not exist!  I'm looking for " + self.argParser.getArg("--outputdir") + "/tree" + t.__str__() + "/node" + node.__str__() + ".dat")
            
            if t == mltree_number:
                os.system("cp " + self.argParser.getArg("--outputdir") + "/tree" + t.__str__() + "/node" + node.__str__() + ".dat " + self.argParser.getArg("--outputdir") + "/ancestor-ml.dat")
                
            fin = open(self.argParser.getArg("--outputdir") + "/tree" + t.__str__() + "/node" + node.__str__() + ".dat", "r")    
            sequence = ""
            for l in fin.readlines():
                tokens = l.split()
                sequence += tokens[1]
            fin.close()
            trees_sequences[t] = sequence
        return trees_sequences
       
    def calculateMapAncestralSequence(self, ingroup, outgroup):
        self.data.calculateTreePosteriorProbabilities()
        mapAnc = self.data.calculateMapAncestralSequence(ingroup, outgroup)
        sequence = ""
        fout = open( self.argParser.getArg("--outputdir") + "/ancestor-map.dat", "w" )
        for s in mapAnc.sites.keys():
            x = mapAnc.sites[s].toString()
            fout.write(x + "\n")
            tokens = x.split()
            sequence += tokens[1]
        fout.close()
        return sequence