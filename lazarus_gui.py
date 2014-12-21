#!/usr/bin/python

#########################################################
#
# Lazarus terminal interface
#
# Victor Hanson-Smith
# victorhs@cs.uoregon.edu
#
# If you run this script, you will see a PhyML-like terminal interface.
#
##########################################################
import os, re, sys
import re
import sys
from includeMe import *
from util import *
from version import *

orderedOptions = ["D", "A", "T", "P", "M", "B", "V", "Z", "O", "I", "U"]

optionNames = {"D":"datatype                ",
               "O":"output to this location ",
               "A":"sequence alignment file ",
               "T":"newick tree file        ",
               "B":"branch lengths          ",
               "V":"evo. rate heterogeneity ",
               "Z":"alpha (for +G model)    ",
               "P":"path to folder w/ models",
               "M":"evolutionary model      ",
               "I":"place ancestral gaps    ",
               "U":"outgroup taxa           "
               }

# some options have presets.  For example, we have several built-in amino acid models
aamodels = ["JTT", "WAG", "MtMam", "MtREV", "MtArt", "Dayhoff", "cpREV"]
nucmodels = ["JC69", "K80", "F81", "HKY85", "T92", "TN93", "REV", "UNREST", "REVu", "UNRESTu"]

presetOptions = {"M": aamodels,
                 "D":["nucleotide", "amino acid"],
                 "B":["find ML estimate", "use fixed values in tree file"],
                 "V":["use fixed rate across sites", "use 4 gamma-distributed categories"],
                 "I":["Yes, parsimoniously place gaps.", "No, I'll manually do it later."]}
                 
optionValues = {"D":None,
                "O":os.getcwd(),
                "A":None,
                "T":None,
                "B":"find ML estimate",
                "V":"use 4 gamma-distributed categories",
                "Z":"find ML estimate",
                "M":None,
                "I":"No, I'll manually do it later.",
                "U":None,
                "P":"current directory",
                "L":"$LAZARUS"
                }

optionErrorMessages = {"D":[None],
                       "O":[None],
                      "A":[None],
                      "T":[None],
                      "B":[None],
                      "M":[None],
                      "P":[None],
                      "V":[None],
                      "Z":[None],
                      "G":[None],
                      "I":[None],
                      "U":[None]}
                        

optionTipsMessages = {"D":[None],
                      "O":[None],
                      "A":[None],
                      "T":[None],
                      "B":[None],
                      "M":[None],
                      "P":[None],
                      "V":[None],
                      "Z":[None],
                      "I":[None],
                      "U":[None]}

# returns the next option to display
def cycleToNextMethod(key):
    if optionValues[key] == None:
        currentOptionIndex = 0
    else:
        count = 0
        for i in presetOptions[ key ]:
            if i == optionValues[key]:
                currentOptionIndex = count
            count += 1
        if currentOptionIndex == presetOptions[key].__len__() - 1:
            currentOptionIndex = 0
        else:
            currentOptionIndex += 1
    optionValues[key] = presetOptions[key][currentOptionIndex]
    return optionValues[key]

##############################
# verification methods
###############################

def checkFilePath(key):
    if optionValues[key] == None:
        return False
    if os.path.exists( optionValues[key] ) == False:
        optionErrorMessages[key][0] = "The file " + optionValues[key] + " does not appear to exist!  Perhaps it is located in a different folder?  Did you check your spelling?"
        return False
    else:
        optionErrorMessages[key] = [None]
        return True

def checkDirectoryPath(key):
    if optionValues[key] == None:
        return False
    if os.path.exists( optionValues[key] ) == False:        
        optionErrorMessages[key][0] = "The directory " + optionValues[key] + " does not appear to exist!"
        return False
    else:
        optionErrorMessages[key] = [None]
        return True

def checkModel(key):
    if optionValues["P"] == None:
        optionErrorMessages["P"][0] = "You need to specify a folder containing models."
        return False
    if False == os.path.exists(optionValues[key]):
        optionErrorMessages["P"][0] = "I can't find your folder containing models at " + optionValues[key]
        return False
    
    # turn relative path into absolute path
    # this transformation must also happen in the CLI/batch version of Lazarus.
    r = relative_path_to_abs_path(optionValues[key])
    if r == False:
        optionErrorMessages["P"][0] = "Hmm, I can't find the folder containing models at " + optionValues[key] + ".  Try adding ./ before the path, or provide an absolute path from root /"
        return False        

    files = os.listdir(optionValues[key])
    count = 0
    for f in files:
        if f.__contains__(".dat") or f.__contains__(".rates"):
            count += 1
    optionTipsMessages["P"][0] = "O.K. I found " + count.__str__() + " models named with *.dat or *.rates"
    optionErrorMessages["P"] = [None]
    return True
    

def checkAlignment():
    if optionValues["A"] == None:
        return False
    if checkFilePath("A") == False:
        return False
    try:
        al = LoadSeqs(optionValues["A"], aligned = True, moltype = PROTEIN)
        optionTipsMessages["A"][0] = "O.K. I found " + al.getNumSeqs().__str__() + " taxa."
        optionErrorMessages["A"] = [None]
        return True
    except parse.record.FileFormatError:
        optionErrorMessages["A"][0] = "Your alignment file does not appear to be a supported format."
        return False
    except core.alignment.DataError:
        optionErrorMessages["A"][0] = "Your sequences are not all the same length!"
        return False

def checkTrees():
    if optionValues["T"] == None:
        return False
    if checkFilePath("T") == False:
        return False
    fin = open(optionValues["T"], "r")
    lines = fin.readlines()
    if lines.__len__() < 1:
        optionErrorMessages["T"][0] = "Your tree file does not appear to contain any trees!"
        fin.close()
        return False
    else:  
        countTrees = 0
        for line in lines:
            try:
                thisTree = LoadTree(treestring = line.strip())
            except parse.record.FileFormatError:
                optionErrorMessages["T"][0] = "I had some trouble reading your tree file at this line:\n" + line
                fin.close()
                return False
        optionTipsMessages["T"][0] = "O.K. I found " + lines.__len__().__str__() + " tree(s)."
        optionErrorMessages["T"] = [None]
        fin.close()
        return True

def parseOutgroup():
    if optionValues["U"] == None:
        return False
    s = optionValues["U"].strip()
    s = re.sub(" ", "", s)
    optionValues["U"] = s
        

# assumes you've already called crossCheckAlignmentTrees
def crossCheckOutgroupAlignment():
    if optionValues["I"] == "No, I'll manually do it later.":
        return True
    else:
        al = LoadSeqs(optionValues["A"], aligned = True, moltype = PROTEIN)
        outgroup = optionValues["U"]
        if outgroup == None:
            optionTipsMessages["I"].append("To place ancestral gaps, you need to specify one or more outgroup taxa (option 'U').")
            return False
        tokens = outgroup.split(",")
        print tokens
        return
        bad_outgroup = []
        for taxon in tokens:
            if False == al.getSeqNames().__contains__( taxon ):
                bad_outgroup.append( taxon )
        if bad_outgroup.__len__() > 0:
            optionErrorMessages["G"][0] = "Your outgroup contains taxon ", bad_outgroup, " which I cannot find in your alignment." 
            return False
        else:
            return True

def crossCheckAlignmentTrees():
    if checkAlignment() and checkTrees():
        al = LoadSeqs(optionValues["A"], aligned = True, moltype = PROTEIN)
        f = open(optionValues["T"], "r")
        
        lines = f.readlines()   
        countTrees = 0
        for line in lines:
            thisTree = LoadTree(treestring = line.strip())
            
            treeNames = []
            for t in thisTree.tips():
                treeNames.append(t.Name)
            
            countTrees += 1
            # are all the tips of the tree in our alignment?
            for tip in treeNames:
                if al.getSeqNames().__contains__(tip.__str__()) == False:
                    optionErrorMessages["G"].append("Your tree contains taxon " + tip.__str__() + ", which I cannot find in your alignment.")
                    f.close()
                    return False
                    
            # the reverse question: are all the alignment taxa in this tree?
            for taxon in al.getSeqNames():
                if treeNames.__contains__(taxon) == False:
                    optionErrorMessages["G"].append("Your alignment contains taxon " + taxon + ", which I cannot find in tree # " + countTrees.__str__())  
                    f.close()
                    return False

# have all options been completed?
def checkOptionsComplete():
    for i in optionValues:
        if optionValues[i] == None:
            optionErrorMessages["G"].append("You need to specify an option for (" + i + ")")
            return False
    return True 

####################
# paint methods
#####################

def paintHeader():
    os.system("clear")
    print "=============================================="
    print "  LAZARUS for PAML " + get_version_string()
    print "  questions? contact Victor Hanson-Smith" 
    print "          at victorhansonsmith@gmail.com"
    print "=============================================="

def paintOption(key):
    value = optionValues[key]
    if value == None:
        value = "PLEASE SPECIFY"
    print "  " + key + "    " + optionNames[key] + " => " + value
    for i in optionTipsMessages[key]:
        if optionTipsMessages[key][0] != None:
            print "            => " + optionTipsMessages[key][0].__str__()

def paintErrorMessages():
    for i in optionErrorMessages:
        for j in optionErrorMessages[i]:
            if j != None:
                print "ERROR: " + j

def paint():
    paintHeader()
    print "\nSettings for this run:\n"
    for o in orderedOptions:
        paintOption(o)
    print "\n"
    paintErrorMessages()
    print "\n"

########################################
def runLazarus():
    if optionValues["D"] == "nucleotide":
        command = "python " + optionValues["L"] + "/lazarus.py --baseml"
    else:
        command = "python " + optionValues["L"] + "/lazarus.py --codeml"
    command += " --outputdir " + optionValues["O"] 
    command += " --alignment " + optionValues["A"]
    command += " --tree " + optionValues["T"]
    
    mdir = "./"
    if optionValues["P"] != "current directory":
        mdir = optionValues["P"]
    
    if optionValues["D"] == "nucleotide":
        command += " --model " + optionValues["M"]
    else:
        command += " --model " + mdir + "/" + aa_modelPaths[ optionValues["M"] ]
    
    if optionValues["B"] == "find ML estimate":
        command += " --branch_lengths estimate"
    elif optionValues["B"] == "use fixed values in tree file":
        command += " --branch_lengths fixed"
    
    if optionValues["V"] == "use 4 gamma-distributed categories":
        command += " --asrv 4"
        if optionValues["Z"] != "find ML estimate":
            command += " --fix_asrv True"
            command += " --alpha " + optionValues["Z"]
        else:
            command += " --alpha 1.0"
    elif optionValues["V"] == "use fixed rate across sites":
        command += " --asrv 0"
    
    if optionValues["I"] == "Yes, parsimoniously place gaps.":
        command += " --gapcorrect"
    if optionValues["U"] != None:
        command += " --outgroup [" + optionValues["U"] + "]"
    
    print "Invoking Lazarus with this command:"
    print command
    os.system(command)
    exit(1)
    #fout = open("scratch.sh", "w")
    #fout.write("#!/bin/bash\n")
    #fout.write(command)
    #fout.close()
    #os.system("export lazarus='python /Users/victor/Documents/EclipseWork/Lazarus/lazarus.py'")
    #os.system("source ./scratch.sh")
    #os.system("bash scratch.sh")
    #os.system("rm scratch.sh")

##########################################
# Main
##########################################
choice = None
while choice != "X":
    paint()
    
    optionErrorMessages["G"] = []
    choice = raw_input('Are these settings correct? (type Y or letter to change option)')
    choice = choice.upper()
    if choice == "M":
        cycleToNextMethod("M")
    elif choice == "A":
        p = raw_input('Which alignment file? (Fasta or Phylip formats) >')
        optionValues[choice] = p
        checkAlignment()
    elif choice == "T":
        p = raw_input('Which Newick tree file? >')
        optionValues[choice] = p
        checkTrees()
    elif choice == "O":
        p = raw_input('Where is the output directory? >')
        optionValues[choice] = p
        checkDirectoryPath(choice)
    elif choice == "P":
        p = raw_input('Where is the directory containing models? >')
        optionValues[choice] = p
        checkModel(choice)
    elif choice == "L":
        p = raw_input('Where is the directory containing lazarus.py? >')
        optionValues[choice] = p
        if False == os.path.exists(p + "/lazarus.py"):
            optionErrorMessages["L"][0] = "I can't find lazarus.py in the directory " + p
    elif choice == "D":
        cycleToNextMethod("D")
        if optionValues["D"] == "nucleotide":
            presetOptions["M"] = nucmodels
            if aamodels.__contains__( optionValues["M"] ):
                optionValues["M"] = None
        elif optionValues["D"] == "amino acid":
            presetOptions["M"] = aamodels
            if nucmodels.__contains__( optionValues["M"] ):
                optionValues["M"] = None
    elif choice == "I":
        cycleToNextMethod("I")
        #if optionValues["I"] == "Yes, parsimoniously place gaps.":
        #    if optionValues["U"] == None:
        #        optionErrorMessages["G"].append("To place ancestral gaps, you need to specify outgroup taxa (option 'U').")
    elif choice == "U":
        p = raw_input('Which taxa are the outgroup? (specify as comma-separated list) >')
        optionValues[choice] = p
        parseOutgroup()
    elif choice == "B":
        cycleToNextMethod("B")
    elif choice == "V":
        cycleToNextMethod("V")
    elif choice == "Z":
        p = raw_input('What is the alpha value for the gamma-distributed model? >')
        optionValues[choice] = p
        parseOutgroup()
    elif choice == "Y":
        if checkOptionsComplete():
            runLazarus()
    crossCheckAlignmentTrees()
    crossCheckOutgroupAlignment()
    
    workspace.save
    