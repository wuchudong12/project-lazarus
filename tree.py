from node import Node
from cogent import *
import os, re

class Tree:
    def __init__(self, id):
        self.treeNumber = id
        self.newickString = None
        self.treeviewString = None
        self.lnL = None
        self.nodes = {} # key = node number, value = a Node object
                        # the nodes in self.nodes are not hierarchically arranged.
                        # use the self.treeviewString to get at topological information.
        self.exclude = False # in some cases, we'll want to exclude certain trees
        
        self.temporary_tree_object = None # in some cases -- i.e. ancestral gap placement -- we'll
                                    # want to store a temporary PyCogent tree object instead
                                    # of re-instantiating the object over several iterations.
        self.temporary_states = {} # key = node name, value: 0 indel, 1 no indel, 0.5 maybe an indel (for parsimonious indel placement)
                                                                 
    #
    #
    #
    def addNode(self, nodeNumber):
        self.nodes[nodeNumber] = Node(nodeNumber)
            
    #
    # This method returns a list of strings, representing the descendant node
    # names of the ancestral 'nodeNumber'
    #
    def getDescendantsOfNode(self, nodeNumber):
        self.root = LoadTree(treestring = self.treeviewString) # this invokes a method from PyCogent
        anc = self.root.getNodeMatchingName( nodeNumber )
        dn = anc.tips()
        
        returnArray = []
        
        for n in dn:   
            returnArray.append(n.Name)
            
        return returnArray

    #
    # This is a helper method for 'toString'
    #        
    def recurRename(self, namedRoot, unnamedRoot):
        if namedRoot.isTip():
            return [namedRoot, unnamedRoot]
        
        unnamedRoot.Name = namedRoot.Name
        
        #print "namedRoot's children: " + namedRoot.Children.__str__()
        #print "unnnamedRoot's children: " + unnamedRoot.Children.__str__()
        
        nChildren = namedRoot.Children
        
        for i in range(0, nChildren.__len__() ):
            #print "i = " + i.__str__()
            #print "this named Child = " + namedRoot.Children[i].Name
            #print "this unnamed Child = " + unnamedRoot.Children[i].Name
            [namedRoot.Children[i], unnamedRoot.Children[i] ] = self.recurRename(namedRoot.Children[i], unnamedRoot.Children[i])
        
        return [namedRoot, unnamedRoot]
        
    
    #
    # This method returns a string containing a Newick-formatted version of this tree.
    # If self.treeviewString != None, the returned string will contain internal node labels.
    #
    def toString(self):
        r = LoadTree(treestring = self.treeviewString)
        #t = LoadTree(treestring = self.newickString)
        
        # the 'r' tree has XX_name, whereas the 't' tree has just name
        # for the tip nodes:
        #for c in r.tips():
        #    c.Name = re.sub("^\d+\s+", "", c.Name)
        
        #[r,t] = self.recurRename(r, t)
        #print t.getNewick()
        return r.getNewick()

    #
    #
    # 10/9/2008:
    # NOTE! This method is now depracated.  Please use getLCA instead!
    #
    # Returns the node number of the last-shared ancestor of 'descendants',
    # where 'ingroup' is a list of ingroup taxon IDs.
    # 'outgroup' is a list of outgroup taxon IDs.
    #
    # I assume that ingroup and outgroup are mutually exclusive.
    def getAncestralNode(self, ingroup, outgroup):
        self.root = LoadTree(treestring = self.treeviewString) # this invokes a method from PyCogent

        #print "DEBUG: getAncestralNode: " + self.root.__str__()

        ingroupNodes = {}
        for i in ingroup:
            n = self.root.getNodeMatchingName(i)
            n.params["mark"] = 1
            ingroupNodes[i] = n

        outgroupNodes = {}
        for o in outgroup:
            n = self.root.getNodeMatchingName(o)
            outgroupNodes[o] = n
        
        potentialAncestralNodes = [] # a list of PhyloNodes which could potentially be the root.

        #
        # walk up the tree, marking ingroup nodes
        #
        for n in ingroupNodes:
            next = ingroupNodes[n].Parent
            
            while next != None:
                if next.params.__contains__("ingroup"):
                    if next.params["ingroup"] == 1:
                        break
                else:
                    next.params["ingroup"] = 1
                    #print "Ingroup Marking node " + next.Name
                next = next.Parent
                            
        #
        # Clean-up the markings from the previous uphill walk. . .
        # This will seem confusing, but the uphill traversal for the first node will leave a mark
        # on all internal nodes between the leaf and the root.  We want to mark only nodes whose
        # descendants are ALL ingroup taxon. 
        #

        m = self.unmarkRecur(self.root)
        #print "m = " + m.Name
        if m.isroot() == False:
            potentialAncestralNodes.append( m )
            #print "DEBUG: potential ancestral node: " + m.Name
        
        #
        # walk up the tree, marking outgroup nodes
        #
        for n in outgroupNodes:
            next = outgroupNodes[n].Parent
            #print "examing node " + next.Name
            
            while next != None:
                if next.params.__contains__("outgroup"):
                    if next.params["outgroup"] == 1:
                        potentialAncestralNodes.append(next)
                        #print "DEBUG: potential ancestral node: " + next.Name
                        break
                next = next.Parent
                
        if potentialAncestralNodes.__len__() != 1:
            return False
        
        #print "ancestral node = " + potentialAncestralNodes[0].Name
        return potentialAncestralNodes[0].Name
        
    #
    # This is a helper method for 'getAncestralNode'
    # if 'node' contains only one child with the "mark" param, then node has its mark removed
    # and we recur down the tree
    #
    # returns the first found PhyloNode object which contains two more more marked children
    #
    def unmarkRecur(self, node, mark):
        markedChildren = []
        for c in node.Children:
            if c.params.__contains__(mark):
                if c.params[mark] == 1:
                    markedChildren.append(c)           
        if markedChildren.__len__() > 1 or markedChildren.__len__() == 0:
            # we're the bottom!
            return node
        else:
            c = markedChildren[0]
            node.params.__delitem__(mark)
            #print "tree.py 179: unmarking node " + node.Name
            return self.unmarkRecur( markedChildren[0], mark )    
    #
    # Recur down the tree.
    # Find the root of the ingroup, based on the outgroup spec.
    #
    def findRootDeep(self, node):
        if node.istip():
            return False      
        for c in node.Children:
            # leaving the ingroup:
            if True == node.params.__contains__("ingroup"):
                if False == c.params.__contains__("ingroup"):
                    if True == self.containsMark(c, "outgroup"):
                        return node
            # entering the ingroup:
            elif False == node.params.__contains__("ingroup"):
                if True == c.params.__contains__("ingroup"):
                    if False == self.containsMark(c, "outgroup"):
                        return c
            
            x = self.findRootDeep(c)
            if x != False:
                return x
        # if we didn't find an ingroup root:
        return False
    
    #
    # Finds the node(s) that are ancestors of the ingroup, connected by one edge
    # to an ancestor of the outgroup
    #
    def findIngroupOutgroupSplit(self, node):
        potentialNodes = []
        for c in node.Children:
            if True == node.params.__contains__("ingroup"):
                if True == c.params.__contains__("outgroup"):
                    potentialNodes.append(node)
            elif True == node.params.__contains__("outgroup"):
                if True == c.params.__contains__("ingroup"):
                    potentialNodes.append(c)
            x = self.findIngroupOutgroupSplit(c)
            for n in x:
                potentialNodes.append(n)
        return potentialNodes
    
    #
    # Does this node contain any descendants with the 'mark'?
    #
    def containsMark(self, node, mark):
        if True == node.params.__contains__(mark):
            return True
        elif node.istip():
            return False
        else:
            for c in node.Children:
                if self.containsMark(c, mark):
                    return True
            return False
                
    #
    # Recur through the tree, return False if we find a node marked with both "ingroup"
    # and "outgroup"
    #
    def isIngroupOutgroupMonophyletic(self, node):
        if node.params.__contains__("ingroup") and node.params.__contains__("outgroup"):
            return False
        else:
            for c in node.Children:
                if False == self.isIngroupOutgroupMonophyletic(c):
                    return False
            return True
            
        
    #
    # Does a leaf taxa named 'taxaname' exist in the tree?
    # If yes, this method returns True
    # If no, this method returns False
    #
    def doesExtantTaxaExist(self, taxaname):
        self.root = LoadTree(treestring = self.newickString)
        c = self.root.getNodeNames()
        return c.__contains__(taxaname)
    
    # resets all the node params, including "ingroup" and "outgroup"
    def reset_node_params(self, node_names):
        for o in node_names:
            n = self.root.getNodeMatchingName(o)
            n.params = {}

    #
    # This method return the PhyloNode object corresponding to the last-common ancestor
    # of the nodes named in 'descendants'.  This method assumes that 'descendants' 
    # is a monophyletic group.
    #
    # INPUT:
    # 'descendants' is a list of node names.  This list must be complete! (i.e., if a taxon is not
    # in the list, then we'll assume it's in the outgroup!)
    #
    # OUTPUT:
    # a PhyloNode object (See the PyCogent API), or False if the LCA is ambiguous
    #
    def getLCA(self, descendants, outgroup):
        self.root = LoadTree(treestring = self.treeviewString)
        
        # reset any a priori labels:
        self.reset_node_params(descendants + outgroup)
        
        # Get the node objects corresponding to the ingroup...  
        descNodes = {}
        for o in descendants:
            n = self.root.getNodeMatchingName(o)
            n.params["ingroup"] = 1
            descNodes[o] = n
        # for each node in the ingroup....
        for n in descNodes:
            #print "tree.py 294: Starting mark chain at ", n
            # get the node's parent
            next = descNodes[n].Parent
            
            while next != None:
                if next.params.__contains__("ingroup"):
                    # we've interested another climbing path, so stop.
                    #print "tree.py 301: Stopping at " + next.Name
                    break
                else:
                    # mark the parent
                    next.params["ingroup"] = 1
                    #print "Marking node " + next.Name
                next = next.Parent
        
        # clean-up our markings, for the case where we marked too-high up the tree.
        self.unmarkRecur(self.root, "ingroup")
        

        # get the outgroup nodes:
        outNodes = {}
        for o in outgroup:
            n = self.root.getNodeMatchingName(o)
            n.params["outgroup"] = 1
            outNodes[o] = n
        # for each outgroup node, walk up the tree and mark nodes
        for n in outNodes:
            next = outNodes[n].Parent
            while next != None:
                if next.params.__contains__("outgroup"):
                    #if next.params["mark"] == 1:
                    break
                elif next.params.__contains__("outgroup"):
                    break
                else:
                    next.params["outgroup"] = 1
                    #print "Marking node " + next.Name
                next = next.Parent
                
        self.unmarkRecur(self.root, "outgroup")
        
        # at this point, the spanning trees for the ingroup and the outgroup subtrees
        # are marked.
        # Sanity Check:
        if False == self.isIngroupOutgroupMonophyletic(self.root):
            return False
            
        #s = self.findIngroupOutgroupSplit(self.root)
        return self.findRootDeep(self.root)
            #

            