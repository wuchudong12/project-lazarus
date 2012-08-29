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

def test_addNode():
     print "\n--> test_addNode:"   
     t = Tree(1)
     
     # 
     sys.stdout.write("test 1. . .") 
     t.addNode(123)
     n = t.nodes[123].nodeNumber
     if n == 123:
         print "SUCCESS."
     else:
        print "failed!"
        print "This test returned: " + n.__str__()
        raise AssertionError

def test_site_fromString():
    print "\n--> test_fromString:"
    
    sys.stdout.write("test 1. . .")
    l = "1  A 1.000"
    s = Site(None)
    s.fromString(l)
    if False == s.stateProbs.__contains__("A"):
        print "failed!"
        raise AssertionError   
    if s.stateProbs["A"] != 1.0:
        print "failed!"
        raise AssertionError
    print "SUCCESS."
    
    sys.stdout.write("test 2. . .")
    l = "5  A 1.000 B 0.9 C 2.5"
    s = Site(None)
    s.fromString(l)
    if False == s.stateProbs.__contains__("A"):
        print "failed!"
        raise AssertionError   
    if s.stateProbs["A"] != 1.0:
        print "failed!"
        raise AssertionError
    if False == s.stateProbs.__contains__("B"):
        print "failed!"
        raise AssertionError
    if s.stateProbs["B"] != 0.9:
        print "failed!"
        raise AssertionError
    if False == s.stateProbs.__contains__("C"):
        print "failed!"
        raise AssertionError
    if s.stateProbs["C"] != 2.5:
        print "failed!"
        raise AssertionError
    print "SUCCESS."


def test_getAncestralNode():
    print "\n--> test_getAncestralNode:"
    
    t = Tree(1)
    t.newickString = "((A,B),(C,D));"
    t.treeviewString = "((A,B)E,(C,D)F)G;" 
    
    # the ingroup and outgroup are both monophyletic
    sys.stdout.write("test 1. . .")
    x = t.getAncestralNode(["A", "B"], ["C", "D"])
    if x != "E":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
    
    # the ingroup and outgroup are not monophyletic
    sys.stdout.write("test 2. . .")
    x = t.getAncestralNode(["A", "C"], ["B", "D"])
    if x == False:
        print "SUCCESS."
    else:
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
        
    # the ingroup straddles the root, the outgroup is monophyletic
    sys.stdout.write("test 3. . .")
    x = t.getAncestralNode(["A", "C"], ["D"])
    print x.Name
    if x != "F":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    t.treeviewString = "((A,(B1,B2)B)E,(C,D)F)G;" 
    # the outgroup straddles the root, the ingroup is monophyletic
    sys.stdout.write("test 4. . .")
    x = t.getAncestralNode(["B1", "B2"], ["A", "D"])
    if x != "B":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."

def test_getLCA():
    t = Tree(1)
    
    print "\n--> test_LCA:"
    
    t = Tree(1)
    t.newickString = "((A,B),(C,D));"
    t.treeviewString = "((A,B)E,(C,D)F)G;"

    # the ingroup and outgroup are both monophyletic
    sys.stdout.write("test 1. . .")
    x = t.getLCA(["C", "D"], ["A", "B"])
    if x == False or x.Name != "F":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    sys.stdout.write("test 2. . .")
    x = t.getLCA(["A", "C"], ["B"])
    print x.Name
    if x == False or x.Name != "E":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    sys.stdout.write("test 3. . .")
    x = t.getLCA(["A", "B", "C"], ["D"])
    if x == False or x.Name != "F":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    sys.stdout.write("test 4. . .")
    t.newickString = "(A,B,C,D);"
    t.treeviewString = "(A,B,C,D)E;"
    x = t.getLCA(["A", "B", "C"], ["D"])
    if x == False or x.Name != "E":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    sys.stdout.write("test 5. . .")
    x = t.getLCA(["A", "B"], ["C"])
    if x == False or x.Name != "E":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."
        
    sys.stdout.write("test 6. . .")
    t.newickString = "((((A,B)),C), (D,E));"
    t.treeviewString = "((((A,B)AB),C)ABC, (D,E)DE)ABCDE;"
    x = t.getLCA(["A", "C"], ["B"])
    if x == False or x.Name != "AB":
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."

    sys.stdout.write("test 7. . .")
    t.newickString = "((((A,B)),C), (D,E));"
    t.treeviewString = "((((A,B)AB),C)ABC, (D,E)DE)ABCDE;"
    x = t.getLCA(["C", "B"], ["A", "E"])
    if x != False:
        print "failed!"
        print "This test returned: " + x.__str__()
        raise AssertionError
    else:
        print "SUCCESS."



#
# main(). . .
#
def test_tree():
    test_site_fromString()
    test_addNode()
    test_getLCA()
    # getAncestralNode is depricated. . .
    #test_getAncestralNode()   

    
