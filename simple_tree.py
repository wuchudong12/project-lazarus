from cogent import *

class SimpleTree:
    def __init__(self):
        self.tree = None
        self.nodes = {}
        self.edges = {} # 2-d hash: (key) origin node name, (key) source node name, (value) edge object
        self.tips = {} # key = node name, value = SimpleNode object 

    def build_tree(self, cogent_tree):
        """Build a SimpleTree from a PyCogent tree"""
        pass
    
    def recursively_add_nodes(self, cogent_node):
        children = []
        if cogent_node.ittip
        for c in cogent_node.Children:
            children.append(c.Name)
            recursively_add_nodes(c)
        add_node(cogent_node.Name, cogent_node.Parent, children)
    
    def add_node(self, name, parent_name, children_name_array):
        """name is a string, parent is the name of a node, children_array is an array of node names."""
        if parent != None:
            p = self.nodes[parent]
        else:
            p = None
        c = []
        for n in children_name_array:
            c.append( self.nodes[n] )
        new_node = SimpleNode(name, p, c)
        self.nodes[name] = new_node
        
        if parent != None:
            new_edge = SimpleEdge(left_node = new_node, right_node = p)
            if False == self.edges.__contains__(name):
                self.edges[name] = {}
            self.edges[name][parent] = new_edge
            if False == self.edge.__contains__(parent):
                self.edges[parent] = {}
            self.edges[parent][name] = new_edge
        
    
    def reroot(self, new_root):
        """new_root can be a SimpleNode, in which case the rooted tree will have polytomy at the root, 
        or new_root can be a SimpleEdge, in which case the rooted tree will bifurcate at the root."""
        pass

class SimpleNode:
    def __init__(self, name = None, parent = None, children = [None]):
        self.name = name
        self.parent = parent
        self.children = children

class SimpleEdge:
    def __init__(self, left_node = None, right_node = None):
        self.left_node = left_node
        self.right_node = right_node