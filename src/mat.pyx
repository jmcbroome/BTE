cimport mat
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef class MATNode:
    """
    A wrapper around the MAT node class. Has an identifier, mutations, parent, and child attributes.
    """
    cdef mat.Node* n

    cdef from_node(self, mat.Node* n):
        '''
        Load a node object from a MAT node. Attributes specific to the node such as mutations and identifier
        will be loaded automatically into python attributes and relational attributes to other nodes can be accessed through fetch methods 
        to the C++ class attributes.
        '''
        self.n = n

    def is_leaf(self):
        return self.n.is_leaf()

    def get_id(self):
        return self.n.identifier.decode("UTF-8")

    def get_parent(self):
        pn = MATNode()
        pn.from_node(self.n.parent)
        return pn

    def get_children(self):
        cnv = []
        for n in self.n.children:
            cn = MATNode()
            cn.from_node(n)
            cnv.append(cn)
        return cnv

    def get_mutations(self):
        return [m.get_string().decode("UTF-8") for m in self.n.mutations]

    def get_annotations(self):
        return [m.decode("UTF-8") for m in self.n.clade_annotations]

cdef class MATree:
    """
    A wrapper around the MAT Tree class. Includes functions to save and load from parsimony .pb files or a newick. Includes 
    numerous functions for tree traversal including both breadth-first and depth-first, the ability to search for a node by name and the
    ability to traverse from a specified node back to the root node.
    """
    cdef mat.Tree t

    def __init__(self, fpath=""):
        if fpath != "":
            if fpath[-3:] == ".pb" or fpath[-6:] == ".pb.gz":
                self.from_pb(fpath)
            elif fpath[-3:] == ".nwk":
                self.from_newick(fpath)
            else:
                raise Exception("Invalid file type. Must be .pb or .nwk")
        else:
            self.t = mat.Tree()

    def from_pb(self,file):
        self.t = mat.load_mutation_annotated_tree(file.encode("UTF-8"))

    def save_pb(self,file):
        mat.save_mutation_annotated_tree(self.t,file.encode("UTF-8"))

    def from_newick(self,nwk):
        self.t = mat.create_tree_from_newick(nwk.encode("UTF-8"))

    def get_parsimony_score(self):
        return self.t.get_parsimony_score()

    def get_node(self,name):
        return MATNode().from_node(self.t.get_node(name.encode("UTF-8")))

    cdef dfe_helper(self, mat.Node* node):
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.depth_first_expansion(node)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def depth_first_expansion(self, nid = ""):
        if nid == "":
            return self.dfe_helper(self.t.root)
        else:
            return self.dfe_helper(self.t.get_node(nid.encode("UTF-8")))

    cdef bfe_helper(self, string nid):
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.breadth_first_expansion(nid.encode("UTF-8"))
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def breadth_first_expansion(self,nid=""):
        return self.bfe_helper(nid)

    def get_newick_string(self,print_internal=False,print_branch_len=False,retain_original_branch_len=True,uncondense_leaves=False):
        return mat.get_newick_string(self.t,print_internal,print_branch_len,retain_original_branch_len,uncondense_leaves)

    cdef rsearch_helper(self, string nid, bool include_self):
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.rsearch(nid.encode("UTF-8"),include_self)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def rsearch(self,nid,include_self=False):
        return self.rsearch_helper(nid,include_self)

