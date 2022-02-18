cimport mat
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef class MATNode:
    cdef mat.Node* n

    cdef from_node(self, mat.Node* n):
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

cdef class MATree:
    cdef mat.Tree t

    def load_from(self,file):
        self.t = mat.load_mutation_annotated_tree(file.encode("UTF-8"))

    def save_to(self,file):
        mat.save_mutation_annotated_tree(self.t,file.encode("UTF-8"))

    def create_from_newick(self,nwk):
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
