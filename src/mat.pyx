cimport mat
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair

class AAChange:
    """
    Pure python class containing translation information for a specific node's nucleotide mutation. Largely for convenience.
    """
    def __init__(self, gene, aa, nuc, codons):
        self.gene = gene        
        self.aa = aa
        self.original_aa = aa[0]
        self.aa_index = int(aa[1:-1])
        self.alternative_aa = aa[-1]
        self.nuc = nuc
        self.original_nt = nuc[0]
        self.nt_index = int(nuc[1:-1])
        self.alternative_nt = nuc[-1]
        self.mutation_type = nuc[0] + ">" + nuc[-1]
        self.original_codon = codons.split(">")[0]
        self.alternative_codon = codons.split(">")[1]

    def is_synonymous(self):
        if self.original_aa == self.alternative_aa:
            return True
        else:
            return False

cdef class MATNode:
    """
    A wrapper around the MAT node class. Has an identifier, mutations, parent, and child attributes.
    """
    cdef mat.Node* n
    cdef object translations

    def __init__(self):
        self.translations = []

    property translation:
        #read-only to public.
        def __get__(self):
            return self.translations

    cdef from_node(self, mat.Node* n):
        '''
        Load a node object from a MAT node. Attributes specific to the node such as mutations and identifier
        will be loaded automatically into python attributes and relational attributes to other nodes can be accessed through fetch methods 
        to the C++ class attributes.
        '''
        self.n = n

    def apply_translation(self, transtr):
        '''
        Store an amino acid translation string as a list attribute containing AAChange objects.
        '''
        #do nothing if an empty string is passed.
        if transtr != "":    
            aachanges = []
            #process the transtr in relevant pieces
            aagenes, nucs, codons = [v.split(";") for v in transtr.split('\t')]
            assert len(aagenes) == len(nucs)
            for i in range(len(aagenes)):
                gene, aa = aagenes[i].split(":")
                #sometimes two adjacent nucleotide mutations contribute to the same amino acid change, making them annoying to store correctly
                #make a copy of the change in another class object even though there's only one actual amino acid change here. Edge case to keep an eye on.
                mnucs = nucs[i].split(",")
                for mn in mnucs:
                    aachanges.append(AAChange(gene, aa, mn, codons[i]))
            self.translations = aachanges

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
    cdef object translation_table

    property translation_table:
        def __get__(self):
            return self.translation_table
        def __set__(self, value):
            self.translation_table = value

    def __init__(self, fpath="", uncondense=True):
        self.translation_table = {}
        if fpath != "":
            if fpath[-3:] == ".pb" or fpath[-6:] == ".pb.gz":
                self.from_pb(fpath,uncondense)
            elif fpath[-3:] == ".nwk":
                self.from_newick(fpath)
            else:
                raise Exception("Invalid file type. Must be .pb or .nwk")
        else:
            self.t = mat.Tree()

    cdef uncondense(self):
        '''
        Uncondense the tree.
        '''
        self.t.uncondense_leaves()

    cdef condense(self):
        '''
        Condense the tree.
        '''
        self.t.condense_leaves([])

    cdef assign_tree(self, mat.Tree t):
        self.t = t

    def from_pb(self,file,uncondense=True):
        self.t = mat.load_mutation_annotated_tree(file.encode("UTF-8"))
        if uncondense:
            self.uncondense()

    def save_pb(self,file,condense=True):
        if condense:
            self.condense()
        mat.save_mutation_annotated_tree(self.t,file.encode("UTF-8"))
        #uncondense again afterwards if it was condensed for saving.
        if condense:
            self.uncondense()

    def from_newick(self,nwk):
        self.t = mat.create_tree_from_newick(nwk.encode("UTF-8"))

    def get_parsimony_score(self):
        return self.t.get_parsimony_score()

    def get_node(self,name):
        nc = MATNode().from_node(self.t.get_node(name.encode("UTF-8")))
        if len(self.translation_table) > 0:
            print("Propagating translation.")
            nc.apply_translation(self.translation_table.get(name,""))
        return nc

    cdef dfe_helper(self, mat.Node* node, bool get_translation=False):
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.depth_first_expansion(node)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            if get_translation:
                nodec.apply_translation(self.translation_table.get(nvec[i].identifier.decode("UTF-8"),""))
            pynvec.append(nodec)
        return pynvec

    def depth_first_expansion(self, nid = ""):
        hastrans = len(self.translation_table) > 0
        if hastrans:
            print("Propagating translations to node search...")
        if nid == "":
            return self.dfe_helper(self.t.root,hastrans)
        else:
            return self.dfe_helper(self.t.get_node(nid.encode("UTF-8")),hastrans)

    cdef bfe_helper(self, string nid, bool hastrans):
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.breadth_first_expansion(nid.encode("UTF-8"))
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            if hastrans:
                nodec.apply_translation(self.translation_table.get(nvec[i].identifier.decode("UTF-8"),""))
            pynvec.append(nodec)
        return pynvec

    def breadth_first_expansion(self,nid=""):
        hastrans = len(self.translation_table) > 0
        if hastrans:
            print("Propagating translations to node search...")
        return self.bfe_helper(nid,hastrans)

    def get_newick_string(self,print_internal=False,print_branch_len=False,retain_original_branch_len=True,uncondense_leaves=False):
        return mat.get_newick_string(self.t,print_internal,print_branch_len,retain_original_branch_len,uncondense_leaves)

    cdef rsearch_helper(self, string nid, bool include_self):
        hastrans = len(self.translation_table) > 0
        if hastrans:
            print("Propagating translations to node search...")
        pynvec = []
        cdef vector[mat.Node*] nvec = self.t.rsearch(nid.encode("UTF-8"),include_self)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            if hastrans:
                nodec.apply_translation(self.translation_table.get(nvec[i].identifier.decode("UTF-8"),""))
            pynvec.append(nodec)
        return pynvec

    def rsearch(self,nid,include_self=False):
        return self.rsearch_helper(nid,include_self)

    cdef get_clade_samples(self, string clade_id):
        '''
        Return samples from the selected clade.
        '''
        cdef vector[string] samples = mat.get_clade_samples(&self.t, clade_id)
        return samples

    cdef get_mutation_samples(self, string mutation):
        '''
        Return samples containing the selected mutation.
        '''
        cdef vector[string] samples = mat.get_mutation_samples(&self.t, mutation.encode("UTF-8"))
        return samples

    cdef get_subtree(self, vector[string] samples):
        '''
        Return a subtree representing samples just from the selection. Does not retain translations computed on the original tree.
        '''
        cdef mat.Tree subtree = mat.filter_master(self.t, samples, False, True)
        subt = MATree()
        subt.assign_tree(subtree)
        return subt

    def get_clade(self, clade_id):
        '''
        Return a subtree representing the selected clade.
        '''
        print("Getting clade: " + clade_id)
        cdef vector[string] samples = self.get_clade_samples(clade_id.encode("UTF-8"))
        if samples.size() == 0:
            print("Error: requested clade not found.")
            return None
        print("Successfully found {} samples.".format(len(samples)))
        return self.get_subtree(samples)

    def with_mutation(self, mutation):
        '''
        Return a subtree of samples containing the indicated mutation.
        '''
        cdef vector[string] samples = self.get_mutation_samples(mutation.encode("UTF-8"))
        return self.get_subtree(samples)

    def count_mutations(self, subroot=""):
        '''
        Cython-only function which computes the counts of individual mutation types across the tree. Can take a specific node to start from.
        '''
        mcount = {}
        cdef Node* target_n = self.t.root
        if subroot != "":
            target_n = self.t.get_node(subroot.encode("UTF-8"))
        cdef vector[Node*] nvec = self.t.depth_first_expansion(target_n)
        for i in range(nvec.size()):
            for m in nvec[i].mutations:
                pym = m.get_string().decode("UTF-8")
                pym_type = pym[0] + ">" + pym[-1]
                if pym_type in mcount:
                    mcount[pym_type] += 1
                else:
                    mcount[pym_type] = 1
        return mcount

    def translate(self,gtf_file,fasta_file):
        """
        Translate amino acid changes across the tree and return the results as a dictionary of node IDs and change strings as returned from matUtils translate. 
        """
        cdef vector[pair[string,string]] changes = mat.do_translation(&self.t,gtf_file.encode("UTF-8"),fasta_file.encode("UTF-8"))
        for i in range(changes.size()):
            self.translation_table[changes[i].first.decode("UTF-8")] = changes[i].second.decode("UTF-8")
