cimport bte
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set as cset
from libcpp.map cimport map
from libc.stdint cimport *

cdef class MATNode:
    """
    A wrapper around the MAT node class. Has an identifier, mutations, parent, and child attributes.
    """
    cdef bte.Node* n

    cdef from_node(self,bte.Node* n):
        '''
        Load a node object from a MAT node. Attributes specific to the node such as mutations and identifier
        will be loaded automatically into python attributes and relational attributes to other nodes can be accessed through fetch methods 
        to the C++ class attributes.
        '''
        self.n = n
        return self

    def is_leaf(self):
        return self.n.is_leaf()

    @property
    def id(self):
        return self.n.identifier.decode("UTF-8")

    @property
    def parent(self):
        pn = MATNode()
        pn.from_node(self.n.parent)
        return pn

    @property
    def children(self):
        cnv = []
        for n in self.n.children:
            cn = MATNode()
            cn.from_node(n)
            cnv.append(cn)
        return cnv

    @property
    def mutations(self):
        return [m.get_string().decode("UTF-8") for m in self.n.mutations]

    @property
    def annotations(self):
        return [m.decode("UTF-8") for m in self.n.clade_annotations]

cdef complement(int8_t input):
    if input == 0b1:
        return 0b1000
    elif input == 0b1000:
        return 0b1
    elif input == 0b10:
        return 0b100
    elif input == 0b100:
        return 0b10
    else:
        return input

cdef class MATree:
    """
    A wrapper around the MAT Tree class. Includes functions to save and load from parsimony .pb files or a newick. Includes 
    numerous functions for tree traversal including both breadth-first and depth-first, the ability to search for a node by name and the
    ability to traverse from a specified node back to the root node.
    """
    cdef bte.Tree t

    def __init__(self, pbf=None, uncondense=True, nwk=None, vcf=None):
        if pbf != None:
            if nwk != None or vcf != None:
                print("WARNING: nwk and vcf arguments are exclusive with pbf. Ignoring nwk and vcf.")
            if pbf[-3:] == ".pb" or pbf[-6:] == ".pb.gz":
                self.from_pb(pbf,uncondense)
            else:
                raise Exception("Invalid file type. Must be .pb or .pb.gz")
        else:
            if nwk != None:
                if vcf == None:
                    raise Exception("Must provide a VCF file if loading from a newick.")
                self.from_newick_and_vcf(nwk,vcf)
            elif vcf != None:
                raise Exception("Loading from VCF requires a newick file.")
            else:
                self.t = bte.Tree()

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

    cdef assign_tree(self, bte.Tree t):
        self.t = t

    cdef resolve_all_polytomies(self):
        self.t = resolve_all_polytomies(self.t)

    def from_pb(self,file,uncondense=True):
        self.t = bte.load_mutation_annotated_tree(file.encode("UTF-8"))
        if uncondense:
            self.uncondense()
    
    def from_newick_and_vcf(self,nwk,vcf):
        self.t = bte.create_tree_from_newick(nwk.encode("UTF-8"))
        cdef vector[Missing_Sample] missing
        bte.read_vcf(&self.t,vcf.encode("UTF-8"),missing,False)

    def save_pb(self,file,condense=True):
        if condense:
            self.condense()
        bte.save_mutation_annotated_tree(self.t,file.encode("UTF-8"))
        #uncondense again afterwards if it was condensed for saving.
        if condense:
            self.uncondense()

    def from_newick(self,nwk):
        self.t = bte.create_tree_from_newick(nwk.encode("UTF-8"))

    def get_parsimony_score(self):
        return self.t.get_parsimony_score()

    def get_node(self,name):
        nc = MATNode()
        nc.from_node(self.t.get_node(name.encode("UTF-8")))
        return nc

    def get_root(self):
        nc = MATNode()
        nc.from_node(self.t.root)
        return nc

    def get_internal_node_descendents(self, name=""):
        cdef vector[string] all_desc = self.t.get_leaves_ids(name.encode("UTF-8"))
        for i in range(all_desc.size()):
            yield all_desc[i].decode("UTF-8")

    cdef dfe_helper(self, bte.Node* node):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.depth_first_expansion(node)
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

    def get_leaves(self, nid = ""):
        cdef vector[Node*] leaves = self.t.get_leaves(nid.encode("UTF-8"))
        wrappers = []
        for i in range(leaves.size()):
            nodec = MATNode().from_node(leaves[i])
            wrappers.append(nodec)
        return wrappers

    def get_leaves_ids(self, nid = ""):
        cdef vector[string] leaves = self.t.get_leaves_ids(nid.encode("UTF-8"))
        names = []
        for i in range(leaves.size()):
            names.append(leaves[i].decode("UTF-8"))
        return names

    cdef bfe_helper(self, string nid):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.breadth_first_expansion(nid.encode("UTF-8"))
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def breadth_first_expansion(self,nid=""):
        return self.bfe_helper(nid)

    def get_newick_string(self,print_internal=False,print_branch_len=False,retain_original_branch_len=True,uncondense_leaves=False):
        return bte.get_newick_string(self.t,print_internal,print_branch_len,retain_original_branch_len,uncondense_leaves)

    cdef rsearch_helper(self, string nid, bool include_self):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.rsearch(nid,include_self)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def rsearch(self,nid,include_self=False):
        return self.rsearch_helper(nid.encode("UTF-8"),include_self)

    cdef get_clade_samples(self, string clade_id):
        '''
        Return samples from the selected clade.
        '''
        cdef vector[string] samples = bte.get_clade_samples(&self.t, clade_id)
        return samples

    cdef get_mutation_samples(self, string mutation):
        '''
        Return samples containing the selected mutation.
        '''
        cdef vector[string] samples = bte.get_mutation_samples(&self.t, mutation)
        return samples

    cdef get_subtree(self, vector[string] samples):
        '''
        Return a subtree representing samples just from the selection.
        '''
        cdef bte.Tree subtree = bte.filter_master(self.t, samples, False, True)
        subt = MATree()
        subt.assign_tree(subtree)
        return subt
    
    def subtree(self, samples):
        cdef vector[string] samples_vec = [s.encode("UTF-8") for s in samples]
        return self.get_subtree(samples_vec)

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

    cdef get_regex_match(self, string regexstr):
        '''
        Return a subtree containing all samples matching the regular expression.
        '''
        #the C++ allows for a preselection of samples, but we don't use that option.
        cdef vector[string] to_check = []
        cdef vector[string] samples = bte.get_sample_match(&self.t, to_check, regexstr)
        return samples
    
    def get_regex(self, regexstr):
        cdef vector[string] samples = self.get_regex_match(regexstr.encode("UTF-8"))
        if samples.size() == 0:
            print("Error: requested regex does not match any samples.")
            return None
        print("Successfully found {} samples.".format(len(samples)))
        return self.get_subtree(samples)

    def get_random(self, size, current_samples = [], lca_limit = False):
        '''
        Select a random subtree of the selected size. Optionally, pass a list of samples to include.
        If the list of samples to include is larger than the target size, random samples will be removed from the list.
        Set lca_limit to True to limit random selection to below the common ancestor of the current selection.
        Selects as many as possible if not enough are available.
        '''
        if len(current_samples) == 0 and lca_limit:
            Exception("LCA limit requires a selection of samples to be passed in.")
        cdef vector[string] starting_samples = current_samples
        cdef size_t target_size = size
        cdef bool lcal = lca_limit
        cdef vector[string] final_samples = bte.fill_random_samples(&self.t, starting_samples, target_size, lcal)
        return self.get_subtree(final_samples)

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

    def count_leaves(self, subroot = ""):
        cdef Node* target_n = self.t.root
        if subroot != "":
            target_n = self.t.get_node(subroot.encode("UTF-8"))
        return self.t.get_num_leaves(target_n)

    cdef accumulate_mutations(self, string sample):
        '''
        Return the set of mutations the indicated sample has with respect to the root. 
        '''
        cdef vector[Node*] ancestors = self.t.rsearch(sample, True)
        allm = set()
        for i in range(ancestors.size()):
            #proceed in reverse order e.g. root to sample.
            for m in ancestors[ancestors.size()-i-1].mutations:
                #check that the opposite of this mutation is not in the set
                #if it is, instead delete it and skip this entry (as they negate each other with respect to the sample's final genome)
                mname = m.get_string().decode("UTF-8")
                oppo = mname[-1] + mname[1:-1] + mname[0]
                if oppo in allm:
                    allm.remove(oppo)
                else:
                    allm.add(mname)
        return allm

    def count_haplotypes(self):
        '''
        Return a dictionary of unique haplotypes and their counts. Used for nucleotide diversity estimates.
        '''
        cdef vector[string] samples = self.t.get_leaves_ids("".encode("UTF-8"))
        divtrack = {}
        for i in range(samples.size()):
            smset = self.accumulate_mutations(samples[i])
            smkey = tuple(sorted(smset))
            if smkey not in divtrack:
                divtrack[smkey] = 0
            divtrack[smkey] += 1
        return divtrack

    def compute_nucleotide_diversity(self):
        '''
        Cython-only function which computes the nucleotide diversity of the tree.
        This is defined as the mean number of pairwise differences in nucleotides between any two leaves
        of the tree. 
        '''
        def count_differences(l1,l2):
            total_matched = 0
            sl2 = set(l2)
            for element in l1:
                if element in sl2:
                    total_matched += 1
            return (len(l1)-total_matched) + (len(l2)-total_matched)
        #compute the set of mutations belonging to each sample in the tree, then compute pi from the resulting frequencies
        #since each sample is individually rsearched, this implementation is less efficient than an informed traversal, but still fast enough for most purposes.
        divtrack = self.count_haplotypes()
        total_seq = sum(divtrack.values())
        div = 0
        for g1 in divtrack.keys():
            g1_freq = divtrack[g1] / total_seq
            for g2 in divtrack.keys():
                if g1 != g2:
                    g2_freq = divtrack[g2] / total_seq
                    #pair_diff = len(set(g1).symmetric_difference(set(g2)))
                    pair_diff = count_differences(g1,g2)
                    div += pair_diff * g1_freq * g2_freq
        #multiply the final result to guarantee an unbiased estimator (see wikipedia entry?)
        return div * (total_seq / (total_seq - 1))

    def simple_parsimony(self, leaf_assignments):
        '''
        This function is an implementation of the small parsimony problem (Fitch algorithm) for a single set of states.
        It takes as input a dictionary mapping leaf names to character states and returns a dictionary mapping both leaf and internal node names to inferred character states.
        '''
        #this algorithm traverses the tree in postorder (reverse depth-first)
        #it requires that the tree be fully resolved and bifurcating, so that's the first step.
        #note on efficiency- we're using Python dicts to do most of the set logic, which are slower than using native C++ objects
        print("Resolving polytomies...")
        self.resolve_all_polytomies()
        #initialize node assignments with leaf states
        node_assignment_set = {l:set(v) for l,v in leaf_assignments.items()}
        print("{} initial assignments".format(len(node_assignment_set)))
        #first, a postorder traversal. We implement this here by generating nodes in depth-first order and proceeding 
        #to iterate through their indeces in reverse.
        cdef vector[Node*] nodes = self.t.depth_first_expansion(self.t.root)
        cdef Node* cnode
        cdef vector[Node*] children
        for i in range(nodes.size()):
            cnode = nodes[nodes.size()-i-1]
            if not cnode.is_leaf():
                #if it is correctly resolved, this node will have exactly two children.
                children = cnode.children
                assert children.size() == 2
                left = node_assignment_set[children[0].identifier.decode("UTF-8")]
                right = node_assignment_set[children[1].identifier.decode("UTF-8")]
                state = left & right
                if not state:
                    state = left | right
                node_assignment_set[cnode.identifier.decode("UTF-8")] = state
        #we then traverse again through the same sets of nodes in preorder/depth-first, finalizing character states based on 
        #the state of the parent.
        assert len(node_assignment_set) == nodes.size()
        final_node_assignment = leaf_assignments
        cdef Node* parent
        for i in range(nodes.size()):
            cnode = nodes[i]
            if not cnode.is_leaf():
                current_state = node_assignment_set[cnode.identifier.decode("UTF-8")]
                if len(current_state) == 1:
                    final_node_assignment[cnode.identifier.decode("UTF-8")] = list(current_state)[0]
                elif cnode.is_root():
                    final_node_assignment[cnode.identifier.decode("UTF-8")] = '0'
                else:
                    parent = cnode.parent
                    #in depth-first, this should always be accessible.
                    parent_state = final_node_assignment[parent.identifier.decode("UTF-8")]
                    if parent_state not in current_state:
                        #if the parent state is not part of the multiple options for this node, just take the first one.
                        final_node_assignment[cnode.identifier.decode("UTF-8")] = list(current_state)[0]
                    else:
                        final_node_assignment[cnode.identifier.decode("UTF-8")] = parent_state
        print("{} final assignments".format(len(final_node_assignment)))
        return final_node_assignment

    def ladderize(self):
        '''
        Sort the branches of the tree according to the size of each partition.
        '''
        self.t.rotate_for_consistency()

    def reverse_strand(self, genome_size = 29903):
        '''
        Inverts the tree representation of mutations such that all mutations are with respect to the reverse strand of the reference.
        All bases are complemented and indeces are reversed. The tree structure itself and parsimony scores are unaffected.
        '''
        cdef vector[Node*] nodes = self.t.depth_first_expansion(self.t.root)
        cdef vector[Mutation] nmv 
        cdef Mutation newmut
        for i in range(nodes.size()):
            nmv.clear()
            for mutation in nodes[i].mutations:
                newmut = mutation.copy()
                newmut.ref_nuc = complement(mutation.ref_nuc)
                newmut.par_nuc = complement(mutation.par_nuc)
                newmut.mut_nuc = complement(mutation.mut_nuc)
                newmut.position = genome_size - mutation.position - 1
                nmv.push_back(newmut)
            #since nodes[i] is a pointer to the original node object, we can simply edit it inplace.
            nodes[i].mutations = nmv
        
