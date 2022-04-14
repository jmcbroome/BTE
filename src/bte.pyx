#cython: embedsignature=True
#cython: annotation_typing=True
cimport bte
import cython
from cython.operator cimport dereference, postincrement, preincrement
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set as cset
from libcpp.map cimport map
from libc.stdint cimport *
from libcpp cimport bool as cbool
import functools
import time
from typing import Optional, Union

def _timer(func, *args, **kwargs):
    @functools.wraps(func)
    def timer_wrap(*args,**kwargs):
        start = time.perf_counter()
        retval = func(*args,**kwargs)
        end = time.perf_counter()
        rt = end-start
        print(f"Finished {func.__name__!r} in {round(rt,4)} seconds")
        return retval
    return timer_wrap

cdef class MATNode:
    """
    A wrapper around the MAT node class. Has an identifier, mutations, parent, and child attributes.
    """
    cdef bte.Node* n

    cdef from_node(self, bte.Node* n):
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
        cdef bte.Node* ancestor = self.n
        cdef size_t anncount = self.n.clade_annotations.size()
        cdef vector[string] annotes 
        cdef size_t k
        annotations = [None for i in range(anncount)]
        while all([a==None for a in annotations]):
            annotes = ancestor.clade_annotations
            for k in range(annotes.size()):
                if annotes[k].size() > 0 and (annotations[k] == None):
                    annotations[k] = annotes[k].decode("UTF-8")
            if ancestor.parent == cython.NULL:
                break
            ancestor = ancestor.parent
        return annotations
            
    def update_mutations(self, mutation_list: list[str]):
        """Take a list of mutations as strings and replace any currently stored mutations on this branch with the new set.
        Mutation strings should be formatted as chro:reflocalt e.g. chr1:A234G. If chromosome is left off, assumes SARS-CoV-2 chromosome.

        Args:
            mutation_list (list[str]): List of mutations to store.
        """        
        self.n.mutations.clear()
        cdef bte.Mutation newmut
        cdef int8_t pn, mn
        cdef char pns, mns
        for mstr in mutation_list:
            if ':' in mstr:
                chro, info = mstr.split(":")
            else:
                chro = "NC_045512v2"
                info = mstr
            assert len(mstr) > 0
            loc = int(mstr[1:-1])
            newmut.chrom = chro.encode("UTF-8")
            #ref_nuc is disregarded with this loading strategy.
            pns = ord(mstr[0])
            pn = bte.get_nuc_id(pns)
            newmut.par_nuc = pn
            
            mns = ord(mstr[-1])
            mn = bte.get_nuc_id(mns)
            newmut.mut_nuc = mn

            newmut.position = loc
            self.n.mutations.push_back(newmut.copy())

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
    numerous functions for tree traversal including breadth-first, depth-first, and traversal from leaf to roots. Also includes
    numerous functions for subtree selection by choosing leaves that match regex patterns, contain specific mutations, or 
    are from a specific clade or lineage.
    """
    cdef bte.Tree t
    cdef public cbool _tree_only
    cdef public cbool _empty

    @_timer
    def __init__(self, pb_file: Optional[str] = None, uncondense: bool = True, nwk_file: Optional[str] = None, nwk_string: Optional[str] = None, vcf_file: Optional[str] = None, json_file: Optional[str] = None) -> None:
        """Instantiate a MATree object.

        Args:
            pb_file (Optional[str], optional): Load from a protobuf. Defaults to None.
            uncondense (bool, optional): Uncondense after loading from a protobuf. Defaults to True.
            nwk_file (Optional[str], optional): Load from a text file containing a newick representation of a tree. Defaults to None.
            nwk_string (Optional[str], optional): Load from a Python string newick. Defaults to None.
            vcf_file (Optional[str], optional): Load genotype information for leaves from a vcf. Defaults to None.
            json_file (Optional[str], optional): Load a complete tree from an Auspice-format JSON. Defaults to None.

        Raises:
            Exception: Invalid file type.
            Exception: Loading from VCF without a newick.
        """        
        self._tree_only = False
        self._empty = False
        if pb_file != None:
            if nwk_file != None or vcf_file != None or nwk_string != None:
                print("WARNING: nwk_file, nwk_string and vcf_file arguments are exclusive with pbf. Ignoring")
            if pb_file[-3:] == ".pb" or pb_file[-6:] == ".pb.gz":
                self.from_pb(pb_file,uncondense)
            else:
                raise Exception("Invalid file type for pb_file argument. Must be .pb or .pb.gz")
        elif json_file != None:
            self.from_json(json_file)
        elif nwk_string != None:
            if vcf_file != None:
                print("WARNING: nwk_string is for tree-only loading and does not use vcf input. Consider using nwk_file.")
            self.t = bte.create_tree_from_newick_string(nwk_string.encode("UTF-8"))
            self._tree_only = True
        else:
            if nwk_file != None:
                if vcf_file == None:
                    self.t = bte.create_tree_from_newick(nwk_file.encode("UTF-8"))
                    self._tree_only = True
                else:
                    self.from_newick_and_vcf(nwk_file,vcf_file)
            elif vcf_file != None:
                raise Exception("Loading from VCF requires a newick file.")
            else:
                self.t = bte.Tree()
    
    def clear(self) -> None:
        '''
        Call this function to explicitly deallocate all tree memory. Use when the tree object is no longer necessary 
        and high memory use is becoming problematic. Automatically called on garbage collection.
        '''
        cdef vector[Node*] nodes 
        if not self._empty:
            nodes = self.t.depth_first_expansion(self.t.root)
            for i in range(nodes.size()):
                nodes[i].mutations.clear()
            bte.clear_tree(self.t)
            self._empty = True

    def __del__(self):
        self.clear()

    @staticmethod
    def _check_newick_only(func, *args, **kwargs):
        '''
        Decorator function applied to class functions which require mutations to exist on the tree. 
        A tree loaded from a newick alone will be unable to apply these functions and this decorator serves
        as a readable way to add checks to these functions.
        '''
        @functools.wraps(func)
        def wrap(self, *args, **kwargs):
            if self._tree_only:
                raise Exception("Tree does not contain explicit mutation information and this function cannot be used. You can add mutations with apply_mutations or load from a vcf or pb.")
            else:
                return func(self, *args,**kwargs)
        return wrap

    def apply_mutations(self, mmap: dict[str,list[str]]) -> None:
        """Pass a set of node:mutation mappings to place into the tree. Current mutations will be replaced.

        Args:
            mmap (dict[str,list[str]]): A dictionary of node:mutation list mappings (e.g. {"node_id":["chro:reflocalt","chro:reflocalt"]}, {"node_1":["chro1:A123G","chro3:T315G"]}
        """        
        for nid, nms in mmap.items():
            node = self.get_node(nid)
            node.update_mutations(nms)
        #the tree is now annotated with mutations and mutation-based functions can be attempted.
        self._tree_only = False

    cdef uncondense(self):
        self.t.uncondense_leaves()

    cdef condense(self):
        self.t.condense_leaves([])

    cdef assign_tree(self, bte.Tree t):
        self.t = t

    cdef resolve_all_polytomies(self):
        self.t = resolve_all_polytomies(self.t)

    @_timer
    def from_pb(self, file: str, uncondense: bool = True) -> None:
        """Load from a protobuf into the initalized wrapper. Includes both tree and mutation information.

        Args:
            file (str): Path to a .pb or .pb.gz file.
            uncondense (bool, optional): Uncondense the tree after loading (split identical samples into individual leaves). Defaults to True.
        """
        cdef string fn = file.encode("UTF-8")
        cdef bte.Tree lt
        with nogil:
            lt = bte.load_mutation_annotated_tree(fn)
        self.t = lt
        if uncondense:
            self.uncondense()
        self._tree_only = False

    @_timer    
    def from_newick_and_vcf(self, nwk: str, vcf: str) -> None:
        """Load from a newick and a vcf. The vcf must contain sample entries (genotype columns) for every leaf in the newick.

        Args:
            nwk (str): Path to a text file containing the tree to load in Newick format.
            vcf (str): Path to a text file containing leaf/sample genotype information in VCF format.
        """        

        self.t = bte.create_tree_from_newick(nwk.encode("UTF-8"))
        cdef vector[Missing_Sample] missing
        bte.read_vcf(&self.t,vcf.encode("UTF-8"),missing,True)
        self._tree_only = False

    def save_pb(self, file: str, condense: bool = True) -> None:
        """Save the tree to a protobuf file.

        Args:
            file (str): Name for the .pb file.
            condense (bool, optional): Condense the tree before saving. Defaults to True.
        """        
        if condense:
            self.condense()
        bte.save_mutation_annotated_tree(self.t,file.encode("UTF-8"))
        #uncondense again afterwards if it was condensed for saving.
        if condense:
            self.uncondense()

    @_timer
    def from_newick(self, nwk_file: str) -> None:
        """Load from a newick file only. The resulting tree will lack mutation information, preventing some functions from being applied.

        Args:
            nwk_file (str): Path to a text file containing a newick representation of the tree.
        """        
        self.t = bte.create_tree_from_newick(nwk_file.encode("UTF-8"))
        self._tree_only = True

    def from_newick_string(self, nwk: str) -> None:
        """Load from a Python string newick. The resulting tree will lack mutation information, preventing some functions from being applied.

        Args:
            nwk (str): A Python string containing a newick representation of the tree.
        """        
        self.t = bte.create_tree_from_newick_string(nwk.encode("UTF-8"))
        self._tree_only = True

    def write_newick(self, subroot: Optional[str] = None, print_internal: bool = True, print_branch_len: bool = True, retain_original_branch_len: bool = True, uncondense_leaves: bool = True) -> str:
        """Write a newick string from the tree.

        Args:
            subroot (Optional[str], optional): Write a newick representing the subtree descended from this node. Defaults to the root.
            print_internal (bool, optional): Include internal node names in the newick output. Defaults to True.
            print_branch_len (bool, optional): Print branch lengths. Defaults to True.
            retain_original_branch_len (bool, optional): Retain the original branch length attribute, if one was provided. Defaults to True.
            uncondense_leaves (bool, optional): Uncondense nodes before writing the newick. Defaults to True.

        Returns:
            str: A newick string representation of the tree.
        """        
        cdef stringstream ss
        cdef Node* sr
        if subroot == None:
            sr = self.t.root
        else:
            se = self.t.get_node(subroot.encode("UTF-8"))
        bte.write_newick_string(ss,self.t,sr,print_internal,print_branch_len,retain_original_branch_len,uncondense_leaves)
        return ss.to_string().decode("UTF-8")

    cpdef vector[string] _read_samples(self,samples):
        '''
        Helper function which converts a Python list of bytes or string samples to a usable C++ string vector.
        Returns all leaf names if samples is empty.
        '''
        cdef vector[string] sample_names
        if len(samples) == 0:
            sample_names = self.t.get_leaves_ids("".encode("UTF-8"))
        else:
            for s in samples:
                if type(s) == str:
                    sample_names.push_back(s.encode("UTF-8"))
                elif type(s) == bytes:
                    sample_names.push_back(s)
                else:
                    raise Exception("Sample names must be strings or bytes.")
        return sample_names

    def write_vcf(self, vcf_file: str , no_genotypes: bool = False, samples: list[str] = []) -> None:        
        """Write a vcf representing the chosen samples to the indicated file. By default, writes a vcf including all samples.

        Args:
            vcf_file (str): Name the output vcf file.
            no_genotypes (bool, optional): Do not include individual genotype information in the output vcf. Defaults to False.
            samples (list[str], optional): Samples to include. Defaults to all samples.
        """        
        sample_names = self._read_samples(samples)
        bte.make_vcf(self.t,vcf_file.encode("UTF-8"),no_genotypes,sample_names)

    def from_json(self, jsonf: str) -> None:
        """Load a mat from a json compatible with the Auspice.us visualization web tool.

        Args:
            jsonf (str): Path to a json file.
        """        
        self.t = bte.load_mat_from_json(jsonf.encode("UTF-8"))

    def write_json(self, jsonf: str, samples: list[str] = [], title: str = "Tree", metafiles: list[str] = []) -> None:
        """Write a json compatible with the Auspice.us visualization web tool containing the indicated samples. Default behavior includes the whole tree.
        You can optionally pass a tsv or csv file or a list of tsv and csv files containing categorical metadata to decorate the json with (one sample per row).

        Args:
            jsonf (str): Name for the JSON output.
            samples (list, optional): Samples to use. Defaults to all samples.
            title (str, optional): Title of the JSON. Defaults to "Tree".
            metafiles (list, optional): Metadata tsv and csv files to use. Defaults to no metadata.
        """
        cdef vector[string] sample_names = self._read_samples(samples)
        cdef cset[string] sample_set
        for i in range(sample_names.size()):
            sample_set.insert(sample_names[i])
        cdef vector[unordered_map[string,unordered_map[string,string]]] catmeta
        if type(metafiles) == str:
            catmeta.push_back(bte.read_metafile(metafiles.encode("UTF-8"),sample_set))
        elif type(metafiles) == list and len(metafiles) > 0:
            for mf in metafiles:
                catmeta.push_back(bte.read_metafile(mf.encode("UTF-8"),sample_set))
        bte.write_json_from_mat(&self.t,jsonf.encode("UTF-8"),&catmeta,title.encode("UTF-8"))

    def get_parsimony_score(self) -> int:
        '''
        Compute the parsimony score of the complete tree.

        Returns:
            int: The parsimony score of the tree.
        '''
        return self.t.get_parsimony_score()

    def get_node(self, name: str) -> MATNode:
        """Create a MATNode class object representing the indicated node.

        Args:
            name (str): ID of the node to fetch.

        Returns:
            Node: MATNode class object representing the indicated node.
        """
        nc = MATNode()
        nc.from_node(self.t.get_node(name.encode("UTF-8")))
        return nc

    @property
    def root(self) -> MATNode:
        """Return a MATNode wrapper representing the root node of the tree.

        Returns:
            MATNode: MATNode wrapper representing the root node of the tree.
        """        
        nc = MATNode()
        nc.from_node(self.t.root)
        return nc

    cdef dfe_helper(self, bte.Node* node, cbool reverse):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.depth_first_expansion(node)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        if reverse:
            pynvec.reverse()
        return pynvec

    @_timer
    def depth_first_expansion(self, nid: Optional[str] = None, reverse: bool = False) -> list[MATNode]:
        """Perform a preorder (depth-first) expansion of the tree, starting from the indicated node. 
        By default, traverses the whole tree. Set reverse to true to traverse in postorder (reverse depth-first) instead.

        Args:
            nid (Optional[str], optional): Node to begin the traversal at. Defaults to the root.
            reverse (bool, optional): Traverse in reverse order. Defaults to False.

        Returns:
            list[MATNode]: List of 
        """        

        if nid == None:
            return self.dfe_helper(self.t.root, reverse)
        else:
            return self.dfe_helper(self.t.get_node(nid.encode("UTF-8")), reverse)

    @_timer
    def get_leaves(self, nid: str = "") -> list[MATNode]:
        """Create a list of MATNode objects representing each leaf descended from the indicated node. By default, returns all leaves on the tree.

        Args:
            nid (str, optional): Node to get leaves descended from. Defaults to the root.

        Returns:
            list[MATNode]: List of MATNode wrappers representing all leaves.
        """        
        cdef vector[Node*] leaves = self.t.get_leaves(nid.encode("UTF-8"))
        wrappers = []
        for i in range(leaves.size()):
            nodec = MATNode().from_node(leaves[i])
            wrappers.append(nodec)
        return wrappers

    @_timer
    def get_leaves_ids(self, nid: str = "") -> list[str]:
        """Return a list of leaf name strings containing all leaves descended from the indicated node. By default, returns all leaves on the tree.

        Args:
            nid (str, optional): Node to get the descendents of. Defaults to the root.

        Returns:
            list[str]: List of leaf names.
        """        

        cdef vector[string] leaves = self.t.get_leaves_ids(nid.encode("UTF-8"))
        names = []
        for i in range(leaves.size()):
            names.append(leaves[i].decode("UTF-8"))
        return names

    cdef bfe_helper(self, string nid, cbool reverse):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.breadth_first_expansion(nid)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        if reverse:
            pynvec.reverse()
        return pynvec

    @_timer
    def breadth_first_expansion(self, nid: str = "", reverse: bool = False) -> list[MATNode]:
        """Perform a level order (breadth-first) expansion starting from the indicated node. Use reverse to traverse in reverse level order (all leaves, then all leaf parents, back to root) instead.

        Args:
            nid (str, optional): Node to begin the traversal at. Defaults to the root.
            reverse (bool, optional): Perform the traversal in reverse. Defaults to False.

        Returns:
            list[MATNode]: List of MATNode wrappers representing the nodes in the traversal.
        """        
        return self.bfe_helper(nid.encode("UTF-8"),reverse)

    cdef rsearch_helper(self, string nid, cbool include_self, cbool reverse):
        pynvec = []
        cdef vector[bte.Node*] nvec = self.t.rsearch(nid,include_self)
        for i in range(nvec.size()):
            nodec = MATNode()
            nodec.from_node(nvec[i])
            pynvec.append(nodec)
        return pynvec

    def rsearch(self, nid: str, include_self: bool = False,reverse: bool = False) -> list[MATNode]:
        """Return a list of MATNode objects representing the ancestors of the indicated node back to the root in order from the node to the root.

        Args:
            nid (str): ID of the node to get the ancestry of.
            include_self (bool, optional): Include the indicated node on the path. Defaults to False.
            reverse (bool, optional): Return the path in reverse order. Defaults to False.

        Returns:
            list[MATNode]: A list of ancestors of the indicated node.

        """ 
        return self.rsearch_helper(nid.encode("UTF-8"),include_self,reverse)

    cdef get_subtree(self, vector[string] samples):
        cdef bte.Tree subtree
        with nogil:
            subtree = bte.filter_master(self.t, samples, False, True)
        subt = MATree()
        subt.assign_tree(subtree)
        return subt

    @_timer    
    def subtree(self, samples: list[Union[str,bytes]]) -> MATree:
        """Retrieve a subtree containing all samples in the input list.

        Args:
            samples (list[Union[str,bytes]]): List of sample names to include in the subtree. Can be bytes or str.

        Returns:
            MATree: the subtree containing all samples in the input list.
        """        
        cdef vector[string] samples_vec = self._read_samples(samples)
        return self.get_subtree(samples_vec)

    cpdef vector[string] get_clade_samples(self, clade_id):
        '''
        Return samples from the selected clade.
        '''
        cdef vector[string] samples = bte.get_clade_samples(&self.t, clade_id.encode("UTF-8"))
        return samples

    def get_clade(self, clade_id: str) -> MATree:
        """Return a subtree representing the selected clade.

        Args:
            clade_id (str): The clade to retrieve.

        Returns:
            MATree: the subtree representing that clade.
        """
        print("Getting clade: " + clade_id)
        cdef vector[string] samples = self.get_clade_samples(clade_id)
        if samples.size() == 0:
            print("Error: requested clade not found.")
            return None
        print("Successfully found {} samples.".format(len(samples)))
        return self.get_subtree(samples)

    cpdef vector[string] get_regex_samples(self, regexstr):
        #the C++ allows for a preselection of samples, but we don't use that option.
        cdef vector[string] to_check = []
        cdef vector[string] samples = bte.get_sample_match(&self.t, to_check, regexstr.encode("UTF-8"))
        return samples
    
    def get_regex(self, regexstr: str) -> MATree:
        """Return a subtree representing all samples matching the regular expression.

        Args:
            regexstr (str): The regular expression to use to query the tree.

        Returns:
            MATree: Subtree containing samples matching the regex.
        """
        cdef vector[string] samples = self.get_regex_samples(regexstr)
        if samples.size() == 0:
            print("Error: requested regex does not match any samples.")
            return None
        print("Successfully found {} samples.".format(len(samples)))
        return self.get_subtree(samples)

    def get_random(self, size: int, current_samples: list = [], lca_limit: bool = False) -> MATree:
        """Select a random subtree of the selected size. Optionally, pass a list of samples to include.
        If the list of samples to include is larger than the target size, random samples will be removed from the list.
        Set lca_limit to True to limit random selection to below the common ancestor of the current selection.
        Selects as many as possible if not enough are available.

        Args:
            size (int): The size of the subtree to select.
            current_samples (list, optional): List of samples to include in the set. Defaults to [].
            lca_limit (bool, optional): Limit randomly selected samples to below the LCA of the input samples. Defaults to False.

        Returns:
            MATree: A subtree containing the selected samples.
        """        
        if len(current_samples) == 0 and lca_limit:
            Exception("LCA limit requires a selection of samples to be passed in.")
        cdef vector[string] starting_samples = current_samples
        cdef size_t target_size = size
        cdef bool lcal = lca_limit
        cdef vector[string] final_samples = bte.fill_random_samples(&self.t, starting_samples, target_size, lcal)
        return self.get_subtree(final_samples)

    cpdef vector[string] get_mutation_samples(self, mutation):
        """Return samples with genotypes containing the selected mutation.

        Args:
            mutation (str): string representation of the mutation in reflocalt format (e.g. "A123C").

        Returns:
            list[str]: samples with the mutation.
        """        
        #can't use the error decorator for this function since it is cpdef.
        if self._tree_only:
            Exception("Cannot find mutations on a tree-only MATree.")
        cdef vector[string] samples = bte.get_mutation_samples(&self.t, mutation.encode("UTF-8"))
        return samples

    @_check_newick_only
    def get_mutation(self, mutation: str) -> MATree:
        """Return a subtree containing samples with genotypes containing the indicated mutation.

        Args:
            mutation (str): string representation of the mutation in reflocalt format (e.g. "A123C").

        Returns:
            MATree: subtree containing samples with the mutation.
        """        
        cdef vector[string] samples = self.get_mutation_samples(mutation)
        return self.get_subtree(samples)

    @_check_newick_only
    def count_mutation_types(self, subroot: Optional[str] = None) -> dict[str,int]:
        """Compute the counts of individual mutation types across the tree. If a subtree root is indicated, it only counts mutations
        descended from that node. By default, this counts across the entire tree.

        Args:
            subroot (Optional[str], optional): Count mutations descended from the indicated node. Defaults to the root.

        Returns:
            dict[str,int]: Dictionary containing mutation counts.
        """
        mcount = {}
        cdef Node* target_n = self.t.root
        if subroot != None:
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

    def count_leaves(self, subroot: Optional[str] = None) -> int:
        """Return the number of leaves descended from the indicated node. By default, counts all leaves on the tree.

        Args:
            subroot (Optional[str], optional): Count leaves descended from the indicated node. Defaults to the root.

        Returns:
            int: The count of leaves.
        """        
        cdef Node* target_n = self.t.root
        if subroot != None:
            target_n = self.t.get_node(subroot.encode("UTF-8"))
        return self.t.get_num_leaves(target_n)

    cdef cset[bte.Mutation] accumulate_mutations(self, string sample):
        cdef vector[bte.Node*] ancestors = self.t.rsearch(sample, True)
        cdef cset[bte.Mutation] allm
        cdef size_t i
        cdef bte.Node* cnode
        cdef bte.Mutation m, om
        cdef cset[bte.Mutation].iterator oml
        cdef size_t ancs = ancestors.size()
        for i in range(ancs):
            #proceed in reverse order e.g. root to sample.
            cnode = ancestors[ancs-i-1]
            for j in range(cnode.mutations.size()):
                m = cnode.mutations[j]
                #check that the opposite of this mutation is not in the set
                #if it is, instead delete it and skip this entry (as they negate each other with respect to the sample's final genome)
                om = m.copy()
                om.mut_nuc = m.par_nuc
                om.par_nuc = m.mut_nuc
                oml = allm.find(om)
                if oml != allm.end():
                    allm.erase(oml)
                else:
                    allm.insert(m)
        return allm

    @_check_newick_only
    def mutation_set(self, nid: str) -> set[str]:
        """Return the complete set of mutations (genotype) the indicated node has with respect to the reference. 

        Args:
            nid (str): The target node to get the genotype for.

        Returns:
            set[str]: the genotype of the node.
        """
        pyset = set()
        cdef cset[bte.Mutation] accm = self.accumulate_mutations(nid.encode("UTF-8"))
        cdef cset[bte.Mutation].iterator accm_it = accm.begin()
        while accm_it != accm.end():
            pyset.add(dereference(accm_it).get_string())
            postincrement(accm_it)
        return pyset

    cdef count_differences(self, cset[bte.Mutation] l1, cset[bte.Mutation] l2):
        cdef size_t total_matched = 0
        cdef bte.Mutation e1
        cdef cset[bte.Mutation].iterator it = l1.begin()
        while it != l1.end():
            e1 = dereference(it)
            if l2.find(e1) != l2.end():
                total_matched += 1
            postincrement(it)
        return l1.size() + l2.size() - (total_matched * 2)

    cdef map[cset[bte.Mutation],size_t] count_haplotypes_c(self):
        '''
        Return a dictionary of unique haplotypes and their counts. Used for nucleotide diversity estimates.
        '''
        cdef vector[string] samples = self.t.get_leaves_ids("".encode("UTF-8"))
        cdef map[cset[bte.Mutation],size_t] divtrack
        cdef map[cset[bte.Mutation],size_t].iterator finder
        cdef size_t i
        cdef cset[bte.Mutation] accum_muts
        for i in range(samples.size()):
            accum_muts = self.accumulate_mutations(samples[i])
            finder = divtrack.find(accum_muts)
            if finder == divtrack.end():
                divtrack[accum_muts] = 1
            else:
                divtrack[accum_muts] = dereference(finder).second + 1
        return divtrack

    def count_clades_inclusive(self, subroot: str = "") -> dict[str,int]:
        """Count the total number of leaves belonging to each clade on the subtree.
        Counts are inclusive (e.g. samples belonging to a clade descended from another clade will count for the ancestor clade as well)
        By default, counts across the whole tree.

        Args:
            subroot (str, optional): Count members of clades descended from this node. Defaults to the root.

        Returns:
            dict[str,int]: Dictionary containing clade counts.
        """
        cdef vector[string] leaves = self.t.get_leaves_ids(subroot.encode("UTF-8"))
        cdef size_t i,j,k
        cdef vector[bte.Node*] ancestors
        cdef vector[string] canns
        cdef string ann
        clade_counts = {}
        for i in range(leaves.size()):
            ancestors = self.t.rsearch(leaves[i], True)
            for j in range(ancestors.size()):
                canns = ancestors[j].clade_annotations
                for k in range(canns.size()):
                    ann = canns[k]
                    if ann.size() > 0:
                        annstr = ann.decode("UTF-8")
                        if annstr not in clade_counts:
                            clade_counts[annstr] = 0
                        clade_counts[annstr] += 1
        return clade_counts

    @_check_newick_only    
    @_timer
    def count_haplotypes(self) -> dict[tuple,int]:
        """Count unique haplotypes from the tree.

        Returns:
            dict[tuple,int]: haplotype counts.
        """        
        cdef map[cset[bte.Mutation],size_t] hmap = self.count_haplotypes_c()
        cdef map[cset[bte.Mutation],size_t].iterator it = hmap.begin()
        cdef cset[bte.Mutation].iterator it2
        pymap = {}
        while it != hmap.end():
            mset = set()
            it2 = dereference(it).first.begin()
            while it2 != dereference(it).first.end():
                mset.add(dereference(it2).get_string())
                postincrement(it2)
            pymap[tuple(mset)] = dereference(it).second
            postincrement(it)
        return pymap

    @_check_newick_only
    def compute_nucleotide_diversity(self) -> float:
        """Function which computes the nucleotide diversity of the tree. This is defined as the mean number of pairwise differences in nucleotides between any two leaves
        of the tree. Computes an unbiased estimator which multiplies the final mean by the total number of sequences divided by the total number of sequences minus one.

        Raises:
            Exception: Can't be computed on an empty tree.

        Returns:
            float: The estimated nucleotide diversity.
        """        

        '''
       '''
        #compute the set of mutations belonging to each sample in the tree, then compute pi from the resulting frequencies
        #since each sample is individually rsearched, this implementation is less efficient than an informed traversal, but still fast enough for most purposes.
        cdef map[cset[bte.Mutation],size_t] divtrack = self.count_haplotypes_c()
        cdef size_t total_seq = self.t.get_leaves_ids(self.t.root.identifier).size()
        if total_seq == 0:
            raise Exception("No sequences found in tree")
        cdef double div = 0
        cdef map[cset[bte.Mutation],size_t].iterator it = divtrack.begin()
        cdef map[cset[bte.Mutation],size_t].iterator it2
        cdef double g1_freq, g2_freq
        cdef size_t pair_diff
        while it != divtrack.end():
            g1_freq = dereference(it).second / total_seq
            assert g1_freq > 0
            it2 = it
            postincrement(it2)
            while it2 != divtrack.end():
                g2_freq = dereference(it2).second / total_seq
                assert g2_freq > 0
                pair_diff = self.count_differences(dereference(it).first,dereference(it2).first)
                div = div + (g1_freq * g2_freq * pair_diff)
                postincrement(it2)
            postincrement(it)
        #multiply the final result to guarantee an unbiased estimator (see wikipedia entry)
        return div * (total_seq / (total_seq - 1))

    def simple_parsimony(self, leaf_assignments: dict[str,str]) -> dict[str,str]:
        """This function is an implementation of the small parsimony problem (Fitch algorithm) for a single set of states.
        It takes as input a dictionary mapping leaf names to character states and returns a dictionary mapping both leaf and internal node names to inferred character states.

        Args:
            leaf_assignments (dict[str,str]): Dictionary mapping leaf names to character states.

        Returns:
            dict[str,str]: Dictionary mapping node names to inferred character states.
        """        
        ##TODO: refactoring this to use more c types would increase efficiency notably. It's fairly niche in application, though.
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
        cdef size_t i
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

    def ladderize(self) -> None:
        '''
        Sort the branches of the tree according to the size of each partition.
        '''
        self.t.rotate_for_consistency()

    @_check_newick_only
    def reverse_strand(self, genome_size: int = 29903) -> None:
        '''
        Inverts the tree representation of mutations such that all mutations are with respect to the reverse strand of the reference.
        All bases are complemented and indeces are reversed. The tree structure itself and parsimony scores are unaffected.

        Args:
            genome_size (int): The size of the genome. Assumes SARS-CoV-2 if not specified.
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
        
    def list_clades(self) -> set[str]:
        """Return a set of all valid clade annotations in the tree that can be used with get_clade and other functions.

        Returns:
            set[str]: Set of all valid clade annotations.
        """        
        clades = set()
        cdef vector[Node*] nodes = self.t.depth_first_expansion(self.t.root)
        cdef vector[string] anns
        for i in range(nodes.size()):
            anns = nodes[i].clade_annotations
            for j in range(anns.size()):
                if anns[j].size() > 0:
                    clades.add(anns[j].decode("UTF-8"))
        return clades
