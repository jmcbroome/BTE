from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from stringstream cimport stringstream

cdef extern from "mutation_annotated_tree.hpp" namespace "Mutation_Annotated_Tree":
    struct Missing_Sample:
        pass
    struct Mutation:
        string chrom
        int position
        int ref_nuc
        int par_nuc
        int mut_nuc
        bool is_missing
        Mutation copy()
        bool is_masked()
        string get_string()
    cppclass Node:
        size_t level
        float branch_length
        string identifier
        vector[string] clade_annotations
        Node* parent
        vector[Node*] children
        vector[Mutation] mutations

        bool is_leaf()
        bool is_root()

        Node() except + 
        Node(string id, float l)
        Node(string id, Node p, float l)

        void add_mutation(Mutation mut)
        void clear_mutations()
        void clear_annotations()
    cppclass Tree:
        Tree() except +
        Node* root
        size_t curr_internal_node
        size_t get_max_level() const
        size_t get_num_annotations() const
        void rename_node(string old_nid, string new_nid)
        vector[Node*] get_leaves(string nid)
        vector[string] get_leaves_ids(string nid)
        size_t get_num_leaves(Node* node)
        Node* create_node (string identifier, float branch_length, size_t num_annotations)
        Node* create_node (string identifier, Node* par, float branch_length)
        Node* create_node (string identifier, string parent_id, float branch_length)
        Node* get_node (string identifier) const
        bool is_ancestor (string anc_id, string nid) const
        vector[Node*] rsearch (const string nid, bool include_self) const
        string get_clade_assignment (const Node* n, int clade_id, bool include_self) const
        void remove_node (string nid, bool move_level)
        void remove_single_child_nodes()
        void move_node (string source, string destination, bool move_level)
        vector[Node*] breadth_first_expansion(string nid)
        vector[Node*] depth_first_expansion(Node* node) const
        size_t get_parsimony_score()

        void condense_leaves(vector[string])
        void uncondense_leaves()
        void collapse_tree()
        void rotate_for_display(bool reverse)
        void rotate_for_consistency()

    string get_newick_string(const Tree T, bool b1, bool b2, bool b3, bool b4)
    string get_newick_string(const Tree T, Node* node, bool b1, bool b2, bool b3, bool b4)
    void write_newick_string (stringstream ss, const Tree T, Node* node, bool b1, bool b2, bool b3, bool b4)
    Tree create_tree_from_newick (string filename)
    Tree create_tree_from_newick_string (string newick_string)
    void string_split(string s, char delim, vector[string] words)
    void string_split(string s, vector[string] words)
    Mutation* mutation_from_string(const string mut_string)

    Tree load_mutation_annotated_tree (string filename)
    void save_mutation_annotated_tree (Tree tree, string filename)

    Tree get_tree_copy(const Tree tree, const string identifier)

    Node* LCA (const Tree tree, const string node_id1, const string node_id2)
    Tree get_subtree (const Tree tree, const vector[string] samples, bool keep_clade_annotations)
    void get_random_single_subtree (Tree* T, vector[string] samples, string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len)
    void get_random_sample_subtrees (Tree* T, vector[string] samples, string outdir, size_t subtree_size, size_t tree_idx, bool use_tree_idx, bool retain_original_branch_len)
    void get_sample_mutation_paths (Tree* T, vector[string] samples, string mutation_paths_filename)
    void clear_tree(Tree tree)

    void read_vcf (Tree* T, string vcf_filename, vector[Missing_Sample] missing_samples, bool create_new_mat)
cdef extern from "mutation_annotated_tree.cpp":
    pass
cdef extern from "usher_graph.hpp":
    pass
cdef extern from "usher_mapper.cpp":
    pass
cdef extern from "parsimony.pb.h":
    pass
cdef extern from "parsimony.pb.cc":
    pass
cdef extern from "matUtils/common.hpp":
    pass
cdef extern from "matUtils/select.cpp":
    vector[string] get_clade_samples(Tree* T, string clade_name)
cdef extern from "matUtils/filter.cpp":
    Tree filter_master(Tree T, vector[string] samples, bool prune, bool keep_clade_annotations)
