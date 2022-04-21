/*this file exists to allow cython to correctly use the timer object extern and provide other C++ objects and functions as needed*/
#include "usher/src/usher_graph.hpp"
#include "usher/src/matUtils/translate.cpp"
Timer timer;

//additionally, define a version of the translation function which yields a vector of pairs instead of writing directly to a file
//this code is used in the cython translation.
std::vector<std::pair<std::string,std::string>> do_translation(MAT::Tree *T, std::string gtf_filename, std::string fasta_filename) {
    std::ifstream fasta_file(fasta_filename);
    if (!fasta_file) {
        fprintf(stderr, "ERROR: Could not open the fasta file: %s!\n", fasta_filename.c_str());
        exit(1);
    }
    std::ifstream gtf_file(gtf_filename);
    if (!gtf_file) {
        fprintf(stderr, "ERROR: Could not open the gtf file: %s!\n", gtf_filename.c_str());
        exit(1);
    }
    std::string reference = build_reference(fasta_file);

    // This maps each position in the reference to a vector of codons.
    // Some positions may be associated with multiple codons (frame shifts).
    // The Codons in the map are updated as the tree is traversed
    std::vector<std::pair<std::string,std::string>> results;
    std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> codon_map = build_codon_map(gtf_file, reference);

    // Traverse the tree in depth-first order. As we descend the tree, mutations at
    // each node are applied to the respective codon(s) in codon_map.
    auto dfs = T->depth_first_expansion();
    MAT::Node *last_visited = nullptr;
    for (auto &node: dfs) {
        std::string mutation_result = "";
        if (last_visited != node->parent) {
            // Jumping across a branch, so we need to revert codon mutations up to
            // the LCA of this node and the last visited node
            MAT::Node *last_common_ancestor = MAT::LCA(*T, node->identifier, last_visited->identifier);
            MAT::Node *trace_to_lca = last_visited;
            while (trace_to_lca != last_common_ancestor) {
                undo_mutations(trace_to_lca->mutations, codon_map);
                trace_to_lca = trace_to_lca->parent;
            }
        } // If we are visiting a child, we can continue mutating
        mutation_result = do_mutations(node->mutations, codon_map, false);
        if (mutation_result != "") {
            std::pair<std::string,std::string> result = std::make_pair(node->identifier, mutation_result);
            results.push_back(result);
        }
        last_visited = node;
    }
    return results;
}