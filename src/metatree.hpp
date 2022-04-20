#include "usher/src/mutation_annotated_tree.hpp"

struct AAMutation {
    std::string chrom;
    std::string gene;
    int codon;
    char ref_aa;
    char par_aa;
    char mut_aa;
    std::vector<Mutation> causal_mutations;
    inline bool operator< (const AAMutation& m) const {
        return ((*this).position < m.position);
    }
    AAMutation () {
        chrom = "";
        gene = "";
        is_missing = false;
    }
    inline std::string get_string() const {
        return gene + ":" + ref_aa + std::to_string(codon) + mut_aa;
    }
}

class MetaNode: public MAT::Node {
    """
    The MetaNode class is a derived MAT::Node which simply adds a vector to store amino acid mutations
    and an unordered map for arbitrary node-level metadata in string format.

    Calling the translation routine from a MATree will cast all nodes to MetaNode. 
    """
    public:
        std::vector<AAMutation> aa_mutations;
        std::unordered_map<string,string> metadata;
        void add_aa_mutation(string aachange);
}

class MetaTree: public MAT::Tree {
    public:
        std::unordered_map<std::string> 

}