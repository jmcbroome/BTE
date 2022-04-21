/*this file exists to allow cython to correctly use the timer object extern and provide other C++ objects and functions as needed*/
#include "usher/src/usher_graph.hpp"
Timer timer;

/*translation functions copied from translate.cpp follow.
It's seemingly impossible to include two protoc compilers in a single project using cython compilers without redefinition errors,
so BTE is incompatible with Taxodium and anything Taxodium related. As several Taxodium functions exist in 
translate.cpp, we've replicated only the ones we need- that aren't taxodium related- here. 

TODO: find some way to include multiple compiled protoc definitions in a single cython project.
Alternatively, refactor matUtils to separate translation and taxodium function into separate source code files.*/

static std::unordered_map<std::string, char> translation_map = {
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'}, {"GCN", 'A'},
    {"TGT", 'C'}, {"TGC", 'C'}, {"TGY", 'C'},
    {"GAT", 'D'}, {"GAC", 'D'}, {"GAY", 'D'},
    {"GAA", 'E'}, {"GAG", 'E'}, {"GAR", 'E'},
    {"TTT", 'F'}, {"TTC", 'F'}, {"TTY", 'F'},
    {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}, {"GGN", 'G'},
    {"CAT", 'H'}, {"CAC", 'H'}, {"CAY", 'H'},
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATH", 'I'},
    {"AAA", 'K'}, {"AAG", 'K'}, {"AAR", 'K'},
    {"TTA", 'L'}, {"TTG", 'L'}, {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'}, {"YTR", 'L'}, {"CTN", 'L'},
    {"ATG", 'M'},
    {"AAT", 'N'}, {"AAC", 'N'}, {"AAY", 'N'},
    {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'}, {"CCN", 'P'},
    {"CAA", 'Q'}, {"CAG", 'Q'}, {"CAR", 'Q'},
    {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'}, {"AGA", 'R'}, {"AGG", 'R'}, {"CGN", 'R'}, {"MGR", 'R'},
    {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'}, {"AGT", 'S'}, {"AGC", 'S'}, {"TCN", 'S'}, {"AGY", 'S'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'}, {"ACN", 'T'},
    {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'}, {"GTN", 'V'},
    {"TGG", 'W'},
    {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAY", 'Y'},
    {"TAG", '*'}, {"TAA", '*'}, {"TGA", '*'}
};

static std::unordered_map<char, char> complement_map = {
    {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'},
    {'M', 'K'}, {'R', 'Y'}, {'W', 'W'}, {'S', 'S'},
    {'Y', 'R'}, {'K', 'M'}, {'V', 'B'}, {'H', 'D'},
    {'D', 'H'}, {'B', 'V'}, {'N', 'N'}
};

struct Codon {
    std::string orf_name;
    std::string nucleotides;
    int codon_number;
    int start_position;
    char protein;

    // Translate codon to amino acid, allowing for ambiguous codons
    inline char translate_codon(std::string nt) {
        auto it = translation_map.find(nt);
        if (it == translation_map.end()) {
            return 'X'; // ambiguous, couldn't resolve aa
        } else {
            return it->second;
        }
    }

    inline void mutate(int nuc_pos, char mutated_nuc) {
        // The nt to mutate is the difference between the
        // genomic coordinate of the mutated nt and the
        // starting coordinate of the codon
        nucleotides[abs(nuc_pos-start_position)] = mutated_nuc;
        protein = translate_codon(nucleotides);
    }

    Codon (std::string _orf_name, int _codon_number, int _start_position, char nt[3]) {
        orf_name = _orf_name;
        start_position = _start_position;
        codon_number = _codon_number;
        nucleotides = "";
        nucleotides += nt[0];
        nucleotides += nt[1];
        nucleotides += nt[2];
        protein = translate_codon(nt);
    }

    inline std::string get_string() const {
        return std::to_string(start_position) + ':'
               + nucleotides[0]
               + nucleotides[1]
               + nucleotides[2]
               + '=' + protein;
    }
};

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

std::string build_reference(std::ifstream &fasta_file) {
    std::string reference_output = "";
    std::string fasta_line;
    size_t line_length;
    while(std::getline(fasta_file, fasta_line)) {
        if (fasta_line[0] == '>' or fasta_line[0] == '\n') {
            continue;
        } else {
            for (auto & c: fasta_line) c = (char)toupper(c);
            line_length = fasta_line.length();
            if (fasta_line[line_length-1] == '\r') {
                fasta_line.erase(line_length-1);
            }
            reference_output += fasta_line;
        }
    }
    return reference_output;
}

char complement(char nt) {
    auto it = complement_map.find(nt);
    if (it == complement_map.end()) {
        return 'N'; // ambiguous, couldn't resolve nt
    } else {
        return it->second;
    }
}


std::string do_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map, bool taxodium_format) {
    std::string prot_string = "";
    std::string nuc_string = "";
    std::string cchange_string = "";
    std::sort(mutations.begin(), mutations.end());
    std::unordered_map<std::string, std::set<MAT::Mutation>> codon_to_nt_map;
    std::unordered_map<std::string, std::string> codon_to_changestring_map;
    std::unordered_map<std::string, char> orig_proteins;
    std::vector<std::shared_ptr<Codon>> affected_codons;

    for (auto &m : mutations) {
        char mutated_nuc = MAT::get_nuc(m.mut_nuc);
        char par_nuc = MAT::get_nuc(m.par_nuc);
        int pos = m.position - 1;
        auto codon_map_it = codon_map.find(pos);
        if (codon_map_it == codon_map.end()) {
            continue; // Not a coding mutation
        } else {
            // Mutate each codon associated with this position
            for (auto codon_ptr : codon_map_it->second) {
                std::string codon_id = codon_ptr->orf_name + ':' + std::to_string(codon_ptr->codon_number+1);
                //first, update the codon to match the parent state instead of the reference state as part of the codon output
                codon_ptr->mutate(pos, par_nuc);
                auto orig_it = orig_proteins.find(codon_id);
                if (orig_it == orig_proteins.end()) {
                    orig_proteins.insert({codon_id, codon_ptr->protein});
                }
                if (std::find(affected_codons.begin(), affected_codons.end(), codon_ptr) == affected_codons.end()) {
                    affected_codons.push_back(codon_ptr);
                }
                std::string original_codon = codon_ptr->nucleotides;
                //then update it again to match the mutated state
                codon_ptr->mutate(pos, mutated_nuc);
                // store a string representing the original and new codons in nucleotides
                // this may incorporate multiple nucleotide mutations, just accounting for the original and end states.
                std::string changestring = original_codon + ">" + codon_ptr->nucleotides;
                codon_to_changestring_map.insert({codon_id, changestring});
                // Build a map of codons and their associated nt mutations
                auto to_nt_it = codon_to_nt_map.find(codon_id);
                if (to_nt_it == codon_to_nt_map.end()) {
                    codon_to_nt_map.insert({codon_id, {m}});
                } else {
                    to_nt_it->second.insert(m);
                }
            }
        }
    }

    for (auto codon_ptr : affected_codons) {
        std::string codon_id = codon_ptr->orf_name + ':' + std::to_string(codon_ptr->codon_number+1);
        char orig_protein = orig_proteins.find(codon_id)->second;
        if (taxodium_format) {
            if (orig_protein == codon_ptr->protein) { // exclude synonymous mutations
                continue;
            }
            prot_string += split(codon_id, ':')[0] + ':' + orig_protein + '_' + split(codon_id, ':')[1] + '_' + codon_ptr->protein + ';';
        } else {
            prot_string += split(codon_id, ':')[0] + ':' + orig_protein + split(codon_id, ':')[1] + codon_ptr->protein + ';';
        }
        auto codon_it = codon_to_nt_map.find(codon_id);
        for (auto &m : codon_it->second) {
            nuc_string += m.get_string() + ",";
        }

        if (!nuc_string.empty() && nuc_string.back() == ',') {
            nuc_string.resize(nuc_string.length() - 1); // remove trailing ','
            nuc_string += ';';
        }
        std::string changestring = codon_to_changestring_map.find(codon_id)->second;
        cchange_string += changestring + ";";
    }

    if (!nuc_string.empty() && nuc_string.back() == ';') {
        nuc_string.resize(nuc_string.length() - 1); // remove trailing ';'
    }
    if (!prot_string.empty() && prot_string.back() == ';') {
        prot_string.resize(prot_string.length() - 1); //remove trailing ';'
    }
    if (!cchange_string.empty() && cchange_string.back() == ';') {
        cchange_string.resize(cchange_string.length() - 1); //remove trailing ';'
    }
    if (nuc_string.empty() || prot_string.empty() || cchange_string.empty()) {
        return "";
    } else if(taxodium_format) { // format string for taxodium pb
        return prot_string;
    } else { // format string for TSV output
        return prot_string + '\t' + nuc_string + '\t' + cchange_string;
    }
}

void undo_mutations(std::vector<MAT::Mutation> &mutations, std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> &codon_map) {
    for (auto &m: mutations) {
        char parent_nuc = MAT::get_nuc(m.par_nuc);
        int pos = m.position - 1;
        auto it = codon_map.find(pos);
        if (it == codon_map.end()) {
            continue;
            // Not a coding mutation
        } else {
            // Revert the mutation by mutating to the parent nucleotide
            for (auto codon_ptr : it->second) {
                codon_ptr->mutate(pos, parent_nuc);
            }
        }
    }
}

// Maps a genomic coordinate to a list of codons it is part of
std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> build_codon_map(std::ifstream &gtf_file, std::string reference) {
    std::unordered_map<int, std::vector<std::shared_ptr<Codon>>> codon_map;
    std::string gtf_line;
    std::vector<std::string> gtf_lines;
    std::vector<std::string> done;

    while (std::getline(gtf_file, gtf_line)) {
        gtf_lines.push_back(gtf_line);
    }
    int curr_line = -1;
    for (std::string line_outer : gtf_lines) {
        curr_line += 1;

        if (line_outer[0] == '#' || line_outer[0] == '\n') {
            continue;
        }
        std::vector<std::string> split_line_outer = split(line_outer, '\t');
        if(split_line_outer.size() <= 1) {
            continue;
        }
        if (split_line_outer[8].substr(0, 7) != "gene_id") {
            fprintf(stderr, "ERROR: GTF file formatted incorrectly. Please see the UShER wiki for details.\n");
            exit(1);
        }
        std::string feature_outer = split_line_outer[2];
        std::string gene_outer = split(split(split_line_outer[8], '\"')[1], '\"')[0];
        char strand_outer = split_line_outer[6][0];

        if (feature_outer == "CDS") {
            bool found = (std::find(done.begin(), done.end(), gene_outer) != done.end());
            if (found) {
                continue;
            } else {
                done.push_back(gene_outer);
            }
            // There may be multiple CDS features per gene. First build codons for the first CDS
            int first_cds_start = std::stoi(split_line_outer[3]); // expect the GTF is ordered by start position
            int first_cds_stop = std::stoi(split_line_outer[4]);
            int codon_counter = 0; // the number of codons we have added so far
            if (strand_outer == '+') {
                for (int pos = first_cds_start - 1; pos < first_cds_stop; pos += 3) {

                    char nt[3] = {
                        reference[pos],
                        reference[pos+1],
                        reference[pos+2]
                    };

                    // Coordinates are 0-based at this point
                    std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt));
                    codon_counter += 1;

                    // The current pos and the next positions
                    // are associated with this codon
                    auto it = codon_map.find(pos);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos, {c}});
                    } else {
                        (it->second).push_back(c);
                    }

                    it = codon_map.find(pos+1);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos+1, {c}});
                    } else {
                        (it->second).push_back(c);
                    }

                    it = codon_map.find(pos+2);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos+2, {c}});
                    } else {
                        (it->second).push_back(c);
                    }
                }
            } else {
                for (int pos = first_cds_stop - 1; pos > first_cds_start; pos -= 3) {

                    char nt[3] = {
                        complement(reference[pos]),
                        complement(reference[pos-1]),
                        complement(reference[pos-2])
                    };

                    // Coordinates are 0-based at this point
                    std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt));
                    codon_counter += 1;

                    // The current pos and the next positions
                    // are associated with this codon
                    auto it = codon_map.find(pos);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos, {c}});
                    } else {
                        (it->second).push_back(c);
                    }

                    it = codon_map.find(pos-1);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos-1, {c}});
                    } else {
                        (it->second).push_back(c);
                    }

                    it = codon_map.find(pos-2);
                    if (it == codon_map.end()) {
                        codon_map.insert({pos-2, {c}});
                    } else {
                        (it->second).push_back(c);
                    }
                }
            }
            for (std::string line_inner : gtf_lines) { // find the rest of the CDS features, assuming they are in position order

                if (line_inner[0] == '#' || line_inner[0] == '\n') {
                    continue;
                }
                std::vector<std::string> split_line_inner = split(line_inner, '\t');
                std::string feature_inner = split_line_inner[2];
                std::string gene_inner = split(split(split_line_inner[8], '\"')[1], '\"')[0];
                if (feature_inner == "CDS" && gene_outer == gene_inner) {
                    int inner_cds_start = std::stoi(split_line_inner[3]);
                    int inner_cds_stop = std::stoi(split_line_inner[4]);
                    char strand_inner = split_line_inner[6][0];
                    if (strand_inner == '+') {
                        if (inner_cds_start != first_cds_start || strand_outer != strand_inner) {
                            for (int pos = inner_cds_start - 1; pos < inner_cds_stop; pos += 3) {
                                char nt[3] = {
                                    reference[pos],
                                    reference[pos+1],
                                    reference[pos+2]
                                };
                                std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt));
                                codon_counter += 1;

                                auto it = codon_map.find(pos);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos+1);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos+1, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos+2);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos+2, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }
                            }
                        }
                    } else {
                        if (inner_cds_start != first_cds_start || strand_outer != strand_inner) {
                            for (int pos = inner_cds_stop - 1; pos > inner_cds_start; pos -= 3) {
                                char nt[3] = {
                                    complement(reference[pos]),
                                    complement(reference[pos-1]),
                                    complement(reference[pos-2])
                                };
                                std::shared_ptr<Codon> c(new Codon(gene_outer, codon_counter, pos, nt));
                                codon_counter += 1;

                                auto it = codon_map.find(pos);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos-1);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos-1, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }

                                it = codon_map.find(pos-2);
                                if (it == codon_map.end()) {
                                    codon_map.insert({pos-2, {c}});
                                } else {
                                    (it->second).push_back(c);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return codon_map;
}

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
