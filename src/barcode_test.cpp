#include "include/rad/rad_headers.h"

std::unordered_set<int64_seq> load_simple_wl(const std::string& filepath) {
    std::unordered_set<int64_seq> barcodes;
    std::string line;
    
    // Check if file is gzipped
    if (filepath.size() >= 3 && filepath.substr(filepath.size() - 3) == ".gz") {
        // Handle gzipped file
        gzFile file = gzopen(filepath.c_str(), "rb");
        if (!file) {
            std::cerr << "Error: Could not open gzipped file " << filepath << std::endl;
            return barcodes;
        }
        
        char buffer[1024];
        std::string current_line;
        
        while (gzgets(file, buffer, sizeof(buffer))) {
            current_line = buffer;
            
            // Remove newline if present
            if (!current_line.empty() && current_line.back() == '\n') {
                current_line.pop_back();
            }
            if (!current_line.empty() && current_line.back() == '\r') {
                current_line.pop_back();
            }
            
            // Trim whitespace and skip empty lines
            current_line.erase(0, current_line.find_first_not_of(" \t"));
            current_line.erase(current_line.find_last_not_of(" \t") + 1);
            
            if (!current_line.empty()) {
                barcodes.emplace(current_line);
            }
        }
        
        gzclose(file);
    } else {
        // Handle regular text file
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filepath << std::endl;
            return barcodes;
        }
        
        while (std::getline(file, line)) {
            // Trim whitespace and skip empty lines
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            if (!line.empty()) {
                barcodes.emplace(line);
            }
        }
    }
    
    return barcodes;
}

// long bit-parallel algo for two strings of different lengths (up to 64)
int bit_ld_variable(int64_t pattern, int64_t text, int pattern_len, int text_len, int max_dist) {
    // Validate inputs
    if (pattern_len <= 0 || text_len <= 0 || pattern_len > 32 || text_len > 32) {
        return -1;
    }
    
    // Build pattern equality vectors for the actual pattern length
    int64_t Peq[4] = {0, 0, 0, 0};
    for (int i = 0; i < pattern_len; i++) {
        int nuc = (pattern >> (2 * i)) & 3;
        Peq[nuc] |= (1LL << i);
    }
    
    // Initialize with pattern length
    int64_t pattern_mask = (1LL << pattern_len) - 1; // Mask for valid pattern bits
    int64_t Pv = pattern_mask; // All 1s for pattern length
    int64_t Mv = 0; // All 0s
    int score = pattern_len; // Start with pattern length (all insertions)
    
    // Process each character in the text
    for (int j = 0; j < text_len; j++) {
        int text_nuc = (text >> (2 * j)) & 3; // Extract nucleotide from text
        int64_t Eq = Peq[text_nuc] & pattern_mask; // Get equality vector, masked to pattern length
        
        // Myers core computation
        int64_t Xv = Eq | Mv;
        int64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        int64_t Ph = Mv | ~(Xh | Pv);
        int64_t Mh = Pv & Xh;
        
        // Apply pattern mask to keep only valid bits
        Ph &= pattern_mask;
        Mh &= pattern_mask;
        
        // Update score based on the last position of the pattern
        if (Ph & (1LL << (pattern_len - 1))) score++;
        if (Mh & (1LL << (pattern_len - 1))) score--;
        
        // Early exit optimization
        // Best possible score is current score minus remaining text characters
        // (assuming all remaining are matches, which would decrease score by 1 each)
        int remaining_text = text_len - j - 1;
        int best_possible = score - remaining_text;
        if (best_possible > max_dist) {
            return -1;
        }
        
        // Update for next iteration
        Ph <<= 1;
        Pv = ((Mh << 1) | ~(Xv | Ph)) & pattern_mask;
        Mv = Ph & Xv;
    }
    
    return score;
}

// calculating levenshtein distance between one query and a set of multiple strings of different lengths (naively, with no wildcard-based spacing)
std::map<int, std::unordered_set<int64_seq>> int64_lvdist_long(const int64_seq &query, const std::unordered_set<int64_seq> &targets, int max_dist = 4) {
    std::map<int, std::unordered_set<int64_seq>> results;
    if (query.bits.empty()){
        return results;
    }
    for (const auto &target : targets) {
        if (target.bits.empty()){
            continue;
        }
        int dist = mutation_tools::bit_partial_match(query.bits[0], target.bits[0], query.length, target.length, max_dist);
        if(dist >= 0){
            results[dist].insert(target);
        }
    }
    return results;
}

//lvdist for long sequences, different length, whole whitelist
std::map<int, std::unordered_set<int64_seq>> int64_lvdist_partial(const int64_seq &query, const std::unordered_set<int64_seq> &targets, int max_dist = 4) {
    std::map<int, std::unordered_set<int64_seq>> results;
    if (query.bits.empty()){
        return results;
    }
    for (const auto &target : targets) {
        if (target.bits.empty()){
            continue;
        }
        
        int dist;
        if (query.length <= target.length) {
            // Query is shorter or equal, search query in target
            dist = mutation_tools::bit_partial_match(query.bits[0], target.bits[0], query.length, target.length, max_dist);
        } else {
            // Target is shorter, search target in query
            dist = mutation_tools::bit_partial_match(target.bits[0], query.bits[0], target.length, query.length, max_dist);
        }
        
        if(dist >= 0){
            results[dist].insert(target);
        }
    }
    return results;
}

void print_lvdist_results(const std::map<int, std::unordered_set<int64_seq>>& results, const std::string& query_seq = "") {
    if (results.empty()) {
        std::cout << "No matches found.\n";
        return;
    }
    
    if (!query_seq.empty()) {
        std::cout << "Query: " << query_seq << "\n";
    }
    std::cout << "Partial match results:\n";
    
    int total_matches = 0;
    for (const auto& [distance, sequences] : results) {
        std::cout << "Distance " << distance << " (" << sequences.size() << " matches):\n";
        for (const auto& seq : sequences) {
            if(distance <= 2){
                std::cout << "  " << seq.bits_to_sequence() << "\n";
            }
        }
        total_matches += sequences.size();
    }
    
    std::cout << "Total matches: " << total_matches << "\n\n";
}

// reference edit distance testing
int reference_edit_distance(const std::string& s1, const std::string& s2) {
    int m = s1.length();
    int n = s2.length();
    
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));
    
    for (int i = 0; i <= m; i++) dp[i][0] = i;
    for (int j = 0; j <= n; j++) dp[0][j] = j;
    
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (s1[i-1] == s2[j-1]) {
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = 1 + std::min({dp[i-1][j], dp[i][j-1], dp[i-1][j-1]});
            }
        }
    }
    
    return dp[m][n];
}

// testing the debug function
void debug_sequence_encoding(const std::string& seq) {
    int64_seq encoded(seq);
    std::string decoded = encoded.bits_to_sequence();
    
    std::cout << "Original: " << seq << "\n";
    std::cout << "Decoded:  " << decoded << "\n";
    std::cout << "Length:   " << encoded.length << "\n";
    std::cout << "Bits:     " << std::bitset<64>(encoded.bits[0]) << "\n";
    std::cout << "Bits hex: 0x" << std::hex << encoded.bits[0] << std::dec << "\n";
    std::cout << "Match:    " << (seq == decoded ? "YES" : "NO") << "\n\n";
}

void test_bit_ld(const std::string& seq1, const std::string& seq2) {
    int64_seq s1(seq1);
    int64_seq s2(seq2);
    
    std::cout << "=== Direct bit_ld test ===\n";
    std::cout << "Seq1: " << seq1 << "\n";
    std::cout << "Seq2: " << seq2 << "\n";
    
    // Check if lengths match
    if (s1.length != s2.length) {
        std::cout << "ERROR: Length mismatch! " << s1.length << " vs " << s2.length << "\n\n";
        return;
    }
    
    int ref_dist = reference_edit_distance(seq1, seq2);
    int bit_dist = mutation_tools::bit_ld(s1.bits[0], s2.bits[0], s1.length, 10);
    int wrapper_dist = mutation_tools::int64_lvdist(s1, s2, 10);
    
    std::cout << "Reference distance: " << ref_dist << "\n";
    std::cout << "bit_ld distance:    " << bit_dist << "\n";
    std::cout << "Wrapper distance:   " << wrapper_dist << "\n";
    std::cout << "Match: " << (ref_dist == bit_dist && ref_dist == wrapper_dist ? "YES" : "NO") << "\n\n";
}

void test_bit_ld_variable(const std::string& seq1, const std::string& seq2) {
    int64_seq s1(seq1);
    int64_seq s2(seq2);
    
    std::cout << "=== Direct bit_ld test ===\n";
    std::cout << "Seq1: " << seq1 << " (length: " << s1.length << ")\n";
    std::cout << "Seq2: " << seq2 << " (length: " << s2.length << ")\n";
    
    // Calculate reference distance (works for any lengths)
    int ref_dist = reference_edit_distance(seq1, seq2);
    std::cout << "Reference distance: " << ref_dist << "\n";
    
    // Test the bit-parallel algorithm for variable lengths
    int dist = bit_ld_variable(s1.bits[0], s2.bits[0], s1.length, s2.length, 10);
    std::cout << "bit_ld_variable:    " << dist << "\n";
    
    // Only test the original bit_ld if lengths are equal
    if (s1.length == s2.length) {
        int bit_dist = mutation_tools::bit_ld(s1.bits[0], s2.bits[0], s1.length, 10);
        int wrapper_dist = mutation_tools::int64_lvdist(s1, s2, 10);
        
        std::cout << "bit_ld distance:    " << bit_dist << "\n";
        std::cout << "Wrapper distance:   " << wrapper_dist << "\n";
        
        bool all_match = (ref_dist == bit_dist && ref_dist == wrapper_dist && ref_dist == dist);
        std::cout << "All algorithms match: " << (all_match ? "YES" : "NO") << "\n";
        
        if (!all_match) {
            std::cout << "MISMATCH DETAILS:\n";
            std::cout << "  Reference vs bit_ld: " << (ref_dist == bit_dist ? "OK" : "FAIL") << "\n";
            std::cout << "  Reference vs wrapper: " << (ref_dist == wrapper_dist ? "OK" : "FAIL") << "\n";
            std::cout << "  Reference vs bit_ld_variable: " << (ref_dist == dist ? "OK" : "FAIL") << "\n";
        }
    } else {
        // Different lengths - only compare reference with bit_ld_variable
        bool match = (ref_dist == dist);
        std::cout << "Variable-length match: " << (match ? "YES" : "NO") << "\n";
        
        if (!match) {
            std::cout << "MISMATCH: Reference=" << ref_dist << ", bit_ld_variable=" << dist << "\n";
        }
        
        std::cout << "Note: Original bit_ld and wrapper only work with equal-length sequences\n";
    }
    
    std::cout << "\n";
}

void test_bit_ld_partial(const std::string & full_seq, const std::string & path) {
    int64_seq full(full_seq);
    std::unordered_set<int64_seq> targets = load_simple_wl(path);
    std::cout << "=== Partial bit_ld test ===\n";
    std::cout << "Full sequence: " << full_seq << " (length: " << full.length << ")\n";
    std::cout << "Targets loaded: " << targets.size() << "\n";  

    // Run the partial bit_ld distance calculation
    std::cout << "Calculating partial distances...\n";
    auto results = int64_lvdist_partial(full, targets);
    print_lvdist_results(results, full_seq);
}

void test_bit_ld_partial_timing(const std::string & full_seq, const std::string & path) {
    auto start_total = std::chrono::high_resolution_clock::now();
    
    int64_seq full(full_seq);
    
    auto start_load = std::chrono::high_resolution_clock::now();
    std::unordered_set<int64_seq> targets = load_simple_wl(path);
    auto end_load = std::chrono::high_resolution_clock::now();
    auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load);
    
    std::cout << "=== Partial bit_ld test ===\n";
    std::cout << "Full sequence: " << full_seq << " (length: " << full.length << ")\n";
    std::cout << "Targets loaded: " << targets.size() << " (took " << load_time.count() << "ms)\n";  
    
    // Run the partial bit_ld distance calculation
    std::cout << "Calculating partial distances...\n";
    auto start_calc = std::chrono::high_resolution_clock::now();
    auto results = int64_lvdist_partial(full, targets, 1);
    auto end_calc = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_calc - start_calc);
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    
    std::cout << "Calculation completed in " << calc_time.count() << "ms\n";
    std::cout << "Total time: " << total_time.count() << "ms\n\n";
    
    print_lvdist_results(results, full_seq);
}

// Helper function to run comprehensive tests
void run_comprehensive_tests() {
    std::cout << "=== ACTUAL DATA TESTS ===\n";
    test_bit_ld_variable("AAGCCAAGAAGATCTG", "CGATCTAAGCCAAGAAGATCTGGTGTCG");   // identical match, subsequence
    test_bit_ld_variable("TAAGCCAAGAAGATCT", "CGATCTAAGCCAAGAAGATCTGGTGTCG");   // identical match, subsequence, part 2
    test_bit_ld_variable("AAACCCAAGAAGATCT", "CGATCTAAGCCAAGAAGATCTGGTGTCG");   // original barcode
    test_bit_ld_variable("CGATCTAAGCCAAGAAGATCTGGTGTCG","AAACCCAAGAAGATCT");   // original barcode
    test_bit_ld_partial_timing("CGATCTAAGCCAAGAAGATCTGGTGTCG", "/Users/cmv/Desktop/rad_paper/cellranger_whitelists/3M-february-2018-3v3.txt.gz"); // partial match test with file
    test_bit_ld_variable("AAACCCACAGATCGTT","AAACCCAAAGATCTTTCGTGTATTTT");
    test_bit_ld_variable("AAACCCAAGAAGATCT","AAACCCAAAGATCTTTCGTGTATTTT");

}

struct BitTrie {
    struct TrieNode {
        TrieNode* children[4];  // A=0, T=1, C=2, G=3
        bool is_end;
        int64_seq sequence;     // Store the complete sequence at leaf nodes
        
        TrieNode() : is_end(false) {
            for (int i = 0; i < 4; i++) {
                children[i] = nullptr;
            }
        }
        
        ~TrieNode() {
            for (int i = 0; i < 4; i++) {
                delete children[i];
            }
        }
    };
    
    TrieNode* root;
    
    BitTrie() {
        root = new TrieNode();
    }
    
    ~BitTrie() {
        delete root;
    }
    
    // Insert a sequence using its bit encoding
    void insert(const int64_seq& seq) {
        if (seq.bits.empty() || seq.length > 32) return;
        
        TrieNode* current = root;
        int64_t bits = seq.bits[0];
        
        // Traverse from least significant bits (position 0) to most significant
        for (int pos = 0; pos < seq.length; pos++) {
            int nucleotide = (bits >> (2 * pos)) & 3;  // Extract 2 bits
            
            if (current->children[nucleotide] == nullptr) {
                current->children[nucleotide] = new TrieNode();
            }
            current = current->children[nucleotide];
        }
        
        current->is_end = true;
        current->sequence = seq;
    }
    
    // Exact search for a sequence
    bool search(const int64_seq& seq) {
        TrieNode* node = find_node(seq);
        return node != nullptr && node->is_end;
    }
    
    // Find all sequences with exact prefix match
    std::vector<int64_seq> prefix_search(const int64_seq& prefix) {
        std::vector<int64_seq> results;
        TrieNode* prefix_node = find_node(prefix);
        
        if (prefix_node != nullptr) {
            collect_all_sequences(prefix_node, results);
        }
        
        return results;
    }
    
    // Approximate search allowing up to max_mismatches
    std::vector<std::pair<int64_seq, int>> approximate_search(const int64_seq& query, int max_mismatches) {
        std::vector<std::pair<int64_seq, int>> results;
        if (query.bits.empty() || query.length > 32) return results;
        
        approximate_search_helper(root, query.bits[0], query.length, 0, 0, max_mismatches, results);
        return results;
    }
    
private:
    // Helper to find a node for exact sequence
    TrieNode* find_node(const int64_seq& seq) {
        if (seq.bits.empty() || seq.length > 32) return nullptr;
        
        TrieNode* current = root;
        int64_t bits = seq.bits[0];
        
        for (int pos = 0; pos < seq.length; pos++) {
            int nucleotide = (bits >> (2 * pos)) & 3;
            
            if (current->children[nucleotide] == nullptr) {
                return nullptr;
            }
            current = current->children[nucleotide];
        }
        
        return current;
    }
    
    // Collect all sequences in subtree
    void collect_all_sequences(TrieNode* node, std::vector<int64_seq>& results) {
        if (node->is_end) {
            results.push_back(node->sequence);
        }
        
        for (int i = 0; i < 4; i++) {
            if (node->children[i] != nullptr) {
                collect_all_sequences(node->children[i], results);
            }
        }
    }
    
    // Recursive approximate search with mismatch counting
    void approximate_search_helper(TrieNode* node, int64_t query_bits, int query_len, 
                                 int pos, int mismatches, int max_mismatches,
                                 std::vector<std::pair<int64_seq, int>>& results) {
        
        // If we've reached the end of the query
        if (pos == query_len) {
            if (node->is_end) {
                results.emplace_back(node->sequence, mismatches);
            }
            return;
        }
        
        // Early termination if too many mismatches
        if (mismatches > max_mismatches) {
            return;
        }
        
        int query_nucleotide = (query_bits >> (2 * pos)) & 3;
        
        // Try all possible nucleotides at this position
        for (int nuc = 0; nuc < 4; nuc++) {
            if (node->children[nuc] != nullptr) {
                int new_mismatches = mismatches + (nuc != query_nucleotide ? 1 : 0);
                approximate_search_helper(node->children[nuc], query_bits, query_len,
                                        pos + 1, new_mismatches, max_mismatches, results);
            }
        }
    }
    
public:
    // Build trie from a set of sequences
    static BitTrie build_from_set(const std::unordered_set<int64_seq>& sequences) {
        BitTrie trie;
        for (const auto& seq : sequences) {
            trie.insert(seq);
        }
        return trie;
    }
    
    // Get statistics about the trie
    void print_stats() {
        int total_nodes = count_nodes(root);
        int leaf_nodes = count_leaves(root);
        std::cout << "Trie stats: " << total_nodes << " total nodes, " 
                  << leaf_nodes << " leaf nodes\n";
    }
    
private:
    int count_nodes(TrieNode* node) {
        if (node == nullptr) return 0;
        
        int count = 1;
        for (int i = 0; i < 4; i++) {
            count += count_nodes(node->children[i]);
        }
        return count;
    }
    
    int count_leaves(TrieNode* node) {
        if (node == nullptr) return 0;
        if (node->is_end) return 1;
        
        int count = 0;
        for (int i = 0; i < 4; i++) {
            count += count_leaves(node->children[i]);
        }
        return count;
    }
};

// Test function using entire whitelist
// Test function using entire whitelist with k-mer approach
void test_bit_trie(const std::string& query_seq, const std::string& whitelist_path, int max_mismatches = 2, int kmer_size = 16) {
    auto start_total = std::chrono::high_resolution_clock::now();
    
    // Load whitelist
    std::cout << "=== BitTrie Test ===\n";
    std::cout << "Query: " << query_seq << "\n";
    std::cout << "Loading whitelist from: " << whitelist_path << "\n";
    
    auto start_load = std::chrono::high_resolution_clock::now();
    auto whitelist = load_simple_wl(whitelist_path);
    auto end_load = std::chrono::high_resolution_clock::now();
    auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load);
    
    std::cout << "Loaded " << whitelist.size() << " sequences (" << load_time.count() << "ms)\n";
    
    // Build trie
    std::cout << "Building trie...\n";
    auto start_build = std::chrono::high_resolution_clock::now();
    BitTrie trie = BitTrie::build_from_set(whitelist);
    auto end_build = std::chrono::high_resolution_clock::now();
    auto build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_build - start_build);
    
    std::cout << "Trie built (" << build_time.count() << "ms)\n";
    trie.print_stats();
    
    // Test k-mer searches
    auto kmers = seq_utils::kmerize(query_seq, kmer_size);
    std::cout << "\nGenerated " << kmers.size() << " k-mers of size " << kmer_size << "\n";
    
    for (const auto& kmer : kmers) {
        std::cout << "Testing k-mer: " << kmer << "\n";
        int64_seq kmer_seq(kmer);
        
        // Exact search for this k-mer
        auto start_exact = std::chrono::high_resolution_clock::now();
        bool exact_found = trie.search(kmer_seq);
        auto end_exact = std::chrono::high_resolution_clock::now();
        auto exact_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_exact - start_exact);
        
        std::cout << "  Exact search: " << (exact_found ? "Found" : "Not found")
                  << " (" << exact_time.count() << "ms)\n";
        
        // Approximate search for this k-mer
        auto start_approx = std::chrono::high_resolution_clock::now();
        auto approx_results = trie.approximate_search(kmer_seq, max_mismatches);
        auto end_approx = std::chrono::high_resolution_clock::now();
        auto approx_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_approx - start_approx);
        
        std::cout << "  Approximate search: " << approx_results.size() << " matches"
                  << " (" << approx_time.count() << "ms)\n";
        
        // Show first few approximate matches
        if (!approx_results.empty()) {
            std::map<int, int> dist_counts;
            for (const auto& [seq, dist] : approx_results) {
                dist_counts[dist]++;
            }
            std::cout << "    By distance: ";
            for (const auto& [dist, count] : dist_counts) {
                std::cout << "d" << dist << "=" << count << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "Total time: " << total_time.count() << "ms\n";
    std::cout << "Breakdown: Load=" << load_time.count() << "ms, Build=" << build_time.count() 
              << "ms\n\n";
}


// Add this to your existing code - Levenshtein Automaton implementation
// Enhanced Myers algorithm with automaton-style optimizations

class LevenshteinAutomaton {
private:
    std::string query;
    int max_edit_distance;
    int query_length;
    
public:
    LevenshteinAutomaton(const std::string& query_str, int max_dist) 
        : query(query_str), max_edit_distance(max_dist), query_length(query_str.length()) {}
    
    // Test if a target string is accepted (within edit distance)
    bool accepts(const std::string& target) const {
        return get_edit_distance(target) >= 0;
    }
    
    // Get exact edit distance if within threshold, -1 otherwise
    int get_edit_distance(const std::string& target) const {
        int target_len = target.length();
        
        // Early termination based on length difference
        if (std::abs(target_len - query_length) > max_edit_distance) {
            return -1;
        }
        
        // Use two vectors for DP computation
        std::vector<int> prev_row(target_len + 1);
        std::vector<int> curr_row(target_len + 1);
        
        // Initialize first row (all insertions)
        for (int j = 0; j <= target_len; j++) {
            prev_row[j] = j;
        }
        
        // Process each character of query
        for (int i = 1; i <= query_length; i++) {
            curr_row[0] = i; // Deletion cost
            int min_in_row = i;
            
            for (int j = 1; j <= target_len; j++) {
                int cost = (query[i-1] == target[j-1]) ? 0 : 1;
                
                curr_row[j] = std::min({
                    prev_row[j] + 1,      // deletion
                    curr_row[j-1] + 1,    // insertion
                    prev_row[j-1] + cost  // substitution/match
                });
                
                min_in_row = std::min(min_in_row, curr_row[j]);
            }
            
            // Early termination: if minimum in row exceeds max distance
            if (min_in_row > max_edit_distance) {
                return -1;
            }
            
            std::swap(prev_row, curr_row);
        }
        
        return (prev_row[target_len] <= max_edit_distance) ? prev_row[target_len] : -1;
    }
};

class MyersAutomaton {
private:
    std::string query;
    int64_t query_bits;
    int query_length;
    int max_edit_distance;
    
    // Pre-computed pattern equality vectors (Myers optimization)
    int64_t Peq[4];
    int64_t pattern_mask;
    
public:
    MyersAutomaton(const std::string& query_str, int max_dist) 
        : query(query_str), max_edit_distance(max_dist) {
        
        int64_seq query_seq(query_str);
        query_bits = query_seq.bits[0];
        query_length = query_seq.length;
        
        // Pre-build pattern equality vectors (do this once, reuse many times)
        Peq[0] = Peq[1] = Peq[2] = Peq[3] = 0;
        for (int i = 0; i < query_length; i++) {
            int nuc = (query_bits >> (2 * i)) & 3;
            Peq[nuc] |= (1LL << i);
        }
        
        pattern_mask = (1LL << query_length) - 1;
    }
    
    // Ultra-fast acceptance test using your Myers algorithm
    bool accepts(const int64_seq& target) const {
        return get_edit_distance_fast(target) >= 0;
    }
    
    // Your Myers algorithm, but optimized for automaton use
    int get_edit_distance_fast(const int64_seq& target) const {
        if (target.bits.empty() || target.length > 32) return -1;
        
        int target_len = target.length;
        int64_t target_bits = target.bits[0];
        
        // Early termination based on length difference (automaton optimization)
        if (std::abs(target_len - query_length) > max_edit_distance) {
            return -1;
        }
        
        // Your Myers bit-parallel core (but optimized for partial matching)
        int64_t Pv = pattern_mask;
        int64_t Mv = 0;
        int score = query_length;
        int min_score = query_length; // Track minimum across all positions
        
        for (int j = 0; j < target_len; j++) {
            int target_nuc = (target_bits >> (2 * j)) & 3;
            int64_t Eq = Peq[target_nuc] & pattern_mask; // Use pre-computed Peq!
            
            // Myers core computation (your existing code)
            int64_t Xv = Eq | Mv;
            int64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
            int64_t Ph = Mv | ~(Xh | Pv);
            int64_t Mh = Pv & Xh;
            
            Ph &= pattern_mask;
            Mh &= pattern_mask;
            
            if (Ph & (1LL << (query_length - 1))) score++;
            if (Mh & (1LL << (query_length - 1))) score--;
            
            min_score = std::min(min_score, score);
            
            // Automaton-style early termination
            if (min_score == 0) return 0; // Perfect match found
            
            // Enhanced early termination: if impossible to get within max_dist
            int remaining = target_len - j - 1;
            if (score - remaining > max_edit_distance) {
                return -1; // Impossible to reach acceptable distance
            }
            
            Ph <<= 1;
            Pv = ((Mh << 1) | ~(Xv | Ph)) & pattern_mask;
            Mv = Ph & Xv;
        }
        
        return (min_score <= max_edit_distance) ? min_score : -1;
    }
    
    // Batch processing optimized for large whitelists
    std::vector<std::pair<int64_seq, int>> find_all_matches(const std::vector<int64_seq>& targets) const {
        std::vector<std::pair<int64_seq, int>> matches;
        matches.reserve(targets.size() / 100); // Rough estimate
        
        for (const auto& target : targets) {
            int dist = get_edit_distance_fast(target);
            if (dist >= 0) {
                matches.emplace_back(target, dist);
                
                // Optional: early exit if you only need a few matches
                // if (matches.size() >= some_limit) break;
            }
        }
        
        return matches;
    }
};

// Debug function to understand the differences
void debug_myers_vs_hybrid_difference(const std::string& query_seq, const std::string& whitelist_path, int max_dist = 2) {
    std::cout << "=== DEBUG: Myers vs Hybrid Differences ===\n";
    std::cout << "Query: " << query_seq << " (length: " << query_seq.length() << ")\n";
    std::cout << "Max distance: " << max_dist << "\n\n";
    
    // Load a small sample for detailed analysis
    auto whitelist = load_simple_wl(whitelist_path);
    std::vector<int64_seq> sample;
    
    // Take first 100 sequences for detailed analysis
    int count = 0;
    for (const auto& seq : whitelist) {
        sample.push_back(seq);
        if (++count >= 100) break;
    }
    
    std::cout << "Analyzing first " << sample.size() << " sequences...\n\n";
    
    int64_seq query_bits(query_seq);
    MyersAutomaton hybrid(query_seq, max_dist);
    
    std::vector<std::string> myers_only;
    std::vector<std::string> hybrid_only;
    std::vector<std::string> both_match;
    
    for (const auto& target : sample) {
        std::string target_str = target.bits_to_sequence();
        
        // Test Pure Myers
        int myers_dist = mutation_tools::bit_partial_match(
            query_bits.bits[0], target.bits[0], 
            query_bits.length, target.length, max_dist
        );
        
        // Test Hybrid
        int hybrid_dist = hybrid.get_edit_distance_fast(target);
        
        bool myers_match = (myers_dist >= 0);
        bool hybrid_match = (hybrid_dist >= 0);
        
        if (myers_match && hybrid_match) {
            both_match.push_back(target_str + " (M:" + std::to_string(myers_dist) + ", H:" + std::to_string(hybrid_dist) + ")");
        } else if (myers_match && !hybrid_match) {
            myers_only.push_back(target_str + " (M:" + std::to_string(myers_dist) + ", len:" + std::to_string(target.length) + ")");
        } else if (!myers_match && hybrid_match) {
            hybrid_only.push_back(target_str + " (H:" + std::to_string(hybrid_dist) + ", len:" + std::to_string(target.length) + ")");
        }
    }
    
    std::cout << "Results:\n";
    std::cout << "Both match: " << both_match.size() << "\n";
    std::cout << "Myers only: " << myers_only.size() << "\n";
    std::cout << "Hybrid only: " << hybrid_only.size() << "\n\n";
    
    // Show examples of differences
    if (!myers_only.empty()) {
        std::cout << "Examples where ONLY Myers matches (partial matching):\n";
        for (size_t i = 0; i < std::min((size_t)5, myers_only.size()); i++) {
            std::cout << "  " << myers_only[i] << "\n";
        }
        std::cout << "\n";
    }
    
    if (!hybrid_only.empty()) {
        std::cout << "Examples where ONLY Hybrid matches (should be rare):\n";
        for (size_t i = 0; i < std::min((size_t)5, hybrid_only.size()); i++) {
            std::cout << "  " << hybrid_only[i] << "\n";
        }
        std::cout << "\n";
    }
    
    if (!both_match.empty()) {
        std::cout << "Examples where both match:\n";
        for (size_t i = 0; i < std::min((size_t)5, both_match.size()); i++) {
            std::cout << "  " << both_match[i] << "\n";
        }
        std::cout << "\n";
    }
    
    // Analyze length distribution of mismatches
    if (!myers_only.empty()) {
        std::cout << "Length analysis of Myers-only matches:\n";
        std::map<int, int> length_counts;
        
        for (const auto& target : sample) {
            int myers_dist = mutation_tools::bit_partial_match(
                query_bits.bits[0], target.bits[0], 
                query_bits.length, target.length, max_dist
            );
            int hybrid_dist = hybrid.get_edit_distance_fast(target);
            
            if (myers_dist >= 0 && hybrid_dist < 0) {
                length_counts[target.length]++;
            }
        }
        
        for (const auto& [length, count] : length_counts) {
            int length_diff = std::abs(static_cast<int>(length) - static_cast<int>(query_seq.length()));
            std::cout << "  Length " << length << " (diff=" << length_diff << "): " << count << " sequences\n";
        }
    }
}

// Fixed hybrid implementation that matches Myers behavior
class MyersAutomatonFixed {
private:
    std::string query;
    int64_t query_bits;
    int query_length;
    int max_edit_distance;
    
    int64_t Peq[4];
    int64_t pattern_mask;
    
public:
    MyersAutomatonFixed(const std::string& query_str, int max_dist) 
        : query(query_str), max_edit_distance(max_dist) {
        
        int64_seq query_seq(query_str);
        query_bits = query_seq.bits[0];
        query_length = query_seq.length;
        
        // Pre-build pattern equality vectors
        Peq[0] = Peq[1] = Peq[2] = Peq[3] = 0;
        for (int i = 0; i < query_length; i++) {
            int nuc = (query_bits >> (2 * i)) & 3;
            Peq[nuc] |= (1LL << i);
        }
        
        pattern_mask = (1LL << query_length) - 1;
    }
    
    // Modified to match the partial matching behavior of Pure Myers
    int get_edit_distance_fast(const int64_seq& target) const {
        if (target.bits.empty() || target.length > 32) return -1;
        
        int target_len = target.length;
        int64_t target_bits = target.bits[0];
        
        // Remove the length difference check to match partial matching behavior
        // if (std::abs(target_len - query_length) > max_edit_distance) {
        //     return -1;
        // }
        
        // Use the EXACT same logic as bit_partial_match
        int64_t Pv = pattern_mask;
        int64_t Mv = 0;
        int score = query_length;
        int min_score = query_length; // This is the key - track minimum across ALL positions
        
        for (int j = 0; j < target_len; j++) {
            int target_nuc = (target_bits >> (2 * j)) & 3;
            int64_t Eq = Peq[target_nuc] & pattern_mask;
            
            // Myers core computation (identical to your bit_partial_match)
            int64_t Xv = Eq | Mv;
            int64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
            int64_t Ph = Mv | ~(Xh | Pv);
            int64_t Mh = Pv & Xh;
            
            Ph &= pattern_mask;
            Mh &= pattern_mask;
            
            if (Ph & (1LL << (query_length - 1))) score++;
            if (Mh & (1LL << (query_length - 1))) score--;
            
            min_score = std::min(min_score, score); // Key: partial matching behavior
            
            if (min_score == 0) return 0; // Perfect match found
            
            Ph <<= 1;
            Pv = ((Mh << 1) | ~(Xv | Ph)) & pattern_mask;
            Mv = Ph & Xv;
        }
        
        return (min_score <= max_edit_distance) ? min_score : -1;
    }
};

// Test function comparing all three approaches
void test_myers_automaton_hybrid(const std::string& query_seq, const std::string& whitelist_path, int max_dist = 2) {
    auto start_total = std::chrono::high_resolution_clock::now();
    
    std::cout << "=== Three-Way Algorithm Comparison ===\n";
    std::cout << "Query: " << query_seq << "\n";
    std::cout << "Max distance: " << max_dist << "\n";
    
    // Load whitelist
    auto start_load = std::chrono::high_resolution_clock::now();
    auto whitelist = load_simple_wl(whitelist_path);
    auto end_load = std::chrono::high_resolution_clock::now();
    auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load);
    
    std::cout << "Loaded " << whitelist.size() << " sequences (" << load_time.count() << "ms)\n";
    
    // Convert to vector for consistent testing
    std::vector<int64_seq> whitelist_vec(whitelist.begin(), whitelist.end());
    std::vector<std::string> whitelist_strings;
    for (const auto& seq : whitelist_vec) {
        whitelist_strings.push_back(seq.bits_to_sequence());
    }
    
    // Test 1: Pure Myers (your current approach)
    std::cout << "\n--- Testing Pure Myers ---\n";
    int64_seq query_bits(query_seq);
    std::vector<std::pair<std::string, int>> myers_matches;
    
    auto start_myers = std::chrono::high_resolution_clock::now();
    
    for (const auto& target_seq : whitelist_vec) {
        int dist = mutation_tools::bit_partial_match(
            query_bits.bits[0], target_seq.bits[0], 
            query_bits.length, target_seq.length, max_dist
        );
        if (dist >= 0) {
            myers_matches.emplace_back(target_seq.bits_to_sequence(), dist);
        }
    }
    
    auto end_myers = std::chrono::high_resolution_clock::now();
    auto myers_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_myers - start_myers);
    
    std::cout << "Myers found " << myers_matches.size() << " matches in " 
              << myers_time.count() << "ms\n";
    
    // Test 2: Myers Automaton Hybrid
    std::cout << "\n--- Testing Myers Automaton Hybrid ---\n";
    MyersAutomaton hybrid(query_seq, max_dist);
    
    auto start_hybrid = std::chrono::high_resolution_clock::now();
    auto hybrid_matches = hybrid.find_all_matches(whitelist_vec);
    auto end_hybrid = std::chrono::high_resolution_clock::now();
    auto hybrid_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_hybrid - start_hybrid);
    
    std::cout << "Hybrid found " << hybrid_matches.size() << " matches in " 
              << hybrid_time.count() << "ms\n";
    
    // Test 3: Traditional DP Automaton (for comparison)
    std::cout << "\n--- Testing Traditional DP Automaton ---\n";
    LevenshteinAutomaton traditional(query_seq, max_dist);
    std::vector<std::pair<std::string, int>> traditional_matches;
    
    auto start_traditional = std::chrono::high_resolution_clock::now();
    
    for (const auto& target : whitelist_strings) {
        int dist = traditional.get_edit_distance(target);
        if (dist >= 0) {
            traditional_matches.emplace_back(target, dist);
        }
    }
    
    auto end_traditional = std::chrono::high_resolution_clock::now();
    auto traditional_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_traditional - start_traditional);
    
    std::cout << "Traditional found " << traditional_matches.size() << " matches in " 
              << traditional_time.count() << "ms\n";
    
    // Performance comparison
    std::cout << "\n--- Performance Comparison ---\n";
    std::cout << "Pure Myers:      " << (myers_time.count() / 1000.0) << "s\n";
    std::cout << "Myers Hybrid:    " << (hybrid_time.count() / 1000.0) << "s\n";
    std::cout << "Traditional DP:  " << (traditional_time.count() / 1000.0) << "s\n";
    
    double hybrid_speedup = (double)myers_time.count() / hybrid_time.count();
    double traditional_speedup = (double)traditional_time.count() / hybrid_time.count();
    
    std::cout << "Hybrid vs Pure Myers: " << hybrid_speedup << "x speedup\n";
    std::cout << "Hybrid vs Traditional: " << traditional_speedup << "x speedup\n";
    
    // Verify results match
    std::sort(myers_matches.begin(), myers_matches.end());
    
    std::vector<std::pair<std::string, int>> hybrid_string_matches;
    for (const auto& match : hybrid_matches) {
        hybrid_string_matches.emplace_back(match.first.bits_to_sequence(), match.second);
    }
    std::sort(hybrid_string_matches.begin(), hybrid_string_matches.end());
    
    bool results_match = (myers_matches.size() == hybrid_string_matches.size());
    std::cout << "Results match Myers: " << (results_match ? "YES" : "NO") << "\n";
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "Total test time: " << total_time.count() << "ms\n\n";
}

// Add this function to test the hybrid approach
void run_hybrid_automaton_tests() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "MYERS AUTOMATON HYBRID PERFORMANCE TESTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::string whitelist_path = "/Users/cmv/Desktop/rad_paper/cellranger_whitelists/3M-february-2018-3v3.txt.gz";
    
    // Test with your typical barcode
    test_myers_automaton_hybrid("AAACCCAAGAAGATCT", whitelist_path, 2);
}

// Whitelist-based Levenshtein Automaton for barcode demultiplexing
// Build automata for each whitelist barcode, then test reads against all of them

class WhitelistAutomaton {
private:
    std::string barcode;
    int64_t barcode_bits;
    int barcode_length;
    int max_edit_distance;
    
    // Pre-computed pattern equality vectors (Myers optimization)
    int64_t Peq[4];
    int64_t pattern_mask;
    
public:
    WhitelistAutomaton(const int64_seq& barcode_seq, int max_dist) 
        : max_edit_distance(max_dist) {
        
        barcode = barcode_seq.bits_to_sequence();
        barcode_bits = barcode_seq.bits[0];
        barcode_length = barcode_seq.length;
        
        // Pre-build pattern equality vectors for this specific barcode
        Peq[0] = Peq[1] = Peq[2] = Peq[3] = 0;
        for (int i = 0; i < barcode_length; i++) {
            int nuc = (barcode_bits >> (2 * i)) & 3;
            Peq[nuc] |= (1LL << i);
        }
        
        pattern_mask = (1LL << barcode_length) - 1;
    }
    
    // Test if a read matches this barcode within edit distance
    bool accepts(const int64_seq& read) const {
        return get_edit_distance(read) >= 0;
    }
    
    // Get exact edit distance if within threshold, -1 otherwise
    int get_edit_distance(const int64_seq& read) const {
        if (read.bits.empty() || read.length > 32) return -1;
        
        int read_len = read.length;
        int64_t read_bits = read.bits[0];
        
        // Myers bit-parallel algorithm optimized for this specific barcode
        int64_t Pv = pattern_mask;
        int64_t Mv = 0;
        int score = barcode_length;
        int min_score = barcode_length; // Track minimum for partial matching
        
        for (int j = 0; j < read_len; j++) {
            int read_nuc = (read_bits >> (2 * j)) & 3;
            int64_t Eq = Peq[read_nuc] & pattern_mask; // Use pre-computed Peq!
            
            // Myers core computation
            int64_t Xv = Eq | Mv;
            int64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
            int64_t Ph = Mv | ~(Xh | Pv);
            int64_t Mh = Pv & Xh;
            
            Ph &= pattern_mask;
            Mh &= pattern_mask;
            
            if (Ph & (1LL << (barcode_length - 1))) score++;
            if (Mh & (1LL << (barcode_length - 1))) score--;
            
            min_score = std::min(min_score, score);
            
            // Early termination optimizations
            if (min_score == 0) return 0; // Perfect match
            if (min_score > max_edit_distance && 
                (read_len - j - 1) < (min_score - max_edit_distance)) {
                return -1; // Impossible to recover
            }
            
            Ph <<= 1;
            Pv = ((Mh << 1) | ~(Xv | Ph)) & pattern_mask;
            Mv = Ph & Xv;
        }
        
        return (min_score <= max_edit_distance) ? min_score : -1;
    }
    
    const std::string& get_barcode() const { return barcode; }
};

// Whitelist manager - builds and manages all barcode automata
class WhitelistAutomatonSet {
private:
    std::vector<WhitelistAutomaton> automata;
    std::vector<int64_seq> original_barcodes;
    int max_edit_distance;
    
public:
    WhitelistAutomatonSet(const std::unordered_set<int64_seq>& whitelist, int max_dist) 
        : max_edit_distance(max_dist) {
        
        std::cout << "Building automata for " << whitelist.size() << " barcodes...\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        automata.reserve(whitelist.size());
        original_barcodes.reserve(whitelist.size());
        
        for (const auto& barcode : whitelist) {
            automata.emplace_back(barcode, max_dist);
            original_barcodes.push_back(barcode);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto build_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Automata built in " << build_time.count() << "ms\n";
    }
    
    // Find the best matching barcode(s) for a read
    struct MatchResult {
        int64_seq barcode;
        int distance;
        size_t barcode_index;
    };
    
    std::vector<MatchResult> find_matches(const int64_seq& read) const {
        std::vector<MatchResult> matches;
        
        for (size_t i = 0; i < automata.size(); i++) {
            int dist = automata[i].get_edit_distance(read);
            if (dist >= 0) {
                matches.push_back({original_barcodes[i], dist, i});
            }
        }
        
        // Sort by distance (best matches first)
        std::sort(matches.begin(), matches.end(), 
                 [](const MatchResult& a, const MatchResult& b) {
                     return a.distance < b.distance;
                 });
        
        return matches;
    }
    
    // Find only the best match (most common use case)
    std::optional<MatchResult> find_best_match(const int64_seq& read) const {
        std::optional<MatchResult> best_match;
        int best_distance = max_edit_distance + 1;
        
        for (size_t i = 0; i < automata.size(); i++) {
            int dist = automata[i].get_edit_distance(read);
            if (dist >= 0 && dist < best_distance) {
                best_distance = dist;
                best_match = {original_barcodes[i], dist, i};
                
                if (dist == 0) break; // Perfect match, can't do better
            }
        }
        
        return best_match;
    }
    
    // Batch processing for multiple reads
    std::vector<std::optional<MatchResult>> batch_process(const std::vector<int64_seq>& reads) const {
        std::vector<std::optional<MatchResult>> results;
        results.reserve(reads.size());
        
        for (const auto& read : reads) {
            results.push_back(find_best_match(read));
        }
        
        return results;
    }
    
    size_t size() const { return automata.size(); }
};

// Performance testing function
void test_whitelist_automaton_performance(const std::string& whitelist_path, 
                                        const std::vector<std::string>& test_reads,
                                        int max_dist) {
    
    std::cout << "=== Whitelist Automaton Performance Test ===\n";
    std::cout << "Max edit distance: " << max_dist << "\n";
    std::cout << "Test reads: " << test_reads.size() << "\n";
    
    // Load whitelist
    auto start_load = std::chrono::high_resolution_clock::now();
    auto whitelist = load_simple_wl(whitelist_path);
    auto end_load = std::chrono::high_resolution_clock::now();
    auto load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load);
    
    std::cout << "Loaded " << whitelist.size() << " barcodes (" << load_time.count() << "ms)\n";
    
    // Build automaton set
    WhitelistAutomatonSet automaton_set(whitelist, max_dist);
    
    // Convert test reads
    std::vector<int64_seq> test_read_seqs;
    for (const auto& read_str : test_reads) {
        test_read_seqs.emplace_back(read_str);
    }
    
    // Test individual read processing
    std::cout << "\n--- Individual Read Processing ---\n";
    auto start_individual = std::chrono::high_resolution_clock::now();
    
    int total_matches = 0;
    std::map<int, int> distance_counts;
    
    for (const auto& read : test_read_seqs) {
        auto best_match = automaton_set.find_best_match(read);
        if (best_match.has_value()) {
            total_matches++;
            distance_counts[best_match->distance]++;
        }
    }
    
    auto end_individual = std::chrono::high_resolution_clock::now();
    auto individual_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_individual - start_individual);
    
    std::cout << "Individual processing: " << individual_time.count() << "ms total\n";
    std::cout << "Average per read: " << (individual_time.count() / test_reads.size()) << "ms\n";
    std::cout << "Matches found: " << total_matches << "/" << test_reads.size() << "\n";
    std::cout << "Distance distribution: ";
    for (const auto& [dist, count] : distance_counts) {
        std::cout << "d" << dist << "=" << count << " ";
    }
    std::cout << "\n";
    
    // Test batch processing
    std::cout << "\n--- Batch Processing ---\n";
    auto start_batch = std::chrono::high_resolution_clock::now();
    auto batch_results = automaton_set.batch_process(test_read_seqs);
    auto end_batch = std::chrono::high_resolution_clock::now();
    auto batch_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - start_batch);
    
    std::cout << "Batch processing: " << batch_time.count() << "ms total\n";
    std::cout << "Average per read: " << (batch_time.count() / test_reads.size()) << "ms\n";
    
    // Performance metrics
    double reads_per_second = (test_reads.size() * 1000.0) / individual_time.count();
    std::cout << "Throughput: " << std::fixed << std::setprecision(0) << reads_per_second << " reads/second\n";
    
    // Show some example matches
    std::cout << "\n--- Example Matches ---\n";
    for (size_t i = 0; i < std::min((size_t)5, batch_results.size()); i++) {
        std::cout << "Read: " << test_reads[i];
        if (batch_results[i].has_value()) {
            std::cout << " -> " << batch_results[i]->barcode.bits_to_sequence() 
                      << " (distance: " << batch_results[i]->distance << ")\n";
        } else {
            std::cout << " -> No match\n";
        }
    }
    std::cout << "\n";
}

// Updated comprehensive test function
void run_whitelist_automaton_tests() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "WHITELIST AUTOMATON PERFORMANCE TESTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::string whitelist_path = "/Users/cmv/Desktop/rad_paper/cellranger_whitelists/3M-february-2018-3v3.txt.gz";
    
    // Create test reads - mix of exact matches, close matches, and no matches
    std::vector<std::string> test_reads = {
        "AAACCCAAGAAGATCT",  // Should match something in whitelist
        "AAACCCAAGAAGATCC",  // 1 edit from above
        "AAACCCAAGAAGATAA",  // 2 edits from first
        "CGATCTAAGCCAAGAAGATCTGGTGTCG", // Your longer test sequence
        "AAGCCAAGAAGATCTG",  // Substring match
        "TTTTTTTTTTTTTTTT",  // Probably no match
        "AAACCCACAGATCGTT",  // Another test case
        "AAACCCAAAGATCTTTCGTGTATTTT", // Different length
    };
    
    // Test with different edit distances
    for (int max_dist = 1; max_dist <= 3; max_dist++) {
        std::cout << "\n" << std::string(40, '-') << "\n";
        std::cout << "Testing with max_dist = " << max_dist << "\n";
        std::cout << std::string(40, '-') << "\n";
        
        test_whitelist_automaton_performance(whitelist_path, test_reads, max_dist);
    }
}

// Generate all sequences within a given Levenshtein distance of a barcode

#include <string>
#include <unordered_set>
#include <vector>
#include <chrono>
#include <iostream>

class LevenshteinGenerator {
public:
    // Generate all sequences within max_distance edits of the input barcode
    static std::unordered_set<std::string> generate_all_edits(const std::string& barcode, int max_distance) {
        std::unordered_set<std::string> all_sequences;
        generate_edits_recursive(barcode, max_distance, all_sequences);
        return all_sequences;
    }

private:
    static void generate_edits_recursive(const std::string& sequence, int remaining_distance,
                                       std::unordered_set<std::string>& results) {
        // Add current sequence to results
        results.insert(sequence);
        
        // Base case: no more edits allowed
        if (remaining_distance == 0) {
            return;
        }
        
        // Generate all possible single edits from current sequence
        std::vector<std::string> single_edits = generate_single_edits(sequence);
        
        // Recursively generate edits from each single edit
        for (const std::string& edited_seq : single_edits) {
            // Only recurse if we haven't seen this sequence before (pruning)
            if (results.find(edited_seq) == results.end()) {
                generate_edits_recursive(edited_seq, remaining_distance - 1, results);
            }
        }
    }
    
    static std::vector<std::string> generate_single_edits(const std::string& sequence) {
        std::vector<std::string> edits;
        const char nucleotides[] = {'A', 'T', 'C', 'G'};
        
        // Substitutions: replace each character with each possible nucleotide
        for (size_t i = 0; i < sequence.length(); i++) {
            for (char nuc : nucleotides) {
                if (nuc != sequence[i]) {  // Only generate actual changes
                    std::string substituted = sequence;
                    substituted[i] = nuc;
                    edits.push_back(substituted);
                }
            }
        }
        
        // Deletions: remove each character
        for (size_t i = 0; i < sequence.length(); i++) {
            std::string deleted = sequence.substr(0, i) + sequence.substr(i + 1);
            edits.push_back(deleted);
        }
        
        // Insertions: insert each nucleotide at each possible position
        for (size_t i = 0; i <= sequence.length(); i++) {
            for (char nuc : nucleotides) {
                std::string inserted = sequence.substr(0, i) + nuc + sequence.substr(i);
                edits.push_back(inserted);
            }
        }
        
        return edits;
    }
};

// Optimized version that avoids string operations using bit encoding
class BitLevenshteinGenerator {
private:
    static const char nucleotides[4];
    
public:
    // Generate all int64_seq within max_distance edits
    static std::unordered_set<int64_seq> generate_all_edits(const int64_seq& barcode, int max_distance) {
        std::unordered_set<int64_seq> all_sequences;
        generate_edits_recursive(barcode, max_distance, all_sequences);
        return all_sequences;
    }

private:
    static void generate_edits_recursive(const int64_seq& sequence, int remaining_distance,
                                       std::unordered_set<int64_seq>& results) {
        results.insert(sequence);
        
        if (remaining_distance == 0) {
            return;
        }
        
        std::vector<int64_seq> single_edits = generate_single_edits(sequence);
        
        for (const int64_seq& edited_seq : single_edits) {
            if (results.find(edited_seq) == results.end()) {
                generate_edits_recursive(edited_seq, remaining_distance - 1, results);
            }
        }
    }
    
    static std::vector<int64_seq> generate_single_edits(const int64_seq& sequence) {
        std::vector<int64_seq> edits;
        
        if (sequence.bits.empty() || sequence.length > 32) {
            return edits;
        }
        
        int64_t bits = sequence.bits[0];
        int len = sequence.length;
        
        // Substitutions
        for (int pos = 0; pos < len; pos++) {
            int current_nuc = (bits >> (2 * pos)) & 3;
            
            for (int new_nuc = 0; new_nuc < 4; new_nuc++) {
                if (new_nuc != current_nuc) {
                    int64_t new_bits = bits;
                    // Clear the 2 bits at this position
                    new_bits &= ~(3LL << (2 * pos));
                    // Set the new nucleotide
                    new_bits |= ((int64_t)new_nuc << (2 * pos));
                    
                    int64_seq new_seq;
                    new_seq.bits.push_back(new_bits);
                    new_seq.length = len;
                    edits.push_back(new_seq);
                }
            }
        }
        
        // Deletions
        for (int pos = 0; pos < len; pos++) {
            if (len > 1) {  // Don't delete if it would make empty sequence
                int64_t new_bits = 0;
                int new_pos = 0;
                
                // Copy all nucleotides except the one at 'pos'
                for (int i = 0; i < len; i++) {
                    if (i != pos) {
                        int nuc = (bits >> (2 * i)) & 3;
                        new_bits |= ((int64_t)nuc << (2 * new_pos));
                        new_pos++;
                    }
                }
                
                int64_seq new_seq;
                new_seq.bits.push_back(new_bits);
                new_seq.length = len - 1;
                edits.push_back(new_seq);
            }
        }
        
        // Insertions
        for (int pos = 0; pos <= len; pos++) {
            for (int new_nuc = 0; new_nuc < 4; new_nuc++) {
                if (len < 32) {  // Don't exceed maximum length
                    int64_t new_bits = 0;
                    int new_pos = 0;
                    
                    // Copy nucleotides before insertion point
                    for (int i = 0; i < pos; i++) {
                        int nuc = (bits >> (2 * i)) & 3;
                        new_bits |= ((int64_t)nuc << (2 * new_pos));
                        new_pos++;
                    }
                    
                    // Insert new nucleotide
                    new_bits |= ((int64_t)new_nuc << (2 * new_pos));
                    new_pos++;
                    
                    // Copy nucleotides after insertion point
                    for (int i = pos; i < len; i++) {
                        int nuc = (bits >> (2 * i)) & 3;
                        new_bits |= ((int64_t)nuc << (2 * new_pos));
                        new_pos++;
                    }
                    
                    int64_seq new_seq;
                    new_seq.bits.push_back(new_bits);
                    new_seq.length = len + 1;
                    edits.push_back(new_seq);
                }
            }
        }
        
        return edits;
    }
};

// Test and benchmark the generator
void test_levenshtein_generator() {
    std::cout << "=== Levenshtein Generator Test ===\n";
    
    std::string test_barcode = "AAACCCAAGAAGATCT";
    
    for (int distance = 0; distance <= 2; distance++) {
        std::cout << "\n--- Distance " << distance << " ---\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Test string version
        auto string_results = LevenshteinGenerator::generate_all_edits(test_barcode, distance);
        
        auto mid_time = std::chrono::high_resolution_clock::now();
        
        // Test bit version
        int64_seq bit_barcode(test_barcode);
        auto bit_results = BitLevenshteinGenerator::generate_all_edits(bit_barcode, distance);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        
        auto string_duration = std::chrono::duration_cast<std::chrono::microseconds>(mid_time - start_time);
        auto bit_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - mid_time);
        
        std::cout << "String version: " << string_results.size() << " sequences in " 
                  << string_duration.count() << "μs\n";
        std::cout << "Bit version: " << bit_results.size() << " sequences in " 
                  << bit_duration.count() << "μs\n";
        std::cout << "Speedup: " << (double)string_duration.count() / bit_duration.count() << "x\n";
        
        // Verify results match
        std::unordered_set<std::string> bit_strings;
        for (const auto& seq : bit_results) {
            bit_strings.insert(seq.bits_to_sequence());
        }
        
        bool results_match = (string_results.size() == bit_strings.size());
        if (results_match) {
            for (const auto& str : string_results) {
                if (bit_strings.find(str) == bit_strings.end()) {
                    results_match = false;
                    break;
                }
            }
        }
        
        std::cout << "Results match: " << (results_match ? "YES" : "NO") << "\n";
        
        // Show examples for each distance
        if (distance <= 2) {
            std::cout << "First 50 examples:\n";
            int shown = 0;
            for (const auto& seq : string_results) {
                if (seq != test_barcode && shown < 50) {
                    std::cout << "  " << seq << "\n";
                    shown++;
                }
            }
        }
    }
}

// Memory usage estimator
void estimate_memory_usage(const std::string& barcode, int max_distance) {
    std::cout << "\n=== Memory Usage Estimation ===\n";
    std::cout << "Barcode: " << barcode << " (length: " << barcode.length() << ")\n";
    std::cout << "Max distance: " << max_distance << "\n\n";
    
    for (int dist = 0; dist <= max_distance; dist++) {
        auto sequences = LevenshteinGenerator::generate_all_edits(barcode, dist);
        
        size_t total_chars = 0;
        for (const auto& seq : sequences) {
            total_chars += seq.length();
        }
        
        size_t estimated_memory = sequences.size() * 64 + total_chars; // rough estimate
        
        std::cout << "Distance " << dist << ":\n";
        std::cout << "  Sequences: " << sequences.size() << "\n";
        std::cout << "  Total characters: " << total_chars << "\n";
        std::cout << "  Estimated memory: ~" << (estimated_memory / 1024) << " KB\n\n";
    }
    
    // Extrapolate to full whitelist
    auto full_sequences = LevenshteinGenerator::generate_all_edits(barcode, max_distance);
    size_t per_barcode_memory = full_sequences.size() * 64;
    size_t total_barcodes = 737000; // Approximate 10x whitelist size
    size_t total_memory_gb = (per_barcode_memory * total_barcodes) / (1024 * 1024 * 1024);
    
    std::cout << "Extrapolation to full 10x whitelist (737K barcodes):\n";
    std::cout << "Total estimated memory: ~" << total_memory_gb << " GB\n";
}

// Main test function
void run_levenshtein_generator_tests() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "LEVENSHTEIN EDIT GENERATOR TESTS\n";
    std::cout << std::string(60, '=') << "\n";
    
    test_levenshtein_generator();
    estimate_memory_usage("AAACCCAAGAAGATCT", 2);
}


// C++ Comparison Function
void compare_mutation_functions(const std::string& test_barcode, int max_rounds = 2, int iterations = 10) {
    std::cout << "=== Mutation Function Comparison ===\n";
    std::cout << "Barcode: " << test_barcode << " (length: " << test_barcode.length() << ")\n";
    std::cout << "Iterations: " << iterations << "\n\n";
    
    int64_seq seq(test_barcode);
    
    for (int rounds = 1; rounds <= max_rounds; ++rounds) {
        std::cout << "--- Rounds " << rounds << " ---\n";
        
        // Test generate_mutated_barcodes
        std::vector<double> mutated_times;
        mutated_times.reserve(iterations);
        size_t mutated_size = 0;
        
        for (int i = 0; i < iterations; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            auto results = mutation_tools::generate_mutated_barcodes(seq, rounds);
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            mutated_times.push_back(duration.count());
            if (i == 0) mutated_size = results.size();
        }
        
        // Test generate_lv_barcodes  
        std::vector<double> lv_times;
        lv_times.reserve(iterations);
        size_t lv_size = 0;
        
        for (int i = 0; i < iterations; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            auto results = mutation_tools::generate_lv_barcodes(seq, rounds);
            auto end = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            lv_times.push_back(duration.count());
            if (i == 0) lv_size = results.size();
        }
        
        // Calculate averages
        double mut_avg = std::accumulate(mutated_times.begin(), mutated_times.end(), 0.0) / iterations;
        double lv_avg = std::accumulate(lv_times.begin(), lv_times.end(), 0.0) / iterations;
        
        std::cout << "generate_mutated_barcodes: " << mutated_size << " sequences, " 
                  << mut_avg << "μs avg\n";
        std::cout << "generate_lv_barcodes: " << lv_size << " sequences, " 
                  << lv_avg << "μs avg\n";
        
        if (lv_avg < mut_avg) {
            std::cout << "Levenshtein is " << (mut_avg / lv_avg) << "x faster\n";
        } else {
            std::cout << "Mutated is " << (lv_avg / mut_avg) << "x faster\n";
        }
        
        std::cout << "Size ratio (LV/Mut): " << ((double)lv_size / mutated_size) << "\n\n";
    }
}

int main() {
    std::cout << "Testing barcode edit distance functions...\n\n";
    
    // Test 1: Single sequence distance
    std::cout << "=== Test 1: Single Sequence Distance ===\n";
    int64_seq query("TATCATCGATCGATCG");
    int64_seq target1("TATCATCGATCGATCG");  // exact match
    int64_seq target2("TATCATCGATCGATCC");  // 1 substitution
    int64_seq target3("TATCGTCGATCGATCG");  // 1 substitution
    int64_seq target4("TATCTTCGATCGATCG");  // 1 substitution
    
    std::cout << "Query: " << query.bits_to_sequence() << "\n";
    std::cout << "Distance to " << target1.bits_to_sequence() << ": " 
              << mutation_tools::int64_lvdist(query, target1, 4) << "\n";
    std::cout << "Distance to " << target2.bits_to_sequence() << ": " 
              << mutation_tools::int64_lvdist(query, target2, 4) << "\n";
    std::cout << "Distance to " << target3.bits_to_sequence() << ": " 
              << mutation_tools::int64_lvdist(query, target3, 4) << "\n";
    std::cout << "Distance to " << target4.bits_to_sequence() << ": " 
              << mutation_tools::int64_lvdist(query, target4, 4) << "\n";
    
    // Test 2: Multiple sequence distance with barcode_entry vector
    std::cout << "\n=== Test 2: Multiple Sequence Distance ===\n";
    int64_seq query2("TAAAACCCGGGTTTAG");
    std::vector<barcode_entry> targets;
    
    // Create some test barcodes
    barcode_entry be1; be1.barcode = int64_seq("TAAAACCCGGGTTTAG");  // exact
    barcode_entry be2; be2.barcode = int64_seq("TAAAACCCGGGTTCAG");  // 1 diff
    barcode_entry be3; be3.barcode = int64_seq("TAAAACCCGGGTCCAG");  // 2 diff
    barcode_entry be4; be4.barcode = int64_seq("TAAAACCCGGGCCCAG");  // 3 diff
    
    targets.push_back(be1);
    targets.push_back(be2);
    targets.push_back(be3);
    targets.push_back(be4);
    
    std::cout << "Query: " << query2.bits_to_sequence() << "\n";
    auto results = mutation_tools::int64_lvdist(query2, targets, 5);
    
    for (const auto& [dist, barcodes] : results) {
        std::cout << "Distance " << dist << ": " << barcodes.size() << " barcode(s)\n";
        for (const auto& bc : barcodes) {
            std::cout << "  - " << bc.bits_to_sequence() << "\n";
        }
    }
    
    // Test 3: Multiple sequence distance with unordered_set
    std::cout << "\n=== Test 3: Unordered Set Distance ===\n";
    int64_seq query3("ATCG");
    std::unordered_set<int64_seq> target_set = {
        int64_seq("ATCG"),  // distance 0
        int64_seq("ATCC"),  // distance 1
        int64_seq("TTCG"),  // distance 1
        int64_seq("ATGG"),  // distance 1
        int64_seq("GGGG"),  // distance 3
    };
    
    std::cout << "Query: " << query3.bits_to_sequence() << "\n";
    auto results3 = mutation_tools::int64_lvdist(query3, target_set, 5);
    
    for (const auto& [dist, barcodes] : results3) {
        std::cout << "Distance " << dist << ": " << barcodes.size() << " barcode(s)\n";
        for (const auto& bc : barcodes) {
            std::cout << "  - " << bc.bits_to_sequence() << "\n";
        }
    }
    

    std::cout << "=== DEBUGGING BARCODE DISTANCE CALCULATION ===\n\n";
    
    // First, let's check if encoding/decoding works correctly
    std::cout << "=== Step 1: Check sequence encoding/decoding ===\n";
    debug_sequence_encoding("ATCG");
    debug_sequence_encoding("TATCGTCGATCGATCG");
    debug_sequence_encoding("TAAAACCCGGGTTCAG");
    
    // Test your problematic cases
    std::cout << "=== Step 2: Test problematic cases ===\n";
    
    // Case 1: These should be 1 substitution but showing 0
    test_bit_ld("AATCGTCGATCGATCG", "TATCGTCGATCGATCG");  // A->T at position 0
    test_bit_ld("AATCGTCGATCGATCG", "TATCTTCGATCGATCG");  // A->T at pos 0, G->T at pos 5
    
    // Case 2: These should be 1, 2, 3 but showing 2
    test_bit_ld("TAAAACCCGGGTTCAG", "TAAAACCCGGGTTCAG");  // exact match (should be 0)
    test_bit_ld("TAAAACCCGGGTTCAG", "TAAAACCCGGGTTCAG");  // 1 diff
    test_bit_ld("TAAAACCCGGGTTCAG", "TAAAACCCGGGTCCAG");  // 2 diff
    test_bit_ld("TAAAACCCGGGTTCAG", "TAAAACCCGGGCCCAG");  // 3 diff
    
    // Test simple cases to verify basic functionality
    std::cout << "=== Step 3: Test simple cases ===\n";
    test_bit_ld("AAAA", "AAAA");  // exact match
    test_bit_ld("AAAA", "AAAT");  // 1 substitution
    test_bit_ld("AAAA", "AATT");  // 2 substitutions
    test_bit_ld("AAAA", "ATTT");  // 3 substitutions
    test_bit_ld("AAAA", "TTTT");  // 4 substitutions

    test_bit_ld("CCCC", "CCCG");  // 1 substitution
    test_bit_ld("CCCC", "CCGG");  // 2 substitutions
    test_bit_ld("CCCC", "CGGG");  // 3 substitutions
    test_bit_ld("CCCC", "GGGG");  // 4 substitutions

    // Test edge cases
    std::cout << "=== Step 4: Test edge cases ===\n";
    test_bit_ld("A", "A");    // single base exact
    test_bit_ld("A", "T");    // single base diff
    test_bit_ld("AT", "AT");  // two base exact
    test_bit_ld("AT", "AC");  // two base diff
    std::cout << "\nTest, Version 1: completed!\n";

    // Run comprehensive tests
    std::cout << "=== Step 5: Run comprehensive tests ===\n";
    run_comprehensive_tests();
    std::cout << "Comprehensive tests completed!\n";

    //test_bit_trie("CGATCTAAGCCAAGAAGATCTGGTGTCG", "/Users/cmv/Desktop/rad_paper/cellranger_whitelists/3M-february-2018-3v3.txt.gz", 2);
    std::cout << "Running automata tests..." << std::endl;
    run_hybrid_automaton_tests();
    std::cout << "Hybrid automata tests completed!" << std::endl;
    std::cout << "Running whitelist automata tests..." << std::endl;
    run_whitelist_automaton_tests();
    std::cout << "Automata tests completed!" << std::endl;

    std::cout << "Running Levenshtein generator tests..." << std::endl;
    run_levenshtein_generator_tests();
    std::cout << "Levenshtein generator tests completed!" << std::endl;

    compare_mutation_functions("TAAAACCCGGGTTCAG");
    
    return 0;
}