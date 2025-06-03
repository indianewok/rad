#pragma once
#include "rad_headers.h"
// set up a counter class to handle atomic operations
// Define an enum with unique indices
enum barcode_counts {
    raw = 0,          // Raw count
    forw = 1,         // Forward count
    forw_concat = 2,  // Forward concatenation count
    rev = 3,          // Reverse count
    rev_concat = 4,   // Reverse concatenation count
    total = 5,         // Total count
    corrected = 6,     // Corrected count
    filtered = 7       // filtered count
};

// building a thread-safe counter to be used during count updating processes
struct counter {
    // Use an array of 8 atomic ints.
    std::array<std::atomic<int>, 8> counts;

    // Default constructor: initialize all counters to 0.
    counter() : counts{{0, 0, 0, 0, 0, 0, 0, 0}} { }

    // A constructor that initializes all counters to a given value.
    counter(int init) : counts{{init, init, init, init, init, init, init, init}} { }

    // Copy constructor: copy each counter's current value.
    counter(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
    }

    // Assignment operator: copy each counter's value.
    counter& operator=(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
        return *this;
    }

    // Increment the counter at the specified index.
    void increment(int index) {
        counts[index].fetch_add(1, std::memory_order_relaxed);
    }

    // Overload: increment using one of the enum values for clarity.
    void increment(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(1, std::memory_order_relaxed);
    }

    // Decrement the counter at the specified slot
    void subtract(int index) {
        counts[index].fetch_add(-1, std::memory_order_relaxed);
    }
    
    // Overload: decrement using an enum value for clarity.
    void subtract(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(-1, std::memory_order_relaxed);
    }

    // Load the value of the counter at the specified index.
    int load(int index) const {
        return counts[index].load(std::memory_order_relaxed);
    }

    // Overload: load using an enum value for clarity.
    int load(barcode_counts index) const {
        return counts[static_cast<int>(index)].load(std::memory_order_relaxed);
    }

    // Load all counter values into a tuple.
    std::tuple<int, int, int, int, int, int, int, int> load_all() const {
        return std::make_tuple(
            counts[barcode_counts::raw].load(std::memory_order_relaxed),
            counts[barcode_counts::forw].load(std::memory_order_relaxed),
            counts[barcode_counts::forw_concat].load(std::memory_order_relaxed),
            counts[barcode_counts::rev].load(std::memory_order_relaxed),
            counts[barcode_counts::rev_concat].load(std::memory_order_relaxed),
            counts[barcode_counts::total].load(std::memory_order_relaxed),
            counts[barcode_counts::corrected].load(std::memory_order_relaxed),
            counts[barcode_counts::filtered].load(std::memory_order_relaxed)
        );
    }
};

// set up the int64_seq class to handle bit-to-sequence & sequence-to-bit
class int64_seq {
    public:
        uint16_t length;            // Total number of bases in the sequence.
        std::vector<int64_t> bits;  // Each int64_t encodes a chunk (2 bits per nucleotide, up to 32 bases per chunk).
        int64_seq() : length(0) {}

        explicit int64_seq(int64_t raw_bits, uint16_t seq_length) : length(seq_length) {
            bits.reserve(1);
            bits.push_back(raw_bits);
        }
    
        explicit int64_seq(const std::string& sequence) {
            sequence_to_bits(sequence);
        }
    
        void sequence_to_bits(const std::string& sequence) {
            length = static_cast<uint16_t>(sequence.size());
            const int chunk_size = 32; // Number of bases per chunk
            std::vector<int64_t> results;
            for (size_t i = 0; i < sequence.size(); i += chunk_size) {
              int64_t result = 0;
              size_t end = std::min(i + chunk_size, sequence.size());
              for (size_t j = i; j < end; ++j) {
                char c = sequence[j];
                result <<= 2;
                switch (c) {
                case 'A': break;
                case 'C': result |= 1; break;
                case 'T': result |= 2; break;
                case 'G': result |= 3; break;
                }
              }
              results.push_back(result);
            }
            bits = results;
          }
    
        std::string bits_to_sequence() const {
            int n = bits.size();
            std::string final_sequence;
            int chunk_size = 32;
            for (int i = 0; i < n; ++i) {
            int64_t int64_code = bits[i];
            std::string result;
            int bases_to_decode = (i < n - 1) ? chunk_size : length % chunk_size;
            bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
            for (int j = 0; j < bases_to_decode; ++j) {
                switch (int64_code & 3) {
                case 0: result.insert(result.begin(), 'A'); break;
                case 1: result.insert(result.begin(), 'C'); break;
                case 2: result.insert(result.begin(), 'T'); break;
                case 3: result.insert(result.begin(), 'G'); break;
                }
                int64_code >>= 2;
            }
            final_sequence += result;
            }
            return final_sequence;
        }
        
        bool operator==(const int64_seq &other) const {
            return (length == other.length) && (bits == other.bits);
        }
    
        void print_int64_seq() {
            std::cout << "Sequence length: " << length << "\n";
            std::cout << "Bit chunks: ";
            for (size_t i = 0; i < bits.size(); ++i) {
                std::cout << bits[i] << " ";
            }
            std::cout << "\n";
        }
    };

// class structure to hold barcode information and its associated count
struct barcode_entry {
    int64_seq barcode;         // barcode sequence
    mutable counter count;     // atomic counter for count
    int edit_dist;             // edit distance from a correct barcode
    bool filtered;             // flag for filtering
    std::string flags;        // additional flags
    
    barcode_entry()
      : barcode(""), 
       count(0),
       edit_dist(0), 
       filtered(false), 
       flags("")
       { }

    bool operator==(const barcode_entry &o) const noexcept {
        return barcode == o.barcode;
    }
};

// setting up hashing for int64_seq & barcode_entry
namespace std {
    template <>
    struct hash<int64_seq> {
        std::size_t operator()(const int64_seq &seq) const {
             // Skip length hashing
            return std::hash<int64_t>{}(seq.bits[0]);
        }
    };

    /* 
    //this version is slower, but takes bit length into account during hashing. will be best used for longer sequences
    template <>
    struct hash<int64_seq> {
        std::size_t operator()(const int64_seq &seq) const {
            std::size_t res = std::hash<uint16_t>()(seq.length);
            for (const auto &b : seq.bits) {
                boost::hash_combine(res, b);
            }
            return res;
        }
    };
    */

    template<>
    struct hash<barcode_entry> {
        size_t operator()(barcode_entry const& bc) const noexcept {
            // hash the barcode field
            return std::hash<int64_seq>()(bc.barcode);
        }
    };
}

namespace mutation_tools {
    // Myers's bit-parallel algorithm for edit distance, followed from the edlib comments for making a better version for shorter strings
    // and taking int64 values as input 
    // p:  bitmask where bit j=1 means pattern[j] matches '1'
    // t:  bitmask for text only used to select p or np per column
    // n:  number of columns (<=64)
    // max_dist: as soon as even the best possible remaining score exceeds max_dist, we return -1 as a sentinel.
    int bit_ld(int64_t p, int64_t t, int n, int max_dist) {
        int64_t np    = ~p;
        int64_t HMASK = 1LL << (n - 1);        // mask for the high bit of the window
        int64_t VP    = (1LL << n) - 1;        // all 1s
        int64_t VN    = 0;                     // all 0s
        int     score = n;                     // start at n
        for (int j = 0; j < n; ++j) {
            // select the Peq mask for column j
            int64_t Bj = ((t >> j) & 1) ? p : np;
            // Myers core
            int64_t D0 = (((VP + (Bj & VP)) ^ VP) | Bj | VN);
            int64_t HN = VP & D0;
            int64_t HP = VN | ~(VP | D0);
            // update score
            if (HP & HMASK) {
                score++;
            }
            else if (HN & HMASK) {
                score--;
            }
            // early cutoff bound: even if we subtract 1 in every remaining column, our best possible = score - (n - j - 1)
            int remaining = n - j - 1;
            if (score - remaining > max_dist)
                return -1;
            // advance VP, VN
            int64_t X = (HP << 1) | 1;
            VN = X & D0;
            VP = (HN << 1) | ~(X | D0);
        }
        return score;
    }

    // calculating levenshtein distance between one query and a set of multiple strings (naively, with no wildcard-based spacing)
    std::map<int, std::unordered_set<int64_seq>> int64_lvdist(const int64_seq &query, const std::unordered_set<int64_seq> &targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (query.bits.empty()){
            return results;
        }
        for (const auto &target : targets) {
            if (target.bits.empty()){
                continue;
            }
            int dist = bit_ld(query.bits[0], target.bits[0], query.length, max_dist);
            if(dist > 0){
                results[dist].insert(target);
            }
        }
        return results;
    }

    std::map<int, std::unordered_set<int64_seq>> int64_lvdist(const int64_seq &query, const std::vector<barcode_entry> &targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (query.bits.empty()){
            return results;
        }
        for (const auto &target : targets) {
            if (target.barcode.bits.empty()){
                continue;
            }
            int dist = bit_ld(query.bits[0], target.barcode.bits[0], query.length, max_dist);
            if(dist > 0){
                results[dist].insert(target.barcode);
            }
        }
        return results;
    }

    int int64_lvdist(const int64_seq &query, const int64_seq &target, int max_dist = 4) {
        int result = -1;
        if (query.bits.empty() || target.bits.empty()){
            return result;
        }
        result = bit_ld(query.bits[0], target.bits[0], query.length, max_dist);
        return result;
    }

    // generate point mutations for a given int64_seq sequence
    std::unordered_set<int64_seq> generate_point_mutations(const int64_seq &seq) {
        std::unordered_set<int64_seq> mutations;
        if (seq.bits.empty()){
            return mutations;
        }
        int sequence_length = seq.length;
        int64_t original_value = seq.bits[0];  // Assumes one chunk
        for (int pos = 0; pos < sequence_length; ++pos) {
            int64_t original_nuc = (original_value >> (2 * pos)) & 0b11;
            for (int64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                if (new_nuc != original_nuc) {
                    int64_t clear_mask = ~(0b11LL << (2 * pos));
                    int64_t mutated_value = (original_value & clear_mask) | (new_nuc << (2 * pos));
                    int64_seq candidate;
                    candidate.length = sequence_length;
                    candidate.bits.push_back(mutated_value);
                    mutations.insert(candidate);
                }
            }
        }
        return mutations;
    } 

    //raw pt mutation gen
    std::unordered_set<int64_t> generate_point_mutations_raw(int64_t original_value, int sequence_length) {
        std::unordered_set<int64_t> mutations;
        mutations.reserve(sequence_length * 3); // Pre-reserve for efficiency
        
        for (int pos = 0; pos < sequence_length; ++pos) {
            int64_t original_nuc = (original_value >> (2 * pos)) & 0b11;
            for (int64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                if (new_nuc != original_nuc) {
                    int64_t clear_mask = ~(0b11LL << (2 * pos));
                    int64_t mutated_value = (original_value & clear_mask) | (new_nuc << (2 * pos));
                    mutations.insert(mutated_value); // Fast int64_t hash and insert
                }
            }
        }
        return mutations;
    }

    // generate mutated barcode candidates across multiple rounds, returning a unique set
    std::unordered_set<int64_seq> generate_mutated_barcodes_depr(const int64_seq &seq, int mutation_rounds = 1) {
        if (seq.bits.empty()) return {};
        int sequence_length = seq.length;
        int64_seq initial = seq;
        //std::unordered_set<int64_seq> all_mutations{ 
        //    initial 
       // };
        //std::unordered_set<int64_seq> current_round{ 
        //    initial
        //};

            // *** HARD-CODED PRE-RESERVATION ***
        size_t expected_total;
        if (mutation_rounds == 1) {
            expected_total = 100;
        } else if (mutation_rounds == 2) {
            expected_total = 2000;
        } else if (mutation_rounds == 3) {
            expected_total = 60000;
        } else {
            expected_total = 1000; // fallback
        }
    
        // Pre-reserve all hash tables
        std::unordered_set<int64_seq> all_mutations;
        all_mutations.reserve(expected_total);
        all_mutations.insert(initial);
        
        std::unordered_set<int64_seq> current_round;
        current_round.reserve(mutation_rounds == 1 ? 100 : 2000);
        current_round.insert(initial);
    
        for (int round = 0; round < mutation_rounds; ++round) {
            std::unordered_set<int64_seq> next_round;
            // Hard-coded reservation for next round
            if (round == 0) {
                next_round.reserve(100);      // Round 1
            } else if (round == 1) {
                next_round.reserve(2000);     // Round 2  
            } else {
                next_round.reserve(60000);    // Round 3+
            }
        }
        // *** END HARD-CODED RESERVATION ***

        for (int round = 0; round < mutation_rounds; ++round) {
            std::unordered_set<int64_seq> next_round;
            for (const auto &candidate : current_round) {
                std::unordered_set<int64_seq> one_round = generate_point_mutations(candidate);
                for (const auto &mutation : one_round) {
                    if (all_mutations.insert(mutation).second) {
                        next_round.insert(mutation);
                    }
                }
            }
            if (next_round.empty()) break;
            current_round = std::move(next_round);
        }
        all_mutations.erase(initial);
        return all_mutations;
    }

    //generating mutated barcodes using raw int64_t values for speed
    std::unordered_set<int64_seq> generate_mutated_barcodes(const int64_seq &seq, int mutation_rounds = 1) {
        if (seq.bits.empty()) return {};
        
        int sequence_length = seq.length;
        int64_t initial_raw = seq.bits[0]; // Extract raw value once
        
        // *** HARD-CODED PRE-RESERVATION ***
        size_t expected_total;
        if (mutation_rounds == 1) {
            expected_total = 100;
        } else if (mutation_rounds == 2) {
            expected_total = 2000;
        } else if (mutation_rounds == 3) {
            expected_total = 60000;
        } else {
            expected_total = 1000; // fallback
        }
        
        // ALL GENERATION USING raw int64_t values
        std::unordered_set<int64_t> all_mutations_raw;
        all_mutations_raw.reserve(expected_total);
        all_mutations_raw.insert(initial_raw);
        
        std::unordered_set<int64_t> current_round_raw;
        current_round_raw.reserve(mutation_rounds == 1 ? 100 : 2000);
        current_round_raw.insert(initial_raw);
        
        // Fast raw mutation generation
        for (int round = 0; round < mutation_rounds; ++round) {
            std::unordered_set<int64_t> next_round_raw;
            
            // Hard-coded reservation for next round
            if (round == 0) {
                next_round_raw.reserve(100); // Round 1
            } else if (round == 1) {
                next_round_raw.reserve(2000); // Round 2  
            } else {
                next_round_raw.reserve(15000); // Round 3+
            }
            
            for (const auto candidate_raw : current_round_raw) {
                // Generate mutations directly as raw int64_t - no function call overhead
                for (int pos = 0; pos < sequence_length; ++pos) {
                    int64_t original_nuc = (candidate_raw >> (2 * pos)) & 0b11;
                    for (int64_t new_nuc = 0; new_nuc < 4; ++new_nuc) {
                        if (new_nuc != original_nuc) {
                            int64_t clear_mask = ~(0b11LL << (2 * pos));
                            int64_t mutated_value = (candidate_raw & clear_mask) | (new_nuc << (2 * pos));
                            
                            // Fast duplicate check and insert using raw int64_t
                            if (all_mutations_raw.insert(mutated_value).second) {
                                next_round_raw.insert(mutated_value);
                            }
                        }
                    }
                }
            }
            
            if (next_round_raw.empty()) break;
            current_round_raw = std::move(next_round_raw);
        }
        
        // Remove the original sequence from raw mutations
        all_mutations_raw.erase(initial_raw);
        
        // *** CONVERT TO int64_seq ONLY AT THE VERY END ***
        std::unordered_set<int64_seq> final_mutations;
        final_mutations.reserve(all_mutations_raw.size());
        
        // Single conversion pass - create int64_seq objects only here
        for (const int64_t raw_mutation : all_mutations_raw) {
            final_mutations.emplace(raw_mutation, sequence_length); // Use emplace for in-place construction
        }
        
        return final_mutations;
    }
    
    //generating shifted barcodes
    /*std::unordered_set<int64_seq> generate_shifted_barcodes(const int64_seq &seq, int shift) {
        std::unordered_set<int64_seq> result;
        // 1) Quick exits
        if (seq.bits.empty() || shift <= 0 || shift >= seq.length){
            return result;
        }
    
        // 2) Get the string
        std::string barcode = seq.bits_to_sequence();
        int L = (int)barcode.size();
    
        // 3) Compute cores
        std::string coreR = barcode.substr(0, L - shift);  // for left shifts
        std::string coreL = barcode.substr(shift);         // for right shifts
    
        const char bases[4] = {'A','C','G','T'};
        std::string buf(shift, ' ');
    
        // 4) Build left shifts
        std::function<void(int)> dfs_left = [&](int pos) {
            if (pos == shift) {
                // Make a fresh int64_seq from the new string
                int64_seq cand;
                cand.sequence_to_bits(buf + coreR);
                result.insert(std::move(cand));
                return;
            }
            for (char b : bases) {
                buf[pos] = b;
                dfs_left(pos + 1);
            }
        };
        dfs_left(0);
    
        // 5) Build right shifts
        std::function<void(int)> dfs_right = [&](int pos) {
            if (pos == shift) {
                int64_seq cand;
                cand.sequence_to_bits(coreL + buf);
                result.insert(std::move(cand));
                return;
            }
            for (char b : bases) {
                buf[pos] = b;
                dfs_right(pos + 1);
            }
        };
        dfs_right(0);
    
        // 6) Remove the original if it was regenerated
        result.erase(seq);
        return result;
    }*/

    // Call this exactly the same way as before, but it will reuse its buffer
    std::unordered_set<int64_seq> generate_shifted_barcodes(const int64_seq &seq, int shift) {
        // 1) Keep one thread-local set so we only pay the bucket allocations once
        static thread_local std::unordered_set<int64_seq> result;
        result.clear();

        // 2) Fast exit if nothing to do
        if (seq.bits.empty() || shift <= 0 || shift >= seq.length) {
            return result;
        }

        // 3) Pre-compute how many combinations we'll get:
        size_t combos = 1ULL << (2 * shift);        // 4^shift
        result.reserve(combos * 2);                 // both left + right shifts

        // 4) Extract the original 2-bit code
        uint64_t orig = seq.bits[0];
        int L = seq.length;
        int core_len = L - shift;

        // 5) Split into “core” and “suffix” bits
        //    - coreR is the top (L-shift) bases
        //    - suffix is the bottom (shift) bases
        uint64_t coreR  = orig >> (2 * shift);
        uint64_t coreL  = orig >> (2 * shift);
        uint64_t suffix = orig & ((1ULL << (2 * shift)) - 1);

            // 6) Build all left shifts: [ new_prefix | coreR ]
            for (uint64_t code = 0; code < combos; ++code) {
                int64_seq cand;
                cand.length = L;
                cand.bits   = {
                    static_cast<int64_t>(  // <-- here
                    (code << (2 * core_len))
                    | coreR
                    )
                };
                result.insert(std::move(cand));
            }
        
            // 7) Build all right shifts: [ coreL | new_suffix ]
            for (uint64_t code = 0; code < combos; ++code) {
                int64_seq cand;
                cand.length = L;
                cand.bits   = {
                    static_cast<int64_t>(
                    (coreL << (2 * shift))
                    | code
                    )
                };
                result.insert(std::move(cand));
            }

        // 8) Drop the original if it snuck back in
        result.erase(seq);
        return result;
    }

};

template<typename key, typename value>
struct bc_multimap {
private:
    // Helper to extract key from either a key or a value (unchanged)
    static inline const key& key_of(key const &k) noexcept { 
        return k;
    }
    static inline const key& key_of(value const &v) noexcept { 
        return v.barcode; 
    }
public:
    // Store unique values - using a set lets us find and reuse identical values
    std::unordered_set<value> unique_values;
    // Maps keys to pointers to values in the unique set
    std::unordered_multimap<key, const value*> associations;

    // Size reporting
    size_t size() const { 
        return associations.size(); 
    }

    bool empty() const { 
        return associations.empty(); 
    }

    void clear() { 
        associations.clear();
        unique_values.clear(); 
    }

    void clear_associations() {
        associations.clear();
    }

    const auto& debug_unique_values() const noexcept { 
        return unique_values; 
    }

    const auto& debug_associations() const noexcept { 
        return associations; 
    }

    // ==== Core Insertion Methods ====
    // Insert a barcode entry with an existing value
    template<typename T>
    void insert_bc_entry(const T &observed, const value &correct) {
        // First, insert or find the value in unique_values
        auto [value_it, value_inserted] = unique_values.insert(correct);
        const value* ptr = &(*value_it);
        // Now map the key to the unique value
        associations.emplace(key_of(observed), ptr);
    }
    // Insert with a moved value
    template<typename T>
    void insert_bc_entry(const T &observed, value &&correct) {
        // When moving, we need to handle things a bit differently
        // First check if an equivalent value already exists
        auto find_it = unique_values.find(correct);
        
        if (find_it != unique_values.end()) {
            // Use existing value
            associations.emplace(key_of(observed), &(*find_it));
        } else {
            // Insert the moved value
            auto [value_it, inserted] = unique_values.insert(std::move(correct));
            associations.emplace(key_of(observed), &(*value_it));
        }
    }
    // Use the key as value too
    template<typename T>
    void insert_bc_entry(const T &observed) {
        value v{};
        v.barcode = key_of(observed);
        insert_bc_entry(observed, std::move(v));
    }
    
    // In-place construction
    template<typename T, typename... Args>
    void emplace_bc_entry(const T &observed, Args&&... args) {
        // Construct the value with the provided arguments
        value v{std::forward<Args>(args)...};
        insert_bc_entry(observed, std::move(v));
    }
    
    // Direct emplace method for compatibility with the original implementation
    template<typename K, typename V>
    void emplace(K&& k, V&& v) {
        // First, insert or find the value in unique_values
        // Use perfect forwarding to preserve move semantics if possible
        auto [value_it, value_inserted] = unique_values.insert(std::forward<V>(v));
        const value* ptr = &(*value_it);
        
        // Now map the key to the unique value
        associations.emplace(std::forward<K>(k), ptr);
    }
    
    // ==== Removal Methods ====
    // Remove a key-value entry
    template<typename T>
    void remove_bc_entry(const T &observed) {
        associations.erase(key_of(observed));
        // A periodic cleanup could remove orphaned values if needed
    }
    
    // ==== Lookup Methods ====
    // Check if whitelist contains this key
    template<typename T>
    bool check_wl_for(T const &x) const {
        if constexpr (std::is_same_v<T,key> || std::is_same_v<T,value>) {
            return associations.find(key_of(x)) != associations.end();
        } else {
            // container of keys
            for (auto const &y : x) {
                if (check_wl_for(y)) return true;
            }
            return false;
        }
    }
    
    // Return matching barcodes
    template<typename T>
    std::unordered_set<key> return_matching_barcodes(T const &x) const {
        std::unordered_set<key> out;
        if constexpr (std::is_same_v<T,key> || std::is_same_v<T,value>) {
            if (check_wl_for(x))
                out.insert(key_of(x));
        } else {
            for (auto const &y : x) {
                if (check_wl_for(y))
                    out.insert(key_of(y));
            }
        }
        return out;
    }
    
    // Return putative correct barcodes
    template<typename T>
    std::unordered_set<key> return_putative_correct_bcs(T const &x) const {
        std::unordered_set<key> out;
        if constexpr (std::is_same_v<T,key> || std::is_same_v<T,value>) {
            auto range = associations.equal_range(key_of(x));
            for (auto it = range.first; it != range.second; ++it)
                out.insert(it->second->barcode);
        } else {
            for (auto const &y : x) {
                auto sub = return_putative_correct_bcs(y);
                out.insert(sub.begin(), sub.end());
            }
        }
        return out;
    }
    
    // ==== Counter Methods ====
    
    // Add a count to a barcode_entry
    template<typename T>
    void update_bc_count(T const &x, barcode_counts slot = total) const {
        auto range = associations.equal_range(key_of(x));
        for (auto it = range.first; it != range.second; ++it) {
            const_cast<value*>(it->second)->count.increment(slot);
        }
    }
    // Subtract a count
    template<typename T>
    void subtract_bc_count(T const &x, barcode_counts slot = total) const {
        auto range = associations.equal_range(key_of(x));
        for (auto it = range.first; it != range.second; ++it) {
            const_cast<value*>(it->second)->count.subtract(slot);
        }
    }
    // Get the count for a barcode entry
    template<typename T>
    int get_bc_count(T const &x, barcode_counts slot = total) const {
        auto range = associations.equal_range(key_of(x));
        if (range.first != range.second)
            return range.first->second->count.load(slot);
        return 0;
    }
    
    // ==== Summary Methods ====
    
    // Summarize counts
    std::vector<std::tuple<std::string,int,int,int,int,int,int,int,int>>
    summarize_counts(const std::unordered_set<key> *filter_keys = nullptr) const {
        using row = std::tuple<std::string,int,int,int,int,int,int,int,int>;
        std::vector<row> rows;
        
        // Build the list of keys to process
        std::vector<key> keys;
        if (filter_keys) {
            keys.reserve(filter_keys->size());
            for (auto const &k : *filter_keys)
                keys.push_back(k);
        } else {
            // Collect distinct keys from the multimap
            std::unordered_set<key> seen;
            seen.reserve(associations.size());
            for (auto const &kv : associations)
                seen.insert(kv.first);
            keys.reserve(seen.size());
            for (auto const &k : seen)
                keys.push_back(k);
        }
        
        rows.reserve(keys.size());
        
        // For each key, find its perfect‐match mapping and pull the counts
        for (auto const &bc : keys) {
            auto range = associations.equal_range(bc);
            int raw=0, tot=0, corr=0, fwd=0, rev=0, fwd_c=0, rev_c=0, filt=0;
            
            if (range.first != range.second) {
                // Prefer the exact key if present
                auto it = range.first;
                for (; it != range.second; ++it) {
                    if (it->first == bc)
                        break;
                }
                if (it == range.second)
                    it = range.first;
                
                raw = it->second->count.load(barcode_counts::raw);
                tot = it->second->count.load(barcode_counts::total);
                corr = it->second->count.load(barcode_counts::corrected);
                fwd = it->second->count.load(barcode_counts::forw);
                rev = it->second->count.load(barcode_counts::rev);
                fwd_c = it->second->count.load(barcode_counts::forw_concat);
                rev_c = it->second->count.load(barcode_counts::rev_concat);
                filt = it->second->count.load(barcode_counts::filtered);
            }
            rows.emplace_back(bc.bits_to_sequence(), raw, tot, corr, fwd, rev, fwd_c, rev_c, filt);
        }
        return rows;
    }
    
    // Write summary to stream
    void write_wl_summary(std::ostream &out, const std::string class_id, 
                        const std::unordered_set<key> *filter_keys = nullptr) const {
        static bool header_printed = false;
        if(!header_printed) {
            out << "class_id,true_barcode,raw_count,total_count,corrected_count,forw_count,rev_count,"
                << "forw_concat_count,rev_concat_count,filtered_count\n";
            header_printed = true;
        }
        
        auto rows = summarize_counts(filter_keys);
        for (auto const & tpl : rows) {
            auto const & bc = std::get<0>(tpl);
            auto raw = std::get<1>(tpl);
            auto tot  = std::get<2>(tpl);
            auto corr = std::get<3>(tpl);
            auto forw = std::get<4>(tpl);
            auto rev  = std::get<5>(tpl);
            auto forw_concat = std::get<6>(tpl);
            auto rev_concat = std::get<7>(tpl);
            auto filtered = std::get<8>(tpl);
            
            // Print the row
            out << class_id << ',' << bc << ','  << raw << "," << tot << ',' << corr << ',' << forw 
                << ',' << rev << ',' << forw_concat << ',' << rev_concat << ',' << filtered << '\n';
        }
    }
    
    // Write summary to file
    void write_wl_summary(const std::string &path = "", 
                        const std::unordered_set<key> *filter_keys = nullptr) const {
        std::ofstream ofs(path);
        if (!ofs) {
            throw std::runtime_error("Failed to open output file: " + path);
        }
        write_wl_summary(ofs, filter_keys);
    }
    
    // ==== Iterator Support ====
    
    // Iterator adapter to dereference pointers automatically
    template<typename BaseIterator>
    class deref_iterator {
    private:
        BaseIterator it;
        
    public:
        // Iterator type definitions
        using difference_type = typename BaseIterator::difference_type;
        using value_type = std::pair<const key&, const value&>;
        using reference = value_type;
        using pointer = void;
        using iterator_category = std::forward_iterator_tag;
        
        deref_iterator(BaseIterator iter) : it(iter) {}
        
        deref_iterator& operator++() { ++it; return *this; }
        deref_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
        
        bool operator==(const deref_iterator& other) const { return it == other.it; }
        bool operator!=(const deref_iterator& other) const { return it != other.it; }
        
        reference operator*() const { 
            return {it->first, *(it->second)}; 
        }
    };
    
    // Iterator type definitions
    using iterator = deref_iterator<typename std::unordered_multimap<key, const value*>::iterator>;
    using const_iterator = deref_iterator<typename std::unordered_multimap<key, const value*>::const_iterator>;
    
    // Iterator methods
    iterator begin() { return iterator(associations.begin()); }
    iterator end() { return iterator(associations.end()); }
    const_iterator begin() const { return const_iterator(associations.begin()); }
    const_iterator end() const { return const_iterator(associations.end()); }
    const_iterator cbegin() const { return const_iterator(associations.cbegin()); }
    const_iterator cend() const { return const_iterator(associations.cend()); }
    
    // Equal range returning dereferencing iterators
    std::pair<iterator, iterator> equal_range(const key& k) {
        auto range = associations.equal_range(k);
        return {iterator(range.first), iterator(range.second)};
    }
    
    std::pair<const_iterator, const_iterator> equal_range(const key& k) const {
        auto range = associations.equal_range(k);
        return {const_iterator(range.first), const_iterator(range.second)};
    }
    
    // Find method with dereferencing iterator
    iterator find(const key& k) {
        return iterator(associations.find(k));
    }
    
    const_iterator find(const key& k) const {
        return const_iterator(associations.find(k));
    }
};

class whitelist {
    public:
        whitelist() = default;
        ~whitelist() = default;

    struct wl_entry {
        std::unordered_set<int64_seq> true_ref;
        bc_multimap<int64_seq, barcode_entry> true_bcs, global_bcs, filter_bcs;
    
        void generate_mismatch_barcodes_(int shift, int mutation_rounds, bool verbose, int nthreads = 1) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for "
                        << true_ref.size() << " true barcodes with shift=" << shift
                        << " and mutation_rounds=" << mutation_rounds << "\n";
                    std::cout << oss.str();
                }
            }

            // build a single hash‐map
            std::unordered_map<int64_seq, barcode_entry> originals;
            originals.reserve(true_bcs.size());
            for (auto const &kv : true_bcs) {
                originals.emplace(kv.first, kv.second);
            }
            //generate shifts & mutations
            for (auto const & [orig_bits, orig_be] : originals) {
                // SHIFTED
                auto shifted = mutation_tools::generate_shifted_barcodes(orig_bits, shift);
                for (auto const &s_bits : shifted) {
                    // skip if it’s already a true barcode or if it equals one of the originals
                    if (global_bcs.check_wl_for(s_bits) || originals.count(s_bits)){
                        continue;
                    }
                    barcode_entry obs = orig_be;
                    obs.barcode = s_bits;
                    true_bcs.insert_bc_entry(obs, orig_be);
                }

                // MUTATED
                auto mutated = mutation_tools::generate_mutated_barcodes(orig_bits, mutation_rounds);
                for (auto const &m_bits : mutated) {
                    if (global_bcs.check_wl_for(m_bits) || originals.count(m_bits)){
                        continue;
                    }
                    barcode_entry obs = orig_be;
                    obs.barcode = m_bits;
                    true_bcs.insert_bc_entry(obs, orig_be);
                }
            }

            size_t missing = 0;
            for (auto const & [orig_bits, orig_be] : originals) {
                if (! true_bcs.check_wl_for(orig_bits)) {
                    ++missing;
                    std::cerr << "[ERROR] Original barcode vanished: "
                              << orig_bits.bits_to_sequence() << "\n";
                }
            }
            if (verbose) {
                if (missing == 0) {
                    std::cout << "[generate_mismatch_barcodes] All "
                              << originals.size()
                              << " original barcodes are present\n";
                } else {
                    std::cerr << "[generate_mismatch_barcodes] "
                              << missing << "/"
                              << originals.size()
                              << " originals missing!\n";
                }
            }
        }

        // Generate shifted and mutated barcodes for the true_ref set
        // also implicitly filters mutations that would be generated against the global bcs set
        // so all barcodes and associations in true bcs are guaranteed to be unique
        void generate_mismatch_barcodes__(int shift, int mutation_rounds, bool verbose, int nthreads = 1) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for "
                        << true_ref.size() << " true barcodes with shift = " << shift
                        << " and mutation_rounds = " << mutation_rounds << "\n";
                    std::cout << oss.str();
                }
            }

            // 1) Collect each original barcode_entry* exactly once
            std::vector<const barcode_entry*> originals;
            originals.reserve(true_bcs.size());
            {
                std::unordered_set<const barcode_entry*> seen;
                seen.reserve(true_bcs.size());
                for (auto const &p : true_bcs.associations) {
                    if (seen.insert(p.second).second)
                    originals.push_back(p.second);
                }
            }
            // 2) For each original, generate shifts & mutations
            for (auto const *orig_be : originals) {
                const auto &orig_bits = orig_be->barcode;

            // SHIFTED
            auto shifted = mutation_tools::generate_shifted_barcodes(orig_bits, shift);
            for (auto const &s_bits : shifted) {
                if (!global_bcs.check_wl_for(s_bits)) {
                // insert without copying orig_be
                    true_bcs.insert_bc_entry(s_bits, *orig_be);
                }
            }

            // MUTATED
            auto mutated = mutation_tools::generate_mutated_barcodes(orig_bits, mutation_rounds);
                for (auto const &m_bits : mutated) {
                    if (!global_bcs.check_wl_for(m_bits)) {
                        true_bcs.insert_bc_entry(m_bits, *orig_be);
                    }
                }
            }

            // 3) Optional sanity check (using true_ref instead of originals)
            size_t missing = 0;
            for (auto const &orig_bits : true_ref) {
                if (!true_bcs.check_wl_for(orig_bits)) {
                    ++missing;
                    std::cerr << "[ERROR] Original barcode vanished: "
                            << orig_bits.bits_to_sequence() << "\n";
                }
            }

            if (verbose) {
                if (missing == 0) {
                    std::cout << "[generate_mismatch_barcodes] All "
                            << true_ref.size()
                            << " original barcodes are present\n";
                } else {
                    std::cerr << "[generate_mismatch_barcodes] "
                            << missing << "/"
                            << true_ref.size()
                            << " originals missing!\n";
                }
            }
        }
    
    void generate_mismatch_barcodes(int shift, int mutation_rounds, bool verbose, int nthreads = 1) {
        using namespace std::chrono;
        auto start_total = high_resolution_clock::now();
        
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for "
                    << true_ref.size() << " true barcodes with shift = " << shift
                    << " and mutation_rounds = " << mutation_rounds << "\n";
                std::cout << oss.str();
            }
        }
        
        // TIMING: Collect original barcode entries
        auto start_collect = high_resolution_clock::now();
        
        std::vector<const barcode_entry*> originals;
        originals.reserve(true_bcs.size());
        {
            std::unordered_set<const barcode_entry*> seen;
            seen.reserve(true_bcs.size());
            for (auto const &p : true_bcs.associations) {
                if (seen.insert(p.second).second)
                    originals.push_back(p.second);
            }
        }
        
        auto end_collect = high_resolution_clock::now();
        double collect_ms = duration<double, std::milli>(end_collect - start_collect).count();
        
        if (verbose) {
            std::cout << "[generate_mismatch_barcodes] Collected " << originals.size() 
                    << " unique originals in " << collect_ms << " ms\n";
        }
        
        // Memory tracking
        size_t initial_true_bcs_size = true_bcs.associations.size();
        size_t initial_memory_estimate = initial_true_bcs_size * 48; // Rough estimate: 48 bytes per association
        
        if (verbose) {
            std::cout << "[generate_mismatch_barcodes] Initial true_bcs size: " << initial_true_bcs_size 
                    << " associations (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
        }
        
        // TIMING: Generate mutations and shifts
        auto start_generation = high_resolution_clock::now();
        
        size_t total_shifts_generated = 0;
        size_t total_mutations_generated = 0;
        size_t shifts_added = 0;
        size_t mutations_added = 0;
        size_t shifts_rejected_global = 0;
        size_t mutations_rejected_global = 0;
        
        // Track timing for individual operations
        double total_shift_generation_ms = 0;
        double total_mutation_generation_ms = 0;
        double total_shift_lookup_ms = 0;
        double total_mutation_lookup_ms = 0;
        double total_shift_insert_ms = 0;
        double total_mutation_insert_ms = 0;
        
        auto start_loop = high_resolution_clock::now();
        
        for (size_t i = 0; i < originals.size(); ++i) {
            auto const *orig_be = originals[i];
            const auto &orig_bits = orig_be->barcode;
            
            // Time individual barcode processing for first few
            auto barcode_start = high_resolution_clock::now();
            
            // === SHIFTED GENERATION ===
            auto shift_gen_start = high_resolution_clock::now();
            auto shifted = mutation_tools::generate_shifted_barcodes(orig_bits, shift);
            auto shift_gen_end = high_resolution_clock::now();
            
            double shift_gen_time = duration<double, std::milli>(shift_gen_end - shift_gen_start).count();
            total_shift_generation_ms += shift_gen_time;
            total_shifts_generated += shifted.size();
            
            // === SHIFTED LOOKUP & INSERT ===
            auto shift_lookup_start = high_resolution_clock::now();
            
            for (auto const &s_bits : shifted) {
                bool in_wl = global_bcs.check_wl_for(s_bits) || true_bcs.check_wl_for(s_bits);
                if (!in_wl) {
                    auto shift_insert_start = high_resolution_clock::now();
                    true_bcs.insert_bc_entry(s_bits, *orig_be);
                    auto shift_insert_end = high_resolution_clock::now();
                    total_shift_insert_ms += duration<double, std::milli>(shift_insert_end - shift_insert_start).count();
                    shifts_added++;
                } else {
                    shifts_rejected_global++;
                }
            }
            
            auto shift_lookup_end = high_resolution_clock::now();
            total_shift_lookup_ms += duration<double, std::milli>(shift_lookup_end - shift_lookup_start).count();
            
            // === MUTATION GENERATION ===
            auto mutation_gen_start = high_resolution_clock::now();
            auto mutated = mutation_tools::generate_mutated_barcodes(orig_bits, mutation_rounds);
            auto mutation_gen_end = high_resolution_clock::now();
            
            double mutation_gen_time = duration<double, std::milli>(mutation_gen_end - mutation_gen_start).count();
            total_mutation_generation_ms += mutation_gen_time;
            total_mutations_generated += mutated.size();
            
            // === MUTATION LOOKUP & INSERT ===
            auto mutation_lookup_start = high_resolution_clock::now();
            for (auto const &m_bits : mutated) {
                auto lookup_start = high_resolution_clock::now();
                bool in_wl = global_bcs.check_wl_for(m_bits) || true_bcs.check_wl_for(m_bits);
                auto lookup_end = high_resolution_clock::now();
                total_mutation_lookup_ms += duration<double, std::milli>(lookup_end - lookup_start).count();
                if (!in_wl) {
                    auto insert_start = high_resolution_clock::now();
                    true_bcs.insert_bc_entry(m_bits, *orig_be);
                    auto insert_end = high_resolution_clock::now();
                    total_mutation_insert_ms += duration<double, std::milli>(insert_end - insert_start).count();
                    mutations_added++;
                } else {
                    mutations_rejected_global++;
                }
            }
            
            auto mutation_lookup_end = high_resolution_clock::now();
            // Detailed timing for first 10 barcodes
            if (verbose && i < 10) {
                auto barcode_end = high_resolution_clock::now();
                double total_time = duration<double, std::milli>(barcode_end - barcode_start).count();
                std::cout << "[generate_mismatch_barcodes] Barcode " << i + 1 << " (" 
                        << orig_bits.bits_to_sequence() << "): "
                        << shifted.size() << " shifts (" << shift_gen_time << "ms gen), "
                        << mutated.size() << " mutations (" << mutation_gen_time << "ms gen), "
                        << "total " << total_time << "ms\n";
            }
            
            // Memory usage check every 100000 barcodes
            if (verbose && i % 1000 == 0 && i > 0) {
                size_t current_true_bcs_size = true_bcs.associations.size();
                size_t growth = current_true_bcs_size - initial_true_bcs_size;
                size_t current_memory_estimate = current_true_bcs_size * 48;
                
                std::cout << "[generate_mismatch_barcodes] After " << i << " barcodes: "
                        << current_true_bcs_size << " associations (+" << growth << "), "
                        << "~" << (current_memory_estimate / 1024 / 1024) << " MB total\n";
            }
        }
        
        auto end_loop = high_resolution_clock::now();
        double loop_ms = duration<double, std::milli>(end_loop - start_loop).count();
        
        auto end_generation = high_resolution_clock::now();
        double generation_ms = duration<double, std::milli>(end_generation - start_generation).count();
        
        // Final memory usage
        size_t final_true_bcs_size = true_bcs.associations.size();
        size_t total_growth = final_true_bcs_size - initial_true_bcs_size;
        size_t final_memory_estimate = final_true_bcs_size * 48;
        
        // TIMING: Sanity check
        auto start_sanity = high_resolution_clock::now();
        
        size_t missing = 0;
        for (auto const &orig_bits : true_ref) {
            if (!true_bcs.check_wl_for(orig_bits)) {
                ++missing;
                std::cerr << "[ERROR] Original barcode vanished: "
                        << orig_bits.bits_to_sequence() << "\n";
            }
        }
        
        auto end_sanity = high_resolution_clock::now();
        double sanity_ms = duration<double, std::milli>(end_sanity - start_sanity).count();
        
        auto end_total = high_resolution_clock::now();
        double total_ms = duration<double, std::milli>(end_total - start_total).count();
        
        if (verbose) {
            std::cout << "\n=== MISMATCH GENERATION TIMING REPORT ===\n";
            std::cout << "Total time: " << total_ms << " ms (" << (total_ms/1000.0) << " s)\n";
            std::cout << "  - Collection phase: " << collect_ms << " ms (" << (collect_ms/total_ms*100) << "%)\n";
            std::cout << "  - Generation loop: " << loop_ms << " ms (" << (loop_ms/total_ms*100) << "%)\n";
            std::cout << "  - Sanity check: " << sanity_ms << " ms (" << (sanity_ms/total_ms*100) << "%)\n";
            
            std::cout << "\n=== DETAILED OPERATION BREAKDOWN ===\n";
            std::cout << "Shift generation: " << total_shift_generation_ms << " ms (" 
                    << (total_shift_generation_ms/total_ms*100) << "%)\n";
            std::cout << "Mutation generation: " << total_mutation_generation_ms << " ms (" 
                    << (total_mutation_generation_ms/total_ms*100) << "%)\n";
            std::cout << "Shift lookups: " << total_shift_lookup_ms << " ms (" 
                    << (total_shift_lookup_ms/total_ms*100) << "%)\n";
            std::cout << "Mutation lookups: " << total_mutation_lookup_ms << " ms (" 
                    << (total_mutation_lookup_ms/total_ms*100) << "%)\n";
            std::cout << "Shift insertions: " << total_shift_insert_ms << " ms (" 
                    << (total_shift_insert_ms/total_ms*100) << "%)\n";
            std::cout << "Mutation insertions: " << total_mutation_insert_ms << " ms (" 
                    << (total_mutation_insert_ms/total_ms*100) << "%)\n";
            
            std::cout << "\n=== GENERATION STATISTICS ===\n";
            std::cout << "Processed " << originals.size() << " original barcodes\n";
            std::cout << "Shifts: generated " << total_shifts_generated << ", added " << shifts_added 
                    << ", rejected " << shifts_rejected_global << "\n";
            std::cout << "Mutations: generated " << total_mutations_generated << ", added " << mutations_added 
                    << ", rejected " << mutations_rejected_global << "\n";
            std::cout << "Average shifts per barcode: " << (double)total_shifts_generated / originals.size() << "\n";
            std::cout << "Average mutations per barcode: " << (double)total_mutations_generated / originals.size() << "\n";
            std::cout << "Time per barcode: " << (loop_ms / originals.size()) << " ms\n";
            
            std::cout << "\n=== MEMORY USAGE ===\n";
            std::cout << "Initial associations: " << initial_true_bcs_size 
                    << " (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
            std::cout << "Final associations: " << final_true_bcs_size 
                    << " (~" << (final_memory_estimate / 1024 / 1024) << " MB)\n";
            std::cout << "Growth: +" << total_growth << " associations (+~" 
                    << ((final_memory_estimate - initial_memory_estimate) / 1024 / 1024) << " MB)\n";
            std::cout << "Memory efficiency: " << (double)total_growth / (total_shifts_generated + total_mutations_generated) * 100 
                    << "% of generated items were unique and added\n";
            
            if (missing == 0) {
                std::cout << "\n[generate_mismatch_barcodes] All " << true_ref.size() 
                        << " original barcodes are present\n";
            } else {
                std::cerr << "\n[generate_mismatch_barcodes] " << missing << "/" 
                        << true_ref.size() << " originals missing!\n";
            }
        }
    }
    };

    std::unordered_map<std::string, std::reference_wrapper<whitelist::wl_entry>> maps;
    std::unordered_map<std::string, whitelist::wl_entry> lists;
    // operator[] for easy insert/access: W["barcode_1"].insert_bc_entry(obs, corr);
    wl_entry & operator[](const std::string &class_id) {
        return maps.at(class_id).get();
    }

    //check all loaded class_ids in an entry
    std::vector<std::string> class_ids() const {
        std::vector<std::string> out;
        out.reserve(maps.size());
        for (auto const &kv : maps)
            out.push_back(kv.first);
        return out;
    }
    
      // for one spec (kit name or file path), load its int64_seq set
    static std::unordered_set<int64_seq> load_barcodes(std::string const &spec, uint16_t default_length, bool verbose) {
        if (verbose) std::cout << "[load_barcodes] Loading barcodes from path/kit: " << spec << "\n";
        // check if spec is a kit name or a file path
        if (verbose) {
            if(whitelist_utils::is_kit(spec)){
                std::cout << "[load_barcodes] Kit name: " << spec << "\n";
                std::cout << "[load_barcodes] Kit path: " << whitelist_utils::kit_to_path(spec) << "\n";
            } else {
                std::cout << "[load_barcodes] File path: " << spec << "\n";
            }
        }
        std::string path = whitelist_utils::kit_to_path(spec);
        // detect bitlist vs char list
        bool isBitlist = whitelist_utils::check_if_bitlist(spec, verbose, /*N=*/10);
        if(verbose){
            std::string islist = isBitlist ? "IS" : "IS NOT";
            std::cout << "[check_if_bitlist] Detected element " << islist << " a bit-configured whitelist!" << std::endl;
        }

        // read all lines (skipping header)  
        auto lines = streaming_utils::import_text(path, SIZE_MAX);
        if (!lines.empty()){
            lines.erase(lines.begin());
        }
        std::unordered_set<int64_seq> out;
        for (auto &ln : lines) {
            if (ln.empty()) continue;
            auto p = ln.find_first_of(",\t");
            std::string tok = (p == std::string::npos ? ln : ln.substr(0,p));

            int64_seq seq;
            if (isBitlist) {
                try {
                    //if(verbose) std::cout << "[load_barcodes] Loading bitlist...\n";
                    int64_t code = std::stoll(tok);
                    seq.length = default_length;
                    seq.bits   = { code };
                } catch (...) { continue; }
            } else {
                //if(verbose) std::cout << "[load_barcodes] Converting characters to bits...\n";
                seq.sequence_to_bits(tok);
                if (seq.bits.empty()) continue;
            }
            out.insert(seq);
        }
        return out;
    }

    // populate a wl_entry with the given barcode set
    void populate_entry(wl_entry &out, std::vector<std::unordered_set<int64_seq>> const &sets){
        if (sets.size() == 1) {
            auto const &A = sets[0];
            for (auto const &seq : A) {
                barcode_entry be;
                be.barcode = seq;
                be.edit_dist = 0;
                be.filtered = false;
                be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
                out.true_ref.insert(seq);
            }
            return;
        }
        // two lists
        auto const &A = sets[0];
        auto const &B = sets[1];

        if ((A.size() != B.size())) {
            // pick smaller -> true, larger -> global
            auto const &small = (A.size() < B.size() ? A : B);
            auto const &large = (A.size() < B.size() ? B : A);

            for (auto const &seq : small) {
                barcode_entry be;
                be.barcode = seq;
                be.edit_dist = 0;
                be.filtered = false;
                be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
                out.true_ref.insert(seq);
            }
            for (auto const &seq : large) {
                if (out.true_ref.count(seq) == 0) {
                    barcode_entry be{};
                    be.barcode   = seq;
                    be.edit_dist = 0;
                    be.filtered  = false;
                    be.flags     = "";
                    out.global_bcs.emplace(seq, be);
                }
            }
        }
        else {
            // fallback: merge both into true_bcs
            for (auto const &seq : A) {
                barcode_entry be;
                be.barcode = seq;
                be.edit_dist = 0;
                be.filtered = false;
                be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
                out.true_ref.insert(seq);
            }
            for (auto const &seq : B) {
                barcode_entry be;
                be.barcode = seq;
                be.edit_dist = 0;
                be.filtered = false;
                be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
                out.true_ref.insert(seq);
            }
        }
    }
    
    //import a whitelist from file--can be true barcodes or not. cheap but easy way to import b/w both global or custom--just set
    //an arbitrary threshold for the number of true barcodes to be kept because ideally after a certain size of barcodes
    //you really just want to treat them both the same way
    wl_entry import_whitelist(std::string const &field, bool verbose, uint16_t default_length = 16) {
        if(verbose) std::cout << "[import_whitelist] Importing whitelist from " << field << "\n";
        // split into 1 or 2 specs
        auto specs = whitelist_utils::parse_whitelist_specs(field);
        if (specs.empty()) {
            std::cerr << "[warning][import_whitelist] empty whitelist spec—returning empty entry\n";
            return {};
        } else {
            for (auto const &spec : specs) {
                if(verbose) std::cout << "[import_whitelist] whitelist kit or path:" << spec << "\n";
            }
        }

        std::vector<std::unordered_set<int64_seq>> sets;
        for (auto const &spec : specs) {
            sets.push_back(load_barcodes(spec, default_length, verbose));
        }

        // assemble the wl_entry
        wl_entry out;
        out.global_bcs.clear();
        out.true_bcs.clear();
        populate_entry(out, sets);
        return out;
    }

};

//needed to include this down here because of the definitions above that I couldn't make work by movign this to the misc utils, though I sure would love for this to be there
namespace bc_mem_utils {
    /// Exact byte‐size of barcode_entry type
    inline constexpr std::size_t get_bc_entry_mem() {
        return sizeof(barcode_entry);
    }
    
    /// Exact byte‐size of int64_seq key
    inline constexpr std::size_t get_int64seq_mem() {
        return sizeof(int64_seq);
    }

    template<typename Key, typename Value>
    std::size_t approx_unique(const bc_multimap<Key,Value>& m) {
        auto const &s = m.unique_values;
        // buckets array
        std::size_t buckets = s.bucket_count() * memory_utils::get_pointer_mem();
        // each node: stored Value + two pointers (next/in‐bucket links)
        std::size_t nodes   = s.size() * (sizeof(Value) + 2*memory_utils::get_pointer_mem());
        return buckets + nodes;
    }

    /// Approximate footprint of the key→pointer associations
    template<typename Key, typename Value>
    std::size_t approx_assoc(const bc_multimap<Key,Value>& m) {
        auto const &mm = m.associations;
        // buckets array
        std::size_t buckets = mm.bucket_count() * memory_utils::get_pointer_mem();
        // each node: Key + Value* + two pointers
        std::size_t nodes   = mm.size() * (
            sizeof(Key)
            + sizeof(const Value*)
            + 2*memory_utils::get_pointer_mem()
        );
        return buckets + nodes;
    }

    /// Total (unique_values + associations)
    template<typename Key, typename Value>
    std::size_t get_full_wl_mem(const bc_multimap<Key,Value>& m) {
        return approx_unique(m) + approx_assoc(m);
    }
};