#pragma once
#include "rad_headers.h"
// set up a counter class to handle atomic operations
// Define an enum with unique indices
enum barcode_counts {
    forw = 0,         // Forward count
    forw_concat = 1,  // Forward concatenation count
    rev = 2,          // Reverse count
    rev_concat = 3,   // Reverse concatenation count
    total = 4,         // Total count
    corrected = 5,     // Corrected count
    filtered = 6
};

// building a thread-safe counter to be used during count updating processes
struct counter {
    // Use an array of 7 atomic ints.
    std::array<std::atomic<int>, 7> counts;

    // Default constructor: initialize all counters to 0.
    counter() : counts{{0, 0, 0, 0, 0, 0, 0}} { }

    // A constructor that initializes all counters to a given value.
    counter(int init) : counts{{init, init, init, init, init, init, init}} { }

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
    std::tuple<int, int, int, int, int, int, int> load_all() const {
        return std::make_tuple(
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
            std::size_t res = std::hash<uint16_t>()(seq.length);
            for (const auto &b : seq.bits) {
                boost::hash_combine(res, b);
            }
            return res;
        }
    };

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
        if (seq.bits.empty()) return mutations;
        int sequence_length = seq.length;
        int64_t original_value = seq.bits[0];  // Assumes one chunk.
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

    // generate mutated barcode candidates across multiple rounds, returning a unique set
    std::unordered_set<int64_seq> generate_mutated_barcodes(const int64_seq &seq, int mutation_rounds = 1) {
        if (seq.bits.empty()) return {};
        int sequence_length = seq.length;
        int64_seq initial = seq;
        std::unordered_set<int64_seq> all_mutations{ 
            initial 
        };
        std::unordered_set<int64_seq> current_round{ 
            initial
        };
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

    //generating shifted barcodes
    std::unordered_set<int64_seq> generate_shifted_barcodes(const int64_seq &seq, int shift) {
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
    }

};

template<typename key, typename value> struct bc_multimap:std::unordered_multimap<key, value> {
    using base = std::unordered_multimap<key, value>;
    using base::base;
    // helper to extract key from either a key or a value
    static inline const key& key_of(key const &k) noexcept { 
        return k;
    }
    static inline const key& key_of(value const &v) noexcept { 
        return v.barcode; 
    }

    // check whitelist for this barcode's existence
    template<typename T> bool check_wl_for(T const &x) const {
        if constexpr (std::is_same_v<T,int64_seq> || std::is_same_v<T,barcode_entry>) {
            return this->find(key_of(x)) != this->end();
        } else {
            // container of keys
            for (auto const &y : x) {
                if (check_wl_for(y)) return true;
            }
            return false;
        }
    }

    // return the matching barcodes from this wl
    template<typename T> std::unordered_set<int64_seq> return_matching_barcodes(T const &x) const {
        std::unordered_set<int64_seq> out;
        if constexpr (std::is_same_v<T,int64_seq> || std::is_same_v<T,barcode_entry>) {
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

    // insert a barcode entry into this wl
    template<typename T> void insert_bc_entry(const T &observed, const value &correct) & {
        this->emplace(key_of(observed), correct);
    }
    
    //insert a barcode entry into this wl, with barcode as key and value
    template<typename T> void insert_bc_entry(const T &observed) & {
        value v{};
        v.barcode = key_of(observed);
        this->emplace(key_of(observed), std::move(v));
    }
    
    //delete a barcode entry from this wl
    template<typename T> void remove_bc_entry(const T &observed) {
        this->erase(key_of(observed));
    }
    
    // return putative correct barcodes for an observed barcode we've seen before
    template<typename T> std::unordered_set<int64_seq> return_putative_correct_bcs(T const &x) const {
        std::unordered_set<int64_seq> out;
        if constexpr (std::is_same_v<T,int64_seq> ||std::is_same_v<T,barcode_entry>) {
            auto range = this->equal_range(key_of(x));
            for (auto it = range.first; it != range.second; ++it)
                out.insert(it->second.barcode);
        } else {
            for (auto const &y : x) {
                auto sub = return_putative_correct_bcs(y);
                out.insert(sub.begin(), sub.end());
            }
        }
        return out;
    }

    // add a count to a barcode_entry in this wl
    template<typename T> void update_bc_count(T const &x, barcode_counts slot = total) const {
        auto it = this->find(key_of(x));
        if (it != this->end())
            it->second.count.increment(slot);
    }

    // subtract a count from a barcode_entry
    template<typename T> void subtract_bc_count(T const &x, barcode_counts slot = total) const {
        auto it = this->find(key_of(x));
        if (it != this->end())
            it->second.count.subtract(slot);
    }

    // get the count for a barcode entry
    template<typename T> int get_bc_count(T const &x, barcode_counts slot = total) const {
        auto it = this->find(key_of(x));
        return it != this->end() ? it->second.count.load(slot) : 0;
    }

    std::vector<std::tuple<std::string,int,int,int,int,int,int,int>>
    summarize_counts(const std::unordered_set<key> *filter_keys = nullptr) const {
            using row = std::tuple<std::string,int,int,int,int,int,int,int>;
            // barcode, total_count, corrected_count, forw_count, rev_count, forw_concat, rev_concat, filtered
            std::vector<row> rows;

            // 1) Build the list of keys to process
            std::vector<key> keys;
            if (filter_keys) {
                keys.reserve(filter_keys->size());
                for (auto const &k : *filter_keys)
                    keys.push_back(k);
            } else {
                // collect distinct keys from the multimap
                std::unordered_set<key> seen;
                seen.reserve(this->size());
                for (auto const &kv : *this)
                    seen.insert(kv.first);
                keys.reserve(seen.size());
                for (auto const &k : seen)
                    keys.push_back(k);
            }

            rows.reserve(keys.size());

            // For each key, find its perfect‐match mapping and pull the counts
            for (auto const &bc : keys) {
                // get the sub‐range in the multimap
                auto range = this->equal_range(bc);
                int tot=0, corr=0, fwd=0, rev=0, fwd_c=0, rev_c=0, filt=0;

                if (range.first != range.second) {
                    // prefer the exact key if present
                    auto it = std::find_if(range.first, range.second, [&](auto const &p){ 
                        return p.first == bc; 
                    });
                    if (it == range.second){ 
                        it = range.first;
                    }
                    tot = it->second.count.load(barcode_counts::total);
                    corr = it->second.count.load(barcode_counts::corrected);
                    fwd = it->second.count.load(barcode_counts::forw);
                    rev = it->second.count.load(barcode_counts::rev);
                    fwd_c = it->second.count.load(barcode_counts::forw_concat);
                    rev_c = it->second.count.load(barcode_counts::rev_concat);
                    filt = it->second.count.load(barcode_counts::filtered);
                }
                rows.emplace_back(bc.bits_to_sequence(), tot, corr, fwd, rev, fwd_c, rev_c, filt);
            }
        return rows;
    }

    void write_wl_summary(std::ostream &out, const std::string class_id, const std::unordered_set<key> *filter_keys = nullptr) const{
        static bool header_printed = false;
        if(!header_printed) {
            out << "class_id,true_barcode,total_count,corrected_count,forw_count,rev_count,forw_concat_count,rev_concat_count,filtered_count\n";
            header_printed = true;
        }
            auto rows = summarize_counts(filter_keys);
            for (auto const & tpl : rows) {
                    auto const & bc = std::get<0>(tpl);
                    auto tot  = std::get<1>(tpl);
                    auto corr = std::get<2>(tpl);
                    auto forw = std::get<3>(tpl);
                    auto rev  = std::get<4>(tpl);
                    auto forw_concat = std::get<5>(tpl);
                    auto rev_concat = std::get<6>(tpl);
                    auto filtered = std::get<7>(tpl);
                    // print the row
                    out << class_id << ',' << bc << ','  << tot << ',' << corr << ',' << forw <<
                    ',' << rev << ',' << forw_concat << ',' << rev_concat << ',' << filtered << '\n';
        }
    }

    void write_wl_summary(const std::string &path = "", const std::unordered_set<key> *filter_keys = nullptr) const {
    std::ofstream ofs(path);
        if (!ofs) {
            throw std::runtime_error("Failed to open output file: " + path);
        }
        write_wl_summary(ofs, filter_keys);
    }
};

class whitelist {
    public:
        whitelist() = default;
        ~whitelist() = default;

    struct wl_entry {
        std::unordered_set<int64_seq> true_ref;
        bc_multimap<int64_seq, barcode_entry> true_bcs, global_bcs, filter_bcs;
    
        void generate_mismatch_barcodes(int shift, int mutation_rounds,bool verbose, int nthreads = 1) {
            if (verbose) {
                #pragma omp critical
                std::cout << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for true barcodes...\n";
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
            std::cout << "Detected element " << islist << " a bit-configured whitelist!" << std::endl;
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
                // decimal -> bits
                try {
                    int64_t code = std::stoll(tok);
                    seq.length = default_length;
                    seq.bits   = { code };
                } catch (...) { continue; }
            } else {
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