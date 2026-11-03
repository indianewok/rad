#pragma once
#include "rad_headers.h"
#include "packed_sequence.hpp"
#include "whitelist_relevance.hpp"

/**
 * @enum barcode_counts
 * @brief enumeration for different barcode count types
 * This enum defines indices for various barcode count categories used in the `counter` class.
 * The categories are raw, total, and corrected assignments used by barcode QC.
 */
enum barcode_counts {
    raw = 0,          // Raw count
    total = 1,        // Total count and QC rejection state
    corrected = 2     // Corrected count
};

/**
 * @struct counter
 * @brief Thread-safe counter class using atomic integers for multiple barcode count types
 * This class provides atomic operations for counting various barcode-related metrics.
 * It uses an array of 3 atomic integers to store counts for different categories defined in the
 * `barcode_counts` enum. The class supports incrementing, decrementing, loading individual counts,
 * and loading all counts as a tuple.
 */
struct counter {
    std::array<std::atomic<int>, 3> counts;

    // Default constructor: initialize all counters to 0.
    counter() : counts{{0, 0, 0}} { }

    // A constructor that initializes all counters to a given value
    counter(int init) : counts{{init, init, init}} { }

    // Copy constructor: copy each counter's current value
    counter(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
    }

    // Assignment operator: copy each counter's value
    counter& operator=(const counter &other) {
        for (size_t i = 0; i < counts.size(); ++i) {
            counts[i].store(other.counts[i].load(std::memory_order_relaxed),
                              std::memory_order_relaxed);
        }
        return *this;
    }

    // Increment the counter at the specified index
    void increment(int index) {
        counts[index].fetch_add(1, std::memory_order_relaxed);
    }

    // Overload: increment using one of the enum values for clarity
    void increment(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(1, std::memory_order_relaxed);
    }

    // Decrement the counter at the specified slot
    void subtract(int index) {
        counts[index].fetch_add(-1, std::memory_order_relaxed);
    }

    // Overload: decrement using an enum value for clarity
    void subtract(barcode_counts index) {
        counts[static_cast<int>(index)].fetch_add(-1, std::memory_order_relaxed);
    }

    // Load the value of the counter at the specified index
    int load(int index) const {
        return counts[index].load(std::memory_order_relaxed);
    }

    // Overload: load using an enum value for clarity
    int load(barcode_counts index) const {
        return counts[static_cast<int>(index)].load(std::memory_order_relaxed);
    }

    // Load all counter values into a tuple
    std::tuple<int, int, int> load_all() const {
        return std::make_tuple(
            counts[barcode_counts::raw].load(std::memory_order_relaxed),
            counts[barcode_counts::total].load(std::memory_order_relaxed),
            counts[barcode_counts::corrected].load(std::memory_order_relaxed)
        );
    }
};

/**
 * @struct barcode_entry
 * @brief Struct representing a barcode entry with its sequence, count, and filtering status
 * @param barcode `int64_seq` representing the barcode sequence
 * @param count `counter` atomic counter for tracking counts
 * @param filtered `bool` flag indicating if the barcode is filtered
 */
struct barcode_entry {
    int64_seq barcode;
    mutable counter count;
    bool filtered;

    barcode_entry()
      : barcode(0,1),
       count(0),
       filtered(false)
       { }

    bool operator==(const barcode_entry &o) const noexcept {
        return barcode == o.barcode;
    }

    bool is_valid() const noexcept {
        return barcode.is_valid();
    }

};

// setting up hashing for int64_seq & barcode_entry
namespace std {
    template<>
    struct hash<barcode_entry> {
        size_t operator()(barcode_entry const& bc) const noexcept {
            // hash the barcode field
            return std::hash<int64_seq>()(bc.barcode);
        }
    };
}

namespace mutation_tools {
    template<class Word>
    using unsigned_word_t = std::conditional_t<(sizeof(Word) > sizeof(uint64_t)), uint128_t, uint64_t>;

    template<class Word>
    constexpr int word_capacity = static_cast<int>(sizeof(unsigned_word_t<Word>) * 4);

    inline uint64_t pattern_mask(int length) noexcept {
        return length == 64 ? ~uint64_t{0} : (uint64_t{1} << length) - 1;
    }

    template<class Word>
    constexpr Word base_mask(int length) noexcept {
        return length == 0 ? Word{0}
             : length >= static_cast<int>(sizeof(Word) * 4) ? ~Word{0}
             : (Word{1} << (2 * length)) - 1;
    }

    template<class Word>
    struct raw_word_hash {
        size_t operator()(Word value) const noexcept {
            uint64_t low = static_cast<uint64_t>(value);
            if constexpr (sizeof(Word) > sizeof(uint64_t)) {
                const uint64_t high = static_cast<uint64_t>(value >> 64);
                low ^= high + 0x9e3779b97f4a7c15ULL + (low << 6) + (low >> 2);
            }
            return std::hash<uint64_t>{}(low);
        }
    };

    // Global edit distance. Packed words use two bits per base; Myers state
    // uses one bit per pattern position and therefore remains uint64_t.
    template<class PatternWord, class TextWord>
    int bit_ld(PatternWord pattern_code, TextWord text_code, int n, int max_dist) {
        if (n <= 0 || n > 64 || n > word_capacity<PatternWord> ||
            n > word_capacity<TextWord> || max_dist < 0) return -1;
        const auto pattern = static_cast<unsigned_word_t<PatternWord>>(pattern_code);
        const auto text = static_cast<unsigned_word_t<TextWord>>(text_code);
        uint64_t peq[4] = {0, 0, 0, 0};
        for (int i = 0; i < n; ++i) {
            const unsigned nuc = static_cast<unsigned>((pattern >> (2 * (n - 1 - i))) & 3);
            peq[nuc] |= uint64_t{1} << i;
        }
        const uint64_t mask = pattern_mask(n);
        const uint64_t high_bit = uint64_t{1} << (n - 1);
        uint64_t pv = mask, mv = 0;
        int score = n;
        for (int j = 0; j < n; ++j) {
            const unsigned nuc = static_cast<unsigned>((text >> (2 * (n - 1 - j))) & 3);
            const uint64_t eq = peq[nuc];
            const uint64_t xv = eq | mv;
            const uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
            uint64_t ph = mv | ~(xh | pv);
            const uint64_t mh = pv & xh;
            if (ph & high_bit) ++score;
            if (mh & high_bit) --score;
            if (score - (n - j - 1) > max_dist) return -1;
            ph = (ph << 1) | uint64_t{1};
            pv = ((mh << 1) | ~(xv | ph)) & mask;
            mv = ph & xv;
        }
        return score <= max_dist ? score : -1;
    }

    // Find a pattern in a text of either packed width. The zero boundary
    // term in ph permits a free text prefix, retaining infix matching.
    template<class PatternWord, class TextWord>
    int bit_partial_match(PatternWord pattern_code, TextWord text_code,
                          int pattern_len, int text_len, int max_dist,
                          int* best_end = nullptr, int* n_best = nullptr) {
        if (pattern_len <= 0 || text_len <= 0 || pattern_len > 64 || text_len > 64 ||
            pattern_len > word_capacity<PatternWord> || text_len > word_capacity<TextWord> ||
            max_dist < 0) return -1;
        const auto pattern = static_cast<unsigned_word_t<PatternWord>>(pattern_code);
        const auto text = static_cast<unsigned_word_t<TextWord>>(text_code);
        uint64_t peq[4] = {0, 0, 0, 0};
        for (int i = 0; i < pattern_len; ++i) {
            const unsigned nuc = static_cast<unsigned>((pattern >> (2 * (pattern_len - 1 - i))) & 3);
            peq[nuc] |= uint64_t{1} << i;
        }
        const uint64_t mask = pattern_mask(pattern_len);
        const uint64_t high_bit = uint64_t{1} << (pattern_len - 1);
        uint64_t pv = mask, mv = 0;
        int score = pattern_len;
        int min_score = pattern_len, min_end = -1, min_count = 0;
        const bool keep_position = best_end != nullptr || n_best != nullptr;
        for (int j = 0; j < text_len; ++j) {
            const unsigned nuc = static_cast<unsigned>((text >> (2 * (text_len - 1 - j))) & 3);
            const uint64_t eq = peq[nuc] & mask;
            const uint64_t xv = eq | mv;
            const uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
            uint64_t ph = (mv | ~(xh | pv)) & mask;
            const uint64_t mh = (pv & xh) & mask;
            if (ph & high_bit) ++score;
            if (mh & high_bit) --score;
            if (score < min_score) {
                min_score = score;
                if (keep_position) {
                    min_end = j;
                    min_count = 1;
                }
            } else if (keep_position && score == min_score && min_end >= 0) {
                ++min_count;
            }
            if (!keep_position && min_score == 0) return 0;
            ph <<= 1;
            pv = ((mh << 1) | ~(xv | ph)) & mask;
            mv = ph & xv;
        }
        if (best_end) *best_end = min_end;
        if (n_best) *n_best = min_count;
        return min_score <= max_dist ? min_score : -1;
    }

    inline int int64_lvdist(const int64_seq& query, const int64_seq& target, int max_dist,
                            int* best_end = nullptr, int* n_best = nullptr) {
        if (!query.is_valid() || !target.is_valid()) return -1;
        const int64_seq& pattern = query.length <= target.length ? query : target;
        const int64_seq& text = query.length <= target.length ? target : query;
        return pattern.visit([&](auto pattern_word) {
            return text.visit([&](auto text_word) {
                return bit_partial_match(pattern_word, text_word, pattern.length, text.length,
                                         max_dist, best_end, n_best);
            });
        });
    }

    inline std::map<int, std::unordered_set<int64_seq>> int64_lvdist(
        const int64_seq& query, const std::unordered_set<int64_seq>& targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (!query.is_valid()) return results;
        for (const auto& target : targets) {
            const int distance = int64_lvdist(query, target, max_dist);
            if (distance >= 0) results[distance].insert(target);
        }
        return results;
    }

    inline std::map<int, std::unordered_set<int64_seq>> int64_lvdist(
        const int64_seq& query, const std::vector<barcode_entry>& targets, int max_dist = 4) {
        std::map<int, std::unordered_set<int64_seq>> results;
        if (!query.is_valid()) return results;
        for (const auto& target : targets) {
            const int distance = int64_lvdist(query, target.barcode, max_dist);
            if (distance >= 0) results[distance].insert(target.barcode);
        }
        return results;
    }

    template<class Word, class Consumer>
    void for_each_point_mutation(Word original_value, int sequence_length, Consumer&& consume) {
        if (sequence_length <= 0 || sequence_length > word_capacity<Word>) return;
        for (int pos = 0; pos < sequence_length; ++pos) {
            const int shift = 2 * pos;
            const Word original_nuc = (original_value >> shift) & Word{3};
            const Word clear_mask = ~(Word{3} << shift);
            for (Word nuc = 0; nuc < 4; ++nuc) {
                if (nuc != original_nuc)
                    consume((original_value & clear_mask) | (nuc << shift));
            }
        }
    }

    template<class Word, std::enable_if_t<std::is_same_v<Word, uint64_t> ||
                                        std::is_same_v<Word, uint128_t>, int> = 0>
    std::unordered_set<Word, raw_word_hash<Word>> generate_point_mutations_raw(
        Word original_value, int sequence_length) {
        std::unordered_set<Word, raw_word_hash<Word>> mutations;
        if (sequence_length <= 0 || sequence_length > word_capacity<Word>) return mutations;
        mutations.reserve(static_cast<size_t>(sequence_length) * 3);
        for_each_point_mutation(original_value, sequence_length, [&](Word raw) { mutations.insert(raw); });
        return mutations;
    }

    // Preserve the legacy raw signed-word API while doing all shifts unsigned.
    inline std::unordered_set<int64_t> generate_point_mutations_raw(int64_t original_value, int sequence_length) {
        const auto raw = generate_point_mutations_raw(static_cast<uint64_t>(original_value), sequence_length);
        return std::unordered_set<int64_t>(raw.begin(), raw.end());
    }

    inline std::unordered_set<int64_seq> generate_point_mutations(const int64_seq& seq) {
        std::unordered_set<int64_seq> mutations;
        if (!seq.is_valid()) return mutations;
        mutations.reserve(static_cast<size_t>(seq.length) * 3);
        seq.visit([&](auto original) {
            for_each_point_mutation(original, seq.length, [&](auto raw) {
                mutations.insert(int64_seq::from_word(raw, seq.length));
            });
        });
        return mutations;
    }

    // Enumerate each set of Hamming substitutions once, independent of width.
    inline std::unordered_set<int64_seq> generate_mutated_barcodes(const int64_seq& seq, int mutation_rounds = 1) {
        std::unordered_set<int64_seq> result;
        if (!seq.is_valid() || mutation_rounds <= 0) return result;
        const int length = seq.length;
        const int max_distance = std::min(mutation_rounds, length);
        size_t expected = static_cast<size_t>(length) * 3;
        if (max_distance >= 2) expected += static_cast<size_t>(length) * (length - 1) * 9 / 2;
        if (max_distance >= 3) expected += static_cast<size_t>(length) * (length - 1) * (length - 2) * 27 / 6;
        result.reserve(expected);
        seq.visit([&](auto original) {
            using Word = std::decay_t<decltype(original)>;
            const auto enumerate = [&](const auto& self, int first, int remaining, Word raw) -> void {
                if (remaining == 0) {
                    result.insert(int64_seq::from_word(raw, seq.length));
                    return;
                }
                const int last = length - remaining;
                for (int position = first; position <= last; ++position) {
                    const int shift = 2 * position;
                    const Word mask = Word{3} << shift;
                    const Word old_nuc = (original & mask) >> shift;
                    for (Word nuc = 0; nuc < 4; ++nuc) {
                        if (nuc != old_nuc)
                            self(self, position + 1, remaining - 1, (raw & ~mask) | (nuc << shift));
                    }
                }
            };
            for (int distance = 1; distance <= max_distance; ++distance)
                enumerate(enumerate, 0, distance, original);
        });
        return result;
    }

    // Retain the legacy prefix/suffix-replacement candidate set.
    inline std::unordered_set<int64_seq> generate_shifted_barcodes(const int64_seq& seq, int shift) {
        static thread_local std::unordered_set<int64_seq> result;
        result.clear();
        if (!seq.is_valid() || shift <= 0 || shift >= seq.length) return result;
        if (shift >= static_cast<int>(sizeof(size_t) * 4))
            throw std::length_error("Shifted-barcode enumeration exceeds addressable size");
        const size_t combinations = size_t{1} << (2 * shift);
        if (combinations > result.max_size() / 2)
            throw std::length_error("Shifted-barcode enumeration exceeds container capacity");
        result.reserve(combinations * 2);
        const int core_length = seq.length - shift;
        seq.visit([&](auto original) {
            using Word = std::decay_t<decltype(original)>;
            const Word core = original >> (2 * shift);
            for (size_t code = 0; code < combinations; ++code) {
                result.insert(int64_seq::from_word((Word{code} << (2 * core_length)) | core, seq.length));
                result.insert(int64_seq::from_word((core << (2 * shift)) | Word{code}, seq.length));
            }
        });
        result.erase(seq);
        return result;
    }

    inline std::unordered_set<int64_seq> generate_lv_barcodes(const int64_seq& seq, int edit_rounds = 1) {
        if (!seq.is_valid()) return {};
        std::unordered_set<int64_seq> all_mutations{seq}, current_round{seq};
        for (int round = 0; round < edit_rounds; ++round) {
            std::unordered_set<int64_seq> next_round;
            const auto add = [&](const int64_seq& candidate) {
                if (all_mutations.insert(candidate).second) next_round.insert(candidate);
            };
            for (const auto& candidate : current_round) {
                const int length = candidate.length;
                candidate.visit([&](auto original) {
                    using Word = std::decay_t<decltype(original)>;
                    for_each_point_mutation(original, length, [&](Word raw) {
                        add(int64_seq::from_word(raw, length));
                    });
                    if (length > 1) {
                        for (int position = 0; position < length; ++position) {
                            const Word low = original & base_mask<Word>(position);
                            const Word high = position + 1 == word_capacity<Word>
                                ? Word{0} : original >> (2 * (position + 1));
                            add(int64_seq::from_word((low | (high << (2 * position))) &
                                                    base_mask<Word>(length - 1), length - 1));
                        }
                    }
                });
                if (length < 64) {
                    const auto insert = [&](auto original) {
                        using Word = std::decay_t<decltype(original)>;
                        original &= base_mask<Word>(length);
                        for (int position = 0; position <= length; ++position) {
                            const Word low = original & base_mask<Word>(position);
                            const Word high = original >> (2 * position);
                            const Word shifted_high = position + 1 == word_capacity<Word>
                                ? Word{0} : high << (2 * (position + 1));
                            for (Word nuc = 0; nuc < 4; ++nuc)
                                add(int64_seq::from_word(low | (nuc << (2 * position)) | shifted_high, length + 1));
                        }
                    };
                    if (length < 32) insert(candidate.as_uint64());
                    else insert(candidate.as_uint128());
                }
            }
            if (next_round.empty()) break;
            current_round = std::move(next_round);
        }
        all_mutations.erase(seq);
        return all_mutations;
    }

    inline bool detect_hp(const std::string& sequence, int max_hp) {
        if (sequence.empty() || max_hp <= 0) return false;
        const int length = static_cast<int>(sequence.size());
        if (length < max_hp) return false;
        bool leading = true, trailing = true;
        for (int i = 1; i < max_hp; ++i) {
            leading = leading && std::toupper(sequence[i]) == std::toupper(sequence.front());
            trailing = trailing && std::toupper(sequence[length - 1 - i]) == std::toupper(sequence.back());
        }
        return leading || trailing;
    }

    inline bool detect_hp(const int64_seq& sequence, int max_hp) {
        return sequence.is_valid() && detect_hp(sequence.bits_to_sequence(), max_hp);
    }
}

#include "barcode_storage.hpp"
using bc_map = bc_multimap<int64_seq, barcode_entry>;

/**
 * @brief Sparse bitmask for validating joint barcode pairs
 *
 * Stores a dense M x N bit matrix indexed by barcode whitelist position.
 * For Visium HD: 3350 x 3350 = ~1.4 MB.  Bit (i, j) is set iff the
 * (BC1_i, BC2_j) pair was observed in the valid-pairs whitelist.
 *
 * The sequence-to-index lookup maps are built once at load time from the
 * ordered whitelist files.  Validation is a single array access.
 */
struct spat_mask {
    std::vector<bool> bits;
    size_t rows = 0;   // BC1 axis size
    size_t cols = 0;   // BC2 axis size

    // sequence -> axis index
    template <class Index> class axis_index {
        using narrow_table = phmap::flat_hash_map<packed_barcode<uint64_t>, Index>;
        using wide_table = phmap::flat_hash_map<packed_barcode<uint128_t>, Index>;
        barcode_storage_detail::table_pair<narrow_table, wide_table> entries_;
    public:
        void reserve(size_t size) { entries_.narrow().reserve(size); }
        void set(const int64_seq& sequence, Index index) {
            if (!sequence.is_valid()) return;
            if (sequence.is_wide()) {
                entries_.ensure_wide()[barcode_storage_detail::pack_key<uint128_t>(sequence)] = index;
            } else {
                entries_.narrow()[barcode_storage_detail::pack_key<uint64_t>(sequence)] = index;
            }
        }
        std::optional<Index> get(const int64_seq& sequence) const {
            if (!sequence.is_valid()) return std::nullopt;
            if (sequence.is_wide()) {
                const auto* table = entries_.wide();
                if (!table) return std::nullopt;
                const auto it = table->find(barcode_storage_detail::pack_key<uint128_t>(sequence));
                if (it != table->end()) return it->second;
            } else {
                const auto& table = entries_.narrow();
                const auto it = table.find(barcode_storage_detail::pack_key<uint64_t>(sequence));
                if (it != table.end()) return it->second;
            }
            return std::nullopt;
        }
        bool has_wide_storage() const { return entries_.wide() != nullptr; }
    };
    axis_index<uint16_t> bc1_to_idx, bc2_to_idx;

    bool empty() const { return bits.empty(); }

    /// Check whether a (BC1, BC2) pair is valid.
    bool check(const int64_seq& bc1, const int64_seq& bc2) const {
        const auto row = bc1_to_idx.get(bc1);
        const auto col = bc2_to_idx.get(bc2);
        if (!row || !col) return false;
        return bits[static_cast<size_t>(*row) * cols + *col];
    }

    /// Number of valid pairs (set bits).
    size_t count_valid() const {
        size_t n = 0;
        for (bool b : bits) n += b;
        return n;
    }

    /// Load a binary matrix CSV and the two axis whitelist files.
    /// The matrix has no header; row = BC1 index, col = BC2 index, values are 0 or 1.
    /// The whitelist files define the index ordering (one barcode per line, first line may be "whitelist_bcs").
    bool load(const std::string& matrix_path,
              const std::string& bc1_wl_path,
              const std::string& bc2_wl_path,
              bool verbose) {

        auto read_ordered_wl = [](const std::string& path) -> std::vector<std::string> {
            std::vector<std::string> out;
            gzFile fp = gzopen(path.c_str(), "r");
            if (!fp) return out;
            char buf[1024];
            bool first = true;
            while (gzgets(fp, buf, sizeof(buf))) {
                std::string line(buf);
                // trim trailing whitespace/newline
                while (!line.empty() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' '))
                    line.pop_back();
                if (line.empty()) continue;
                // strip anything after a comma (some whitelists have extra columns)
                auto comma = line.find(',');
                if (comma != std::string::npos) line = line.substr(0, comma);
                if (first) { first = false; if (line == "whitelist_bcs") continue; }
                out.push_back(line);
            }
            gzclose(fp);
            return out;
        };

        auto bc1_list = read_ordered_wl(bc1_wl_path);
        auto bc2_list = read_ordered_wl(bc2_wl_path);
        if (bc1_list.empty() || bc2_list.empty()) {
            std::cerr << "[spat_mask] Failed to load axis whitelists\n";
            return false;
        }

        rows = bc1_list.size();
        cols = bc2_list.size();

        bc1_to_idx.reserve(rows);
        bc2_to_idx.reserve(cols);
        for (size_t i = 0; i < rows; ++i) {
            int64_seq bits_i(bc1_list[i]);
            bc1_to_idx.set(bits_i, static_cast<uint16_t>(i));
        }
        for (size_t j = 0; j < cols; ++j) {
            int64_seq bits_j(bc2_list[j]);
            bc2_to_idx.set(bits_j, static_cast<uint16_t>(j));
        }

        // Read the matrix CSV
        bits.assign(rows * cols, false);
        std::ifstream in(matrix_path);
        if (!in.is_open()) {
            std::cerr << "[spat_mask] Failed to open matrix: " << matrix_path << "\n";
            return false;
        }

        size_t valid_count = 0;
        std::string line;
        for (size_t r = 0; r < rows && std::getline(in, line); ++r) {
            size_t c = 0;
            size_t pos = 0;
            while (pos < line.size() && c < cols) {
                char ch = line[pos];
                if (ch == '1') {
                    bits[r * cols + c] = true;
                    ++valid_count;
                    ++c;
                } else if (ch == '0') {
                    ++c;
                }
                // skip commas and whitespace
                ++pos;
            }
        }

        if (verbose) {
            std::cout << "[spat_mask] Loaded " << rows << "x" << cols
                      << " matrix (" << valid_count << " valid pairs) from "
                      << matrix_path << "\n";
        }
        return true;
    }
};

/**
 * @brief Whitelist class managing true and filtered barcodes
 * @param wl_entry: Struct containing true barcodes, filtered barcodes, and global barcodes
 */
class whitelist {
    public:
        whitelist() = default;
        ~whitelist() = default;

/**
 * @brief Struct representing a whitelist entry with true barcodes, filtered barcodes, and global barcodes
 * @param true_bcs: `bc_multimap<int64_seq, barcode_entry>` for true barcodes
 * @param filter_bcs: `bc_multimap<int64_seq, barcode_entry>` for filtered barcodes
 * @param global_bcs: `bc_flat_set` for global barcodes
 * @param spat_wl: optional `spat_mask` for validating joint barcode pairs
 */
    struct wl_entry {
        bc_multimap<int64_seq, barcode_entry> true_bcs, filter_bcs;
        bc_flat_set global_bcs;
        std::optional<spat_mask> spat_wl;

/**
 * @brief Seed-based barcode candidate index
 *
 * Splits each barcode into three partitions and indexes their packed values
 * by barcode length. Same-length queries produce a shortlist; full distance
 * scoring and the exhaustive fallback establish accepted correction calls.
 * Keys use uint64_t for seeds up to 22 bases. Values are handles into the
 * owning width-specific whitelist tables.
 */
        struct seed_idx {
            struct key {
                uint64_t bits = 0;
                uint16_t length = 0;
                uint8_t partition = 0;
                bool operator==(const key& other) const noexcept {
                    return bits == other.bits && length == other.length && partition == other.partition;
                }
            };
            struct key_hash {
                size_t operator()(const key& k) const noexcept {
                    return std::hash<uint64_t>{}(k.bits ^
                        (uint64_t(k.length) << 48) ^ (uint64_t(k.partition) << 56));
                }
            };
            std::unordered_map<key, std::vector<barcode_entry_handle>, key_hash> index;
            bool ready = false;

            static std::array<std::pair<size_t, size_t>, 3> partitions(size_t len) {
                std::array<std::pair<size_t, size_t>, 3> windows{};
                size_t start = 0;
                for (size_t i = 0; i < 3; ++i) {
                    const size_t end = start + len / 3 + (i < len % 3 ? 1 : 0);
                    windows[i] = {start, end};
                    start = end;
                }
                return windows;
            }

            static uint64_t extract_bits(const int64_seq& seq, size_t start, size_t len) {
                return seq.slice(start, len).as_uint64();
            }

            void build(const std::vector<barcode_entry_handle>& entries, bool verbose) {
                index.clear();
                ready = false;
                index.reserve(entries.size() * 3);
                for (const auto entry : entries) {
                    if (!entry) continue;
                    const auto barcode = entry.barcode();
                    if (!barcode.is_valid() || barcode.length < 3) continue;
                    const auto windows = partitions(barcode.length);
                    for (size_t pid = 0; pid < windows.size(); ++pid) {
                        const auto [start, end] = windows[pid];
                        key k{extract_bits(barcode, start, end - start), barcode.length,
                              static_cast<uint8_t>(pid)};
                        index[k].push_back(entry);
                    }
                }
                ready = !index.empty();
                if (verbose) {
                    std::cout << "[seed_idx] Built index: " << index.size() << " keys from "
                              << entries.size() << " barcodes (partitions=3)\n";
                }
            }

            std::vector<barcode_entry_handle> query(const int64_seq& barcode) const {
                std::vector<barcode_entry_handle> candidates;
                if (!ready || !barcode.is_valid() || barcode.length < 3) return candidates;
                const auto windows = partitions(barcode.length);
                for (size_t pid = 0; pid < windows.size(); ++pid) {
                    const auto [start, end] = windows[pid];
                    key k{extract_bits(barcode, start, end - start), barcode.length,
                          static_cast<uint8_t>(pid)};
                    const auto it = index.find(k);
                    if (it != index.end()) candidates.insert(candidates.end(), it->second.begin(), it->second.end());
                }
                return candidates;
            }
        };

        seed_idx true_seeds, global_seeds;

        void build_seed_idx(bool verbose) {
            true_seeds.build(true_bcs.get_unique_entries(), verbose);
            constexpr size_t max_global_seed_barcodes = 10000;
            if (global_bcs.unique_val_size() > max_global_seed_barcodes) {
                global_seeds = seed_idx{};
                if (verbose) {
                    std::cout << "[seed_idx] Skipping global index: " << global_bcs.unique_val_size()
                              << " unique barcodes exceed the 10,000-barcode cap\n";
                }
            } else {
                global_seeds.build(global_bcs.get_unique_entries(), verbose);
            }
        }

        bool has_seed_idx() const noexcept {
            return true_seeds.ready || global_seeds.ready;
        }

        std::vector<barcode_entry_handle> query_bc_seeds(
            const int64_seq& query, const std::string& src) const {
            if (src == "true")  return true_seeds.query(query);
            if (src == "global") return global_seeds.query(query);
            auto r = true_seeds.query(query);
            auto g = global_seeds.query(query);
            r.insert(r.end(), g.begin(), g.end());
            return r;
        }

/**
 * @brief Part table of the true whitelist: a lossless filter in front of the exhaustive scan
 *
 * Each barcode is cut into PARTS parts (seed_idx::partitions). An alignment with at most PARTS - 1 edits leaves at
 * least one part unchanged, so a barcode within max_dist <= PARTS - 1 edits of an infix of a longer window has one
 * part that occurs exactly in the window. The table maps the key of each part (its first key_len bases) to the
 * barcodes that hold it, one group for each barcode length. scan() gives the same match set as int64_lvdist over
 * every unique entry. Built once after the whitelist is loaded; read-only afterwards.
 */
        struct part_idx {
            static constexpr int PARTS = 3;    // parts for each barcode: lossless for max_dist <= PARTS - 1
            static constexpr int KEY_MAX = 8;  // bases of a part in its key: bounds a table at 4^8 slots for each part

            struct group {
                int m = 0;                                     // barcode length of the group
                std::vector<barcode_entry_handle> entries;     // the barcodes of this length
                int part_start[PARTS] = {0, 0, 0};
                int key_len[PARTS] = {0, 0, 0};
                std::vector<uint32_t> offset[PARTS];           // rows offset[k][key] .. offset[k][key + 1) of list[k]
                std::vector<uint32_t> list[PARTS];             // indices into entries
            };
            std::vector<group> groups;
            size_t n_entries = 0;                              // unique entries at build time
            bool ready = false;

            /// int64_lvdist can give a distance for this entry (otherwise it returns -1 for every window)
            static bool comparable(barcode_entry_handle e) noexcept {
                return e && e.barcode().is_valid();
            }

            /// bases [start, start + key_len) of a len-base sequence whose first base is in the highest 2 bits
            static uint32_t key_of(const int64_seq& sequence, int start, int key_len) {
                return static_cast<uint32_t>(sequence.slice(start, key_len).as_uint64());
            }

            void build(const std::vector<barcode_entry_handle>& entries, bool verbose) {
                groups.clear();
                n_entries = entries.size();
                ready = false;
                std::map<int, size_t> by_len;
                for (const auto e : entries) {
                    if (!comparable(e)) continue;
                    const int m = e.barcode().length;
                    auto it = by_len.find(m);
                    if (it == by_len.end()) {
                        it = by_len.emplace(m, groups.size()).first;
                        groups.emplace_back();
                        groups.back().m = m;
                    }
                    groups[it->second].entries.push_back(e);
                }
                for (auto& g : groups) {
                    if (g.m < PARTS) {   // fewer bases than parts: use the full loop
                        groups.clear();
                        if (verbose) std::cout << "[part_idx] Skipping build: barcode length < " << PARTS << "\n";
                        return;
                    }
                    const auto windows = seed_idx::partitions(static_cast<size_t>(g.m));
                    for (int k = 0; k < PARTS; ++k) {
                        g.part_start[k] = static_cast<int>(windows[k].first);
                        g.key_len[k] = std::min<int>(static_cast<int>(windows[k].second - windows[k].first), KEY_MAX);
                        const size_t slots = size_t(1) << (2 * g.key_len[k]);
                        std::vector<uint32_t> count(slots + 1, 0);
                        for (const auto e : g.entries) {
                            ++count[key_of(e.barcode(), g.part_start[k], g.key_len[k]) + 1];
                        }
                        for (size_t s = 0; s < slots; ++s) count[s + 1] += count[s];
                        std::vector<uint32_t> fill(count.begin(), count.end() - 1);
                        g.list[k].assign(g.entries.size(), 0);
                        for (uint32_t i = 0; i < g.entries.size(); ++i) {
                            const uint32_t key = key_of(g.entries[i].barcode(), g.part_start[k], g.key_len[k]);
                            g.list[k][fill[key]++] = i;
                        }
                        g.offset[k] = std::move(count);
                    }
                }
                ready = !groups.empty();
                if (verbose) {
                    std::cout << "[part_idx] Built table: " << n_entries << " barcodes, " << groups.size()
                              << " length group(s)";
                    for (const auto& g : groups) {
                        std::cout << " [len=" << g.m << ": " << g.entries.size() << " barcodes, keys "
                                  << g.key_len[0] << "/" << g.key_len[1] << "/" << g.key_len[2] << "]";
                    }
                    std::cout << ", " << (bytes() / 1024) << " KiB\n";
                }
            }

            /// scan() gives the match set of the loop over every unique entry for this cap and this list size
            bool covers(int max_dist, size_t unique_size) const noexcept {
                return ready && max_dist <= PARTS - 1 && unique_size == n_entries;
            }

            size_t bytes() const noexcept {
                size_t b = sizeof(*this) + groups.capacity() * sizeof(group);
                for (const auto& g : groups) {
                    b += g.entries.capacity() * sizeof(barcode_entry_handle);
                    for (int k = 0; k < PARTS; ++k) {
                        b += (g.offset[k].capacity() + g.list[k].capacity()) * sizeof(uint32_t);
                    }
                }
                return b;
            }

            /// The matches of `win` against every entry: the pairs (barcode, distance) that the loop
            /// `for (entry : get_unique_entries()) int64_lvdist(win, entry->barcode, max_dist) >= 0` gives, in any order.
            void scan(const int64_seq& win, int max_dist, std::vector<std::pair<int64_seq, int>>& match) const {
                const int n = win.length;
                const bool packed = win.is_valid();
                std::vector<barcode_entry_handle> hit;   // entries already reported (a candidate can come up twice)
                for (const group& g : groups) {
                    if (!packed || n <= g.m) {   // the window is not longer than the barcode: compare every entry
                        for (const auto e : g.entries) {
                            const int res = mutation_tools::int64_lvdist(win, e.barcode(), max_dist);
                            if (res >= 0) match.emplace_back(e.barcode(), res);
                        }
                        continue;
                    }
                    hit.clear();
                    for (int k = 0; k < PARTS; ++k) {
                        const int kl = g.key_len[k];
                        const uint32_t* off = g.offset[k].data();
                        const uint32_t* rows = g.list[k].data();
                        for (int j = 0; j + kl <= n; ++j) {
                            const uint32_t key = key_of(win, j, kl);
                            for (uint32_t r = off[key], end = off[key + 1]; r < end; ++r) {
                                const auto e = g.entries[rows[r]];
                                if (std::find(hit.begin(), hit.end(), e) != hit.end()) continue;
                                const int res = mutation_tools::int64_lvdist(win, e.barcode(), max_dist);
                                if (res >= 0) {
                                    hit.push_back(e);
                                    match.emplace_back(e.barcode(), res);
                                }
                            }
                        }
                    }
                }
            }
        };

        part_idx true_parts;

        /// Build the part table of the true whitelist (call once, after the last change to true_bcs)
        void build_part_idx(bool verbose) {
            true_parts.build(true_bcs.get_unique_entries(), verbose);
        }

        template<class F> decltype(auto) with_wl(std::string_view src, F&& f) {
            if (src == "global") return f(global_bcs);
            return f(true_bcs);
        }
        template<class F> decltype(auto) with_wl(std::string_view src, F&& f) const {
            if (src == "global") return f(global_bcs);
            return f(true_bcs);
        }
/**
 * @brief Generate mismatch barcodes by applying mutations and shifts to existing true barcodes
 * @param mutation_rounds Number of mutation rounds to apply
*/
        void generate_mismatch_barcodes(int mutation_rounds, bool verbose, int nthreads = 1) {
                using namespace std::chrono;
                auto start_total = high_resolution_clock::now();

                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "[generate_mismatch_barcodes] Generating shifted and mutated barcodes for "
                            << "mutation_rounds = " << mutation_rounds << "\n";
                        std::cout << oss.str();
                    }
                }

                // TIMING: Collect original barcode entries
                auto start_collect = high_resolution_clock::now();

                auto originals = true_bcs.get_unique_entries();

                auto end_collect = high_resolution_clock::now();
                double collect_ms = duration<double, std::milli>(end_collect - start_collect).count();

                if (verbose) {
                    std::cout << "[generate_mismatch_barcodes] Collected " << originals.size()
                            << " unique originals in " << collect_ms << " ms\n";
                }

                // Memory tracking
                size_t initial_true_bcs_size = true_bcs.size();
                size_t initial_memory_estimate = initial_true_bcs_size * 48; // Rough estimate: 48 bytes per association

                if (verbose) {
                    std::cout << "[generate_mismatch_barcodes] Initial true_bcs size: " << initial_true_bcs_size
                            << " associations (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
                }

                // TIMING: Generate mutations and shifts
                auto start_generation = high_resolution_clock::now();

                size_t total_mutations_generated = 0;
                size_t mutations_added = 0;
                size_t mutations_rejected_global = 0;

                // Track timing for individual operations
                double total_mutation_generation_ms = 0;
                double total_mutation_lookup_ms = 0;
                double total_mutation_insert_ms = 0;

                auto start_loop = high_resolution_clock::now();

                for (size_t i = 0; i < originals.size(); ++i) {
                    const auto handle = originals[i];
                    barcode_entry original;
                    original.barcode = handle.barcode();
                    original.count = *handle.counts();
                    original.filtered = handle.filtered();
                    const auto& orig_bits = original.barcode;

                    // Time individual barcode processing for first few
                    auto barcode_start = high_resolution_clock::now();

                    // === LEVENSHTEIN MUTATION GENERATION ===
                    auto mutation_gen_start = high_resolution_clock::now();
                    auto mutated = mutation_tools::generate_lv_barcodes(orig_bits, mutation_rounds);
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
                            true_bcs.insert_bc_entry(m_bits, original);
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
                                << mutated.size() << " mutations (" << mutation_gen_time << "ms gen), "
                                << "total " << total_time << "ms\n";
                    }

                    // Memory usage check every 10000 barcodes
                    if (verbose && i % 1000 == 0 && i > 0) {
                        size_t current_true_bcs_size = true_bcs.size();
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
                size_t final_true_bcs_size = true_bcs.size();
                size_t total_growth = final_true_bcs_size - initial_true_bcs_size;
                size_t final_memory_estimate = final_true_bcs_size * 48;

                // TIMING: Sanity check
                auto start_sanity = high_resolution_clock::now();
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
                    std::cout << "Mutation generation: " << total_mutation_generation_ms << " ms ("
                            << (total_mutation_generation_ms/total_ms*100) << "%)\n";
                    std::cout << "Mutation lookups: " << total_mutation_lookup_ms << " ms ("
                            << (total_mutation_lookup_ms/total_ms*100) << "%)\n";
                    std::cout << "Mutation insertions: " << total_mutation_insert_ms << " ms ("
                            << (total_mutation_insert_ms/total_ms*100) << "%)\n";

                    std::cout << "\n=== GENERATION STATISTICS ===\n";
                    std::cout << "Processed " << originals.size() << " original barcodes\n";
                    std::cout << "Mutations: generated " << total_mutations_generated << ", added " << mutations_added
                            << ", rejected " << mutations_rejected_global << "\n";
                    std::cout << "Average mutations per barcode: " << (double)total_mutations_generated / originals.size() << "\n";
                    std::cout << "Time per barcode: " << (loop_ms / originals.size()) << " ms\n";

                    std::cout << "\n=== MEMORY USAGE ===\n";
                    std::cout << "Initial associations: " << initial_true_bcs_size
                            << " (~" << (initial_memory_estimate / 1024 / 1024) << " MB)\n";
                    std::cout << "Final associations: " << final_true_bcs_size
                            << " (~" << (final_memory_estimate / 1024 / 1024) << " MB)\n";
                    std::cout << "Growth: +" << total_growth << " associations (+~"
                            << ((final_memory_estimate - initial_memory_estimate) / 1024 / 1024) << " MB)\n";
                    std::cout << "Memory efficiency: " << (double)total_growth / (total_mutations_generated) * 100
                            << "% of generated items were unique and added\n";
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
/**
 * @brief Stream parsed barcodes from a specified kit name or file path
 * @param spec: Kit name or file path to load barcodes from
 * @param default_length: Default length of barcodes if not specified
 * @param consume: Callback receiving each parsed barcode
 */
    template <typename Consumer>
    static void stream_barcodes(std::string const &spec, uint16_t default_length, bool verbose,
                                Consumer&& consume) {
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

        // Keep source parsing identical for staged and direct imports.
        std::unique_ptr<std::istream> in;
        if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
            in = std::make_unique<igzstream>(path.c_str());
        } else {
            in = std::make_unique<std::ifstream>(path);
        }
        if (!in || !*in) {
            throw std::runtime_error("Failed to open whitelist: " + path);
        }

        std::string ln;
        while (std::getline(*in, ln)) {
            if (ln.empty()) continue;
            auto p = ln.find_first_of(",\t");
            std::string tok = (p == std::string::npos ? ln : ln.substr(0,p));

            int64_seq seq;
            if (isBitlist) {
                if (!parse_packed_decimal(tok, default_length, seq)) continue;
            } else {
                seq.sequence_to_bits(tok);
                if (seq.empty()) continue;
            }
            consume(seq);
        }
    }

    static std::unordered_set<int64_seq> load_barcodes(std::string const &spec, uint16_t default_length, bool verbose) {
        std::unordered_set<int64_seq> out;
        stream_barcodes(spec, default_length, verbose, [&](const int64_seq& seq) {
            out.insert(seq);
        });
        return out;
    }

    wl_entry import_single_whitelist(std::string const &spec, uint16_t default_length, bool verbose) {
        constexpr size_t global_threshold = 100000;
        wl_entry out;
        std::unordered_set<int64_seq> prefix;
        prefix.reserve(global_threshold);
        bool global = false;
        auto insert_global = [&](const int64_seq& seq) {
            if (seq.empty()) return;
            barcode_entry entry;
            entry.barcode = seq;
            entry.filtered = false;
            out.global_bcs.emplace(seq, std::move(entry));
        };
        stream_barcodes(spec, default_length, verbose, [&](const int64_seq& seq) {
            if (global) {
                insert_global(seq);
                return;
            }
            // Count the same distinct parsed keys as load_barcodes, including
            // its length-zero sentinel, before choosing the whitelist role.
            prefix.insert(seq);
            if (prefix.size() >= global_threshold) {
                global = true;
                for (const auto& buffered : prefix) insert_global(buffered);
                std::unordered_set<int64_seq>().swap(prefix);
            }
        });
        if (!global) {
            for (const auto& seq : prefix) {
                if (seq.empty()) continue;
                barcode_entry entry;
                entry.barcode = seq;
                entry.filtered = false;
                out.true_bcs.emplace(seq, std::move(entry));
            }
        }
        std::cout << "[populate_entry] Single whitelist populated "
                  << (global ? "global_bcs: " + std::to_string(out.global_bcs.size())
                             : "true_bcs: " + std::to_string(out.true_bcs.size()))
                  << " (classification buffer limited to " << global_threshold << " distinct keys).\n";
        return out;
    }
/**
 * @brief Populate a whitelist entry with true and global barcodes based on provided sets
 * @param out: Reference to the `wl_entry` to populate
 * @param sets: Vector of unordered sets of `int64_seq` representing barcode sets
 */
    void populate_entry(wl_entry &out, std::vector<std::unordered_set<int64_seq>> const &sets){
        if (sets.size() == 1) {
            auto const &A = sets[0];
            bool populate_global = false;
            size_t data_size = A.size();
            std::cout << "[populate_entry] Single whitelist with " << data_size << " barcodes detected.\n";
            if (A.size() >= 100000) {
                populate_global = true;
            }
            std::cout << "[populate_entry] Populating " << (populate_global ? "global_bcs" : "true_bcs") << "...\n";
            for (auto const &seq : A) {
                if (seq.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                if(populate_global){
                    out.global_bcs.emplace(seq, std::move(be));
                } else {
                    out.true_bcs.emplace(seq, std::move(be));
                }
            }
            std::cout << "[populate_entry] Population complete. Size of "
                      << (populate_global ? "global_bcs: " + std::to_string(out.global_bcs.size())
                                          : "true_bcs: " + std::to_string(out.true_bcs.size()))
                      << "\n";
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
                if (seq.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                out.true_bcs.emplace(seq, std::move(be));
            }
            for (auto const &seq : large) {
                //this line trims duplicate entries in global bcs
                if (!out.true_bcs.check_wl_for(seq)) {
                    if (seq.empty()) continue;
                    barcode_entry be{};
                    be.barcode   = seq;
                    be.filtered  = false;
                    out.global_bcs.emplace(seq, be);
                }
            }
        }
        else {
            // fallback: merge both into true_bcs
            for (auto const &seq : A) {
                if (seq.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                //be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
            }
            for (auto const &seq : B) {
                if (seq.empty()) continue;
                barcode_entry be;
                be.barcode = seq;
                be.filtered = false;
                //be.flags = "";
                out.true_bcs.emplace(seq, std::move(be));
            }
        }
    }

    void populate_entry(wl_entry &out, std::vector<std::unordered_set<int64_seq>> &&sets) {
        auto drain = [&](std::unordered_set<int64_seq>& source, auto& destination,
                         bool exclude_true = false) {
            while (!source.empty()) {
                auto current = source.begin();
                const auto& seq = *current;
                if (!seq.empty() && (!exclude_true || !out.true_bcs.check_wl_for(seq))) {
                    barcode_entry entry;
                    entry.barcode = seq;
                    entry.filtered = false;
                    destination.emplace(seq, std::move(entry));
                }
                source.erase(current);
            }
            std::unordered_set<int64_seq>().swap(source);
        };
        if (sets.size() == 1) {
            const size_t data_size = sets[0].size();
            const bool global = data_size >= 100000;
            std::cout << "[populate_entry] Single whitelist with " << data_size << " barcodes detected.\n";
            std::cout << "[populate_entry] Populating " << (global ? "global_bcs" : "true_bcs") << "...\n";
            if (global) drain(sets[0], out.global_bcs);
            else drain(sets[0], out.true_bcs);
            std::cout << "[populate_entry] Population complete. Size of "
                      << (global ? "global_bcs: " + std::to_string(out.global_bcs.size())
                                 : "true_bcs: " + std::to_string(out.true_bcs.size()))
                      << "\n";
            return;
        }
        // Choose both roles before draining changes either source cardinality.
        const size_t size_a = sets[0].size();
        const size_t size_b = sets[1].size();
        if (size_a != size_b) {
            const size_t small = size_a < size_b ? 0 : 1;
            const size_t large = 1 - small;
            drain(sets[small], out.true_bcs);
            drain(sets[large], out.global_bcs, true);
        } else {
            drain(sets[0], out.true_bcs);
            drain(sets[1], out.true_bcs);
        }
    }

    //import a whitelist from file--can be true barcodes or not.
    //cheap but easy way to import b/w both global or custom--just set
    //an arbitrary threshold for the number of true barcodes to be kept because ideally after a certain size of barcodes
    //you really just want to treat them both the same way
/**
 * @brief Import a whitelist from a specified field (kit name or file path)
 * @param field: Kit name or file path to import the whitelist from
 * @param default_length: Default length of barcodes if not specified
 * @return `wl_entry` containing the imported whitelist
 * @brief cheap but easy way to import b/w both global or custom--just set
 * an arbitrary threshold for the number of true barcodes to be kept because ideally after a certain size of barcodes
 * you really just want to treat them both the same way
 */
    // Stream a single catalog while retaining the exact global-lookup superset
    // for the discovery pass. At most 100,000 unique source keys are held to
    // establish the ORIGINAL catalog role. A retained subset must never turn
    // a large global catalog into a small true whitelist: those paths differ.
    wl_entry import_relevant_whitelist(const std::string& spec, bool verbose,
                                      uint16_t default_length,
                                      const whitelist_relevance::observations& observed,
                                      int mutation_distance) {
        constexpr size_t global_threshold = 100000;
        const std::string path = whitelist_utils::kit_to_path(spec);
        const bool bitlist = whitelist_utils::check_if_bitlist(spec, verbose, 10);
        std::unique_ptr<std::istream> in;
        if (path.size() > 3 && path.substr(path.size() - 3) == ".gz")
            in = std::make_unique<igzstream>(path.c_str());
        else
            in = std::make_unique<std::ifstream>(path);
        if (!in || !*in) throw std::runtime_error("Failed to open whitelist: " + path);
        whitelist_relevance::candidate_index index(observed, mutation_distance);
        wl_entry out;
        std::unordered_set<int64_seq> prefix;
        prefix.reserve(global_threshold);
        bool global = false;
        size_t rows = 0, parsed_rows = 0, wide_rows = 0;
        auto keep = [&](const int64_seq& sequence) {
            return !sequence.length || sequence.visit([&](auto word) {
                return index.contains(word, sequence.length);
            });
        };
        auto insert_global = [&](const int64_seq& sequence) {
            if (!keep(sequence)) return;
            barcode_entry entry;
            entry.barcode = sequence;
            out.global_bcs.emplace(sequence, std::move(entry));
        };
        std::string line;
        while (std::getline(*in, line)) {
            ++rows;
            if (line.empty()) continue;
            const auto delimiter = line.find_first_of(",\t");
            const std::string_view token(line.data(), delimiter == std::string::npos ? line.size() : delimiter);
            int64_seq sequence;
            if (bitlist) {
                if (!parse_packed_decimal(token, default_length, sequence)) continue;
            } else {
                sequence.sequence_to_bits(token);
            }
            ++parsed_rows;
            wide_rows += sequence.is_wide();
            if (global) {
                insert_global(sequence);
            } else {
                prefix.insert(sequence);
                if (prefix.size() >= global_threshold) {
                    global = true;
                    for (const auto& buffered : prefix) insert_global(buffered);
                    std::unordered_set<int64_seq>().swap(prefix);
                }
            }
        }
        if (in->bad()) throw std::runtime_error("Error while reading whitelist: " + path);
        if (!global) {
            for (const auto& sequence : prefix) {
                barcode_entry entry;
                entry.barcode = sequence;
                out.true_bcs.emplace(sequence, std::move(entry));
            }
        }
        std::cout << "[whitelist-load observed] source=" << path
                  << " rows=" << rows << " parsed_rows=" << parsed_rows
                  << " role=" << (global ? "global" : "true-full")
                  << " retained=" << (global ? out.global_bcs.size() : out.true_bcs.size())
                  << " observed_raw=" << observed.raw_size()
                  << " exact_window_keys=" << observed.exact_size()
                  << " mutation_cap=" << std::max(0, mutation_distance)
                  << " index_bytes=" << index.bytes()
                  << " wide_rows=" << wide_rows << '\n';
        return out;
    }

    wl_entry import_whitelist(std::string const &field, bool verbose, uint16_t default_length = 16,
                              const whitelist_relevance::observations* observed = nullptr,
                              int mutation_distance = 2) {
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

        if (observed && specs.size() == 1 && !observed->unsupported_long_queries)
            return import_relevant_whitelist(specs.front(), verbose, default_length,
                                             *observed, mutation_distance);
        if (observed) {
            std::cout << "[whitelist-load observed] using full loader: "
                      << (specs.size() != 1 ? "dual catalog roles require full source cardinalities"
                                          : "barcode queries longer than 64 bases") << '\n';
        }

        if (specs.size() == 1)
            return import_single_whitelist(specs.front(), default_length, verbose);

        std::vector<std::unordered_set<int64_seq>> sets;
        for (auto const &spec : specs) {
            sets.push_back(load_barcodes(spec, default_length, verbose));
        }

        // assemble the wl_entry
        wl_entry out;
        out.global_bcs.clear();
        out.true_bcs.clear();
        populate_entry(out, std::move(sets));
        return out;
    }

};

/**
 * @namespace bc_mem_utils
 * @brief Namespace for memory utility functions related to barcode containers
 */
namespace bc_mem_utils {
// Table allocations are estimates, not process RSS. Native key and entry
// sizes include the platform's alignment; resident keys have no heap vectors.
inline std::size_t get_int64seq_mem(const int64_seq&) { return sizeof(int64_seq); }
inline std::size_t get_bc_entry_mem(const barcode_entry&) { return sizeof(barcode_entry); }
template<class Word>
inline std::size_t get_bc_entry_mem(const basic_barcode_entry<Word>&) {
    return sizeof(basic_barcode_entry<Word>);
}
template<class Map> std::size_t approx_unique(const Map& map) { return map.unique_storage_bytes(); }
template<class Map> std::size_t approx_assoc(const Map& map) { return map.association_storage_bytes(); }

template<class Map>
void print_submap_distribution(const Map& map) {
    std::cout << "Resident barcode widths: uint64=" << map.narrow_size()
              << ", uint128=" << map.wide_size() << '\n';
}

template<class Map>
void print_memory_report(const Map& map, const std::string& name = "") {
    if (!name.empty()) std::cout << "=== " << name << " ===\n";
    const auto unique = approx_unique(map);
    const auto associations = approx_assoc(map);
    std::cout << "Counts:\n"
              << "  Unique values:    " << map.unique_val_size() << '\n'
              << "  Associations:     " << map.association_size() << '\n'
              << "Estimated allocated storage:\n"
              << "  Unique storage:   " << unique / 1024 << " KiB\n"
              << "  Association map:  " << associations / 1024 << " KiB\n"
              << "  Total:            " << (unique + associations) / 1024 << " KiB\n"
              << "Resident entry sizes (counts and flags included):\n"
              << "  uint64:           " << sizeof(basic_barcode_entry<uint64_t>) << " bytes\n"
              << "  uint128:          " << sizeof(basic_barcode_entry<uint128_t>) << " bytes\n";
    print_submap_distribution(map);
    const auto bad_keys = map.validate_association_keys();
    const auto null_values = map.validate_association_values();
    if (bad_keys || null_values)
        std::cout << "Issues: " << bad_keys << " invalid keys, " << null_values << " null values\n";
}

template<class WLTrue, class WLGlobal>
void print_mem_snapshot(const WLTrue& true_bcs, const WLGlobal& global_bcs, size_t iteration) {
    const auto true_bytes = approx_unique(true_bcs) + approx_assoc(true_bcs);
    const auto global_bytes = approx_unique(global_bcs) + approx_assoc(global_bcs);
    std::cout << "\n--- Iteration " << iteration << " Memory Snapshot ---\n"
              << "true_bcs:   " << true_bcs.unique_val_size() << " unique, "
              << true_bcs.association_size() << " associations, " << true_bytes / 1024 / 1024 << " MiB\n"
              << "global_bcs: " << global_bcs.unique_val_size() << " unique, "
              << global_bcs.association_size() << " associations, " << global_bytes / 1024 / 1024 << " MiB\n"
              << "TOTAL estimated storage: " << (true_bytes + global_bytes) / 1024 / 1024 << " MiB\n";
}
} // namespace bc_mem_utils
