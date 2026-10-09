#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// A lossless prefilter for the *global* whitelist lookup in sigstring.hpp.
// Global lookup reaches raw/RC identities, same-length expanded-window kmers,
// and Hamming neighbours of the raw barcode. The true-whitelist path also
// searches indels exhaustively; small true and dual catalogs remain unfiltered.
// This index is sized by observed barcode sequences, never by catalog size.
namespace whitelist_relevance {

struct packed_barcode {
    uint64_t bits = 0;
    uint16_t length = 0;
    bool operator==(const packed_barcode& o) const noexcept {
        return bits == o.bits && length == o.length;
    }
};
struct packed_hash {
    size_t operator()(const packed_barcode& b) const noexcept {
        uint64_t x = b.bits ^ (uint64_t(b.length) * 0x9e3779b97f4a7c15ULL);
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<size_t>(x ^ (x >> 31));
    }
};
using barcode_set = std::unordered_set<packed_barcode, packed_hash>;

// Match int64_seq's A=00 C=01 T=10 G=11 encoding exactly.
inline bool pack(std::string_view sequence, packed_barcode& out) noexcept {
    out = {};
    if (sequence.empty() || sequence.size() > 32) return false;
    uint64_t bits = 0;
    for (char c : sequence) {
        bits <<= 2;
        switch (c) {
            case 'A': break;
            case 'C': bits |= 1; break;
            case 'T': bits |= 2; break;
            case 'G': bits |= 3; break;
            default: return false;
        }
    }
    out = {bits, static_cast<uint16_t>(sequence.size())};
    return true;
}
inline packed_barcode reverse_complement(packed_barcode b) noexcept {
    uint64_t rc = 0, bits = b.bits;
    for (unsigned i = 0; i < b.length; ++i) {
        rc = (rc << 2) | ((bits & 3U) ^ 2U);
        bits >>= 2;
    }
    return {rc, b.length};
}
inline int hamming_distance(uint64_t a, uint64_t b) noexcept {
    uint64_t different = a ^ b;
    different = (different | (different >> 1)) & 0x5555555555555555ULL;
    return __builtin_popcountll(different);
}

struct observations {
    barcode_set raw;
    barcode_set exact;
    size_t observations_seen = 0;
    size_t unsupported_long_queries = 0;

    // Both arguments are oriented like correct_barcode_lookup: reverse read
    // elements are reverse complemented before adding. Extra RC keys are a
    // harmless superset and protect identity lookup in either orientation.
    void add(std::string_view sequence, std::string_view expanded_window) {
        ++observations_seen;
        if (sequence.size() > 32) {
            ++unsupported_long_queries;
            return;
        }
        packed_barcode p;
        if (pack(sequence, p)) {
            raw.insert(p);
            raw.insert(reverse_complement(p));
            exact.insert(p);
            exact.insert(reverse_complement(p));
        }
        const size_t len = sequence.size();
        if (!len) return;
        for (size_t i = 0; i + len <= expanded_window.size(); ++i) {
            if (pack(expanded_window.substr(i, len), p)) {
                exact.insert(p);
                exact.insert(reverse_complement(p));
            }
        }
    }
    void merge(const observations& other) {
        raw.insert(other.raw.begin(), other.raw.end());
        exact.insert(other.exact.begin(), other.exact.end());
        observations_seen += other.observations_seen;
        unsupported_long_queries += other.unsupported_long_queries;
    }
};
using observation_map = std::unordered_map<std::string, observations>;

class candidate_index {
    struct partition {
        unsigned start = 0, key_length = 0;
        std::vector<size_t> offsets;
        std::vector<uint64_t> words;
    };
    struct length_group {
        unsigned length = 0;
        int max_distance = 0;
        bool all = false;
        std::vector<partition> partitions;
    };
    const observations* observed_;
    std::array<length_group, 33> groups_;

    static uint32_t key(uint64_t word, unsigned length,
                        unsigned start, unsigned count) noexcept {
        const unsigned shift = 2 * (length - start - count);
        return static_cast<uint32_t>((word >> shift) &
                                     ((uint64_t(1) << (2 * count)) - 1));
    }
public:
    candidate_index(const observations& observed, int mutation_distance)
        : observed_(&observed) {
        const int cap = std::max(0, mutation_distance);
        std::array<std::vector<uint64_t>, 33> words;
        for (const auto& b : observed.raw) {
            if (b.length >= 1 && b.length <= 32)
                words[b.length].push_back(b.bits);
        }
        for (unsigned len = 1; len <= 32; ++len) {
            if (words[len].empty() || cap == 0) continue;
            auto& g = groups_[len];
            g.length = len;
            g.max_distance = cap;
            if (cap >= static_cast<int>(len)) { g.all = true; continue; }
            // With <=cap substitutions, at least one of cap+1 disjoint parts
            // is unchanged. Truncating a long part to 10 bases is lossless;
            // confirmation below prevents saturated seeds retaining a row.
            const unsigned n_parts = static_cast<unsigned>(cap) + 1;
            unsigned start = 0;
            for (unsigned part = 0; part < n_parts; ++part) {
                const unsigned part_len = len / n_parts + (part < len % n_parts);
                partition p;
                p.start = start;
                p.key_length = std::min(10U, part_len);
                start += part_len;
                const size_t slots = size_t(1) << (2 * p.key_length);
                p.offsets.assign(slots + 1, 0);
                for (uint64_t word : words[len])
                    ++p.offsets[key(word, len, p.start, p.key_length) + 1];
                for (size_t i = 0; i < slots; ++i) p.offsets[i + 1] += p.offsets[i];
                auto fill = p.offsets;
                p.words.resize(words[len].size());
                for (uint64_t word : words[len])
                    p.words[fill[key(word, len, p.start, p.key_length)]++] = word;
                g.partitions.push_back(std::move(p));
            }
        }
    }
    bool contains(uint64_t bits, uint16_t length) const noexcept {
        if (observed_->exact.count({bits, length})) return true;
        if (!length || length > 32) return false;
        const auto& g = groups_[length];
        if (g.all) return true;
        for (const auto& p : g.partitions) {
            const uint32_t k = key(bits, length, p.start, p.key_length);
            for (size_t i = p.offsets[k]; i < p.offsets[k + 1]; ++i)
                if (hamming_distance(bits, p.words[i]) <= g.max_distance) return true;
        }
        return false;
    }
    size_t bytes() const noexcept {
        size_t n = sizeof(*this);
        for (const auto& g : groups_)
            for (const auto& p : g.partitions)
                n += p.offsets.capacity() * sizeof(size_t) + p.words.capacity() * sizeof(uint64_t);
        return n;
    }
};

} // namespace whitelist_relevance
