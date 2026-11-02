#pragma once

#include "packed_sequence.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// A lossless prefilter for global whitelist lookup. Short observations keep
// their native 64-bit keys; wide observations have a separate, lazy index.
namespace whitelist_relevance {

using packed_barcode = ::packed_barcode<uint64_t>;
using wide_packed_barcode = ::packed_barcode<uint128_t>;
using packed_hash = std::hash<packed_barcode>;
using barcode_set = std::unordered_set<packed_barcode, packed_hash>;
using wide_barcode_set = std::unordered_set<wide_packed_barcode>;

template <typename Word>
inline bool pack(std::string_view sequence, ::packed_barcode<Word>& out) noexcept {
    out = {};
    if (sequence.empty() || sequence.size() > sizeof(Word) * 4) return false;
    Word bits = 0;
    if (!packed_sequence_detail::encode_word(sequence, bits)) return false;
    out = {bits, static_cast<uint16_t>(sequence.size())};
    return true;
}

template <typename Word>
inline ::packed_barcode<Word> reverse_complement(::packed_barcode<Word> barcode) noexcept {
    Word rc = 0;
    for (unsigned i = 0; i < barcode.length; ++i) {
        rc = (rc << 2) | ((barcode.bits & Word{3}) ^ Word{2});
        barcode.bits >>= 2;
    }
    return {rc, barcode.length};
}

template <typename Word>
inline int hamming_distance(Word a, Word b) noexcept {
    Word different = a ^ b;
    Word alternating = static_cast<Word>(UINT64_C(0x5555555555555555));
    if constexpr (sizeof(Word) > sizeof(uint64_t))
        alternating |= alternating << 64;
    different = (different | (different >> 1)) & alternating;
    int distance = __builtin_popcountll(static_cast<uint64_t>(different));
    if constexpr (sizeof(Word) > sizeof(uint64_t))
        distance += __builtin_popcountll(static_cast<uint64_t>(different >> 64));
    return distance;
}

struct observations {
    barcode_set raw;
    barcode_set exact;
    wide_barcode_set raw_wide;
    wide_barcode_set exact_wide;
    size_t observations_seen = 0;
    size_t unsupported_long_queries = 0;

    template <typename Word, typename Set>
    static void add_typed(std::string_view sequence, std::string_view expanded_window,
                          Set& raw_set, Set& exact_set) {
        ::packed_barcode<Word> barcode;
        if (pack(sequence, barcode)) {
            raw_set.insert(barcode);
            raw_set.insert(reverse_complement(barcode));
            exact_set.insert(barcode);
            exact_set.insert(reverse_complement(barcode));
        }
        const size_t length = sequence.size();
        if (!length) return;
        for (size_t i = 0; i + length <= expanded_window.size(); ++i) {
            if (pack(expanded_window.substr(i, length), barcode)) {
                exact_set.insert(barcode);
                exact_set.insert(reverse_complement(barcode));
            }
        }
    }

    void add(std::string_view sequence, std::string_view expanded_window) {
        ++observations_seen;
        if (sequence.size() > 64) {
            ++unsupported_long_queries;
            return;
        }
        if (sequence.size() <= 32)
            add_typed<uint64_t>(sequence, expanded_window, raw, exact);
        else
            add_typed<uint128_t>(sequence, expanded_window, raw_wide, exact_wide);
    }
    void merge(const observations& other) {
        raw.insert(other.raw.begin(), other.raw.end());
        exact.insert(other.exact.begin(), other.exact.end());
        raw_wide.insert(other.raw_wide.begin(), other.raw_wide.end());
        exact_wide.insert(other.exact_wide.begin(), other.exact_wide.end());
        observations_seen += other.observations_seen;
        unsupported_long_queries += other.unsupported_long_queries;
    }
    size_t raw_size() const noexcept { return raw.size() + raw_wide.size(); }
    size_t exact_size() const noexcept { return exact.size() + exact_wide.size(); }
};
using observation_map = std::unordered_map<std::string, observations>;

class candidate_index {
    template <typename Word>
    struct partition {
        unsigned start = 0, key_length = 0;
        std::vector<size_t> offsets;
        std::vector<Word> words;
    };
    template <typename Word>
    struct length_group {
        unsigned length = 0;
        int max_distance = 0;
        bool all = false;
        std::vector<partition<Word>> partitions;
    };
    const observations* observed_;
    std::array<length_group<uint64_t>, 33> groups_;
    using wide_groups = std::array<length_group<uint128_t>, 65>;
    std::unique_ptr<wide_groups> wide_groups_;

    template <typename Word>
    static uint32_t key(Word word, unsigned length, unsigned start, unsigned count) noexcept {
        const unsigned shift = 2 * (length - start - count);
        return static_cast<uint32_t>((word >> shift) & ((Word{1} << (2 * count)) - Word{1}));
    }

    template <typename Word, size_t N, typename Set>
    static void build_groups(const Set& raw, std::array<length_group<Word>, N>& groups, int cap) {
        std::array<std::vector<Word>, N> words;
        for (const auto& barcode : raw)
            if (barcode.length >= 1 && barcode.length < N)
                words[barcode.length].push_back(barcode.bits);
        for (unsigned length = 1; length < N; ++length) {
            if (words[length].empty() || cap == 0) continue;
            auto& group = groups[length];
            group.length = length;
            group.max_distance = cap;
            if (cap >= static_cast<int>(length)) { group.all = true; continue; }
            const unsigned n_parts = static_cast<unsigned>(cap) + 1;
            unsigned start = 0;
            for (unsigned part = 0; part < n_parts; ++part) {
                const unsigned part_length = length / n_parts + (part < length % n_parts);
                partition<Word> p;
                p.start = start;
                p.key_length = std::min(10U, part_length);
                start += part_length;
                const size_t slots = size_t{1} << (2 * p.key_length);
                p.offsets.assign(slots + 1, 0);
                for (Word word : words[length])
                    ++p.offsets[key(word, length, p.start, p.key_length) + 1];
                for (size_t i = 0; i < slots; ++i) p.offsets[i + 1] += p.offsets[i];
                auto fill = p.offsets;
                p.words.resize(words[length].size());
                for (Word word : words[length])
                    p.words[fill[key(word, length, p.start, p.key_length)]++] = word;
                group.partitions.push_back(std::move(p));
            }
        }
    }

    template <typename Word>
    static bool contains_group(Word bits, uint16_t length, const length_group<Word>& group) noexcept {
        if (group.all) return true;
        for (const auto& p : group.partitions) {
            const uint32_t k = key(bits, length, p.start, p.key_length);
            for (size_t i = p.offsets[k]; i < p.offsets[k + 1]; ++i)
                if (hamming_distance(bits, p.words[i]) <= group.max_distance) return true;
        }
        return false;
    }

    template <typename Word, size_t N>
    static size_t groups_bytes(const std::array<length_group<Word>, N>& groups) noexcept {
        size_t bytes = 0;
        for (const auto& group : groups) {
            bytes += group.partitions.capacity() * sizeof(partition<Word>);
            for (const auto& p : group.partitions)
                bytes += p.offsets.capacity() * sizeof(size_t) + p.words.capacity() * sizeof(Word);
        }
        return bytes;
    }

public:
    candidate_index(const observations& observed, int mutation_distance) : observed_(&observed) {
        const int cap = std::max(0, mutation_distance);
        build_groups(observed.raw, groups_, cap);
        if (!observed.raw_wide.empty() && cap > 0) {
            wide_groups_ = std::make_unique<wide_groups>();
            build_groups(observed.raw_wide, *wide_groups_, cap);
        }
    }

    template <typename Word>
    bool contains(Word bits, uint16_t length) const noexcept {
        if (!length || length > 64) return false;
        if (length <= 32) {
            if constexpr (sizeof(Word) > sizeof(uint64_t))
                if ((static_cast<uint128_t>(bits) >> 64) != 0) return false;
            const uint64_t narrow = static_cast<uint64_t>(bits);
            if (observed_->exact.count({narrow, length})) return true;
            return contains_group(narrow, length, groups_[length]);
        }
        uint128_t wide;
        if constexpr (sizeof(Word) <= sizeof(uint64_t))
            wide = static_cast<uint128_t>(static_cast<uint64_t>(bits));
        else
            wide = static_cast<uint128_t>(bits);
        if (observed_->exact_wide.count({wide, length})) return true;
        return wide_groups_ && contains_group(wide, length, (*wide_groups_)[length]);
    }

    size_t bytes() const noexcept {
        size_t bytes = sizeof(*this) + groups_bytes(groups_);
        if (wide_groups_) bytes += sizeof(wide_groups) + groups_bytes(*wide_groups_);
        return bytes;
    }
};

} // namespace whitelist_relevance
