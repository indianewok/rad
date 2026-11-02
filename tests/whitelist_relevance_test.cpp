#include "rad/whitelist_relevance.hpp"
#include <cstdlib>
#include <iostream>
#include <random>

using namespace whitelist_relevance;
static void require(bool pass, const char* message) {
    if (!pass) { std::cerr << message << '\n'; std::exit(1); }
}
static std::string dna(uint64_t word, unsigned length) {
    std::string result(length, 'A');
    const char bases[] = {'A', 'C', 'T', 'G'};
    for (unsigned i = 0; i < length; ++i) {
        result[length - 1 - i] = bases[word & 3];
        word >>= 2;
    }
    return result;
}
static bool brute(const observations& observed, whitelist_relevance::packed_barcode b, int cap) {
    if (observed.exact.count(b)) return true;
    for (const auto& query : observed.raw)
        if (query.length == b.length && hamming_distance(query.bits, b.bits) <= cap) return true;
    return false;
}
int main() {
    std::mt19937_64 random(20261009);
    size_t comparisons = 0;
    // Exhaustively compare all possible catalog sequences of short lengths.
    // This exercises saturated seeds, zero cap, all-match cap and length groups.
    for (unsigned length = 1; length <= 8; ++length) {
        observations observed;
        for (unsigned i = 0; i < 12; ++i) {
            const std::string raw = dna(random(), length);
            observed.add(raw, "ACT" + raw + "GTA");
        }
        for (int cap = 0; cap <= 3; ++cap) {
            candidate_index index(observed, cap);
            for (uint64_t bits = 0; bits < (uint64_t(1) << (2 * length)); ++bits) {
                require(index.contains(bits, length) == brute(observed, {bits, static_cast<uint16_t>(length)}, cap),
                        "indexed short-word set differs from exhaustive Hamming/exact union");
                ++comparisons;
            }
        }
    }
    // Real Stereo lengths and 32-base boundary. Inject several mutations,
    // shifted exact windows, reverse complements and distinct lengths.
    observations observed;
    for (unsigned length : {24U, 25U, 26U, 32U}) {
        for (unsigned i = 0; i < 40; ++i) {
            const auto raw = dna(random(), length);
            observed.add(raw, "GAT" + raw + "TGC");
        }
    }
    observations merge_copy;
    merge_copy.merge(observed);
    for (int cap = 0; cap <= 4; ++cap) {
        candidate_index index(merge_copy, cap);
        for (const auto& query : observed.raw) {
            for (int trial = 0; trial < 50; ++trial) {
                uint64_t target = query.bits;
                const int mutations = trial % 6;
                for (int edit = 0; edit < mutations; ++edit)
                    target ^= (uint64_t(1 + random() % 3) << (2 * (random() % query.length)));
                require(index.contains(target, query.length) == brute(observed, {target, query.length}, cap),
                        "indexed real-length set differs from exhaustive Hamming/exact union");
                ++comparisons;
            }
        }
        for (const auto& exact : observed.exact)
            require(index.contains(exact.bits, exact.length), "lost exact shifted-window or RC key");
    }
    observations invalid;
    invalid.add("NCGT", "AACGTN");
    require(invalid.raw.empty(), "ambiguous barcode became an observed raw identity");
    whitelist_relevance::packed_barcode acgt;
    require(pack("ACGT", acgt), "packing failed");
    candidate_index invalid_index(invalid, 2);
    require(invalid_index.contains(acgt.bits, acgt.length), "valid shifted kmer next to N was lost");
    require(reverse_complement(reverse_complement(acgt)) == acgt, "RC encoding failed");
    invalid.add(std::string(65, 'A'), std::string(71, 'A'));
    require(invalid.unsupported_long_queries == 1, "long-query fallback not recorded");
    observations wide;
    for (unsigned length : {33U, 63U, 64U}) {
        for (unsigned i = 0; i < 8; ++i) {
            std::string raw(length, 'A');
            for (char& base : raw) base = "ACTG"[random() & 3];
            wide.add(raw, "GAT" + raw + "TGC");
        }
    }
    require(wide.raw.empty() && !wide.raw_wide.empty(), "wide observations entered the narrow set");
    observations wide_copy;
    wide_copy.merge(wide);
    for (int cap = 0; cap <= 3; ++cap) {
        candidate_index index(wide_copy, cap);
        for (const auto& query : wide.raw_wide) {
            for (int trial = 0; trial < 30; ++trial) {
                uint128_t target = query.bits;
                for (int edit = 0; edit < trial % 6; ++edit)
                    target ^= uint128_t(1 + random() % 3) << (2 * (random() % query.length));
                bool expected = wide.exact_wide.count({target, query.length}) != 0;
                for (const auto& raw : wide.raw_wide)
                    if (raw.length == query.length && hamming_distance(raw.bits, target) <= cap) expected = true;
                require(index.contains(target, query.length) == expected,
                        "wide index differs from exhaustive Hamming/exact union");
                ++comparisons;
            }
        }
        for (const auto& exact : wide.exact_wide)
            require(index.contains(exact.bits, exact.length), "lost wide exact shifted-window or RC key");
    }
    observations empty;
    candidate_index empty_index(empty, 2);
    require(!empty_index.contains(acgt.bits, 4), "empty discovery retained an unrelated key");
    std::cout << "PASS whitelist relevance: " << comparisons << " indexed/exhaustive comparisons\n";
}
