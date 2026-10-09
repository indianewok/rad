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
static bool brute(const observations& observed, packed_barcode b, int cap) {
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
    packed_barcode acgt;
    require(pack("ACGT", acgt), "packing failed");
    candidate_index invalid_index(invalid, 2);
    require(invalid_index.contains(acgt.bits, acgt.length), "valid shifted kmer next to N was lost");
    require(reverse_complement(reverse_complement(acgt)) == acgt, "RC encoding failed");
    invalid.add(std::string(33, 'A'), std::string(39, 'A'));
    require(invalid.unsupported_long_queries == 1, "long-query fallback not recorded");
    observations empty;
    candidate_index empty_index(empty, 2);
    require(!empty_index.contains(acgt.bits, 4), "empty discovery retained an unrelated key");
    std::cout << "PASS whitelist relevance: " << comparisons << " indexed/exhaustive comparisons\n";
}
