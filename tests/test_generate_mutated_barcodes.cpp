#include "rad/rad_headers.h"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using mutation_set = std::unordered_set<int64_seq>;

mutation_set legacy_generate_mutated_barcodes(
    const int64_seq& seq, int mutation_rounds) {
    if (seq.bits.empty()) {
        return {};
    }

    const int sequence_length = seq.length;
    const uint64_t initial_raw = static_cast<uint64_t>(seq.bits[0]);

    size_t expected_total;
    if (mutation_rounds == 1) {
        expected_total = 100;
    } else if (mutation_rounds == 2) {
        expected_total = 2000;
    } else if (mutation_rounds == 3) {
        expected_total = 60000;
    } else {
        expected_total = 1000;
    }

    std::unordered_set<uint64_t> all_mutations_raw;
    all_mutations_raw.reserve(expected_total);
    all_mutations_raw.insert(initial_raw);

    std::unordered_set<uint64_t> current_round_raw;
    current_round_raw.reserve(mutation_rounds == 1 ? 100 : 2000);
    current_round_raw.insert(initial_raw);

    for (int round = 0; round < mutation_rounds; ++round) {
        std::unordered_set<uint64_t> next_round_raw;
        if (round == 0) {
            next_round_raw.reserve(100);
        } else if (round == 1) {
            next_round_raw.reserve(2000);
        } else {
            next_round_raw.reserve(15000);
        }

        for (const uint64_t candidate_raw : current_round_raw) {
            for (int position = 0; position < sequence_length; ++position) {
                const int shift = 2 * position;
                const uint64_t original_nucleotide =
                    (candidate_raw >> shift) & uint64_t{3};
                const uint64_t position_mask = uint64_t{3} << shift;
                for (uint64_t nucleotide = 0; nucleotide < 4; ++nucleotide) {
                    if (nucleotide == original_nucleotide) {
                        continue;
                    }
                    const uint64_t mutated =
                        (candidate_raw & ~position_mask) |
                        (nucleotide << shift);
                    if (all_mutations_raw.insert(mutated).second) {
                        next_round_raw.insert(mutated);
                    }
                }
            }
        }

        if (next_round_raw.empty()) {
            break;
        }
        current_round_raw = std::move(next_round_raw);
    }

    all_mutations_raw.erase(initial_raw);
    mutation_set result;
    result.reserve(all_mutations_raw.size());
    for (const uint64_t raw : all_mutations_raw) {
        result.emplace(static_cast<int64_t>(raw), sequence_length);
    }
    return result;
}

int hamming_distance(uint64_t lhs, uint64_t rhs, int length) {
    uint64_t differences = lhs ^ rhs;
    int distance = 0;
    for (int position = 0; position < length; ++position) {
        distance += ((differences >> (2 * position)) & uint64_t{3}) != 0;
    }
    return distance;
}

mutation_set exhaustive_oracle(const int64_seq& seq, int mutation_rounds) {
    mutation_set result;
    if (seq.bits.empty() || mutation_rounds <= 0 ||
        seq.length == 0 || seq.length > 10) {
        return result;
    }

    const int length = seq.length;
    const int max_distance = std::min(mutation_rounds, length);
    const uint64_t original = static_cast<uint64_t>(seq.bits[0]);
    const uint64_t used_mask =
        length == 32 ? ~uint64_t{0}
                     : (uint64_t{1} << (2 * length)) - 1;
    const uint64_t unused_bits = original & ~used_mask;
    const uint64_t possible_sequences = uint64_t{1} << (2 * length);

    for (uint64_t candidate = 0; candidate < possible_sequences;
         ++candidate) {
        const uint64_t raw = unused_bits | candidate;
        const int distance = hamming_distance(original, raw, length);
        if (distance >= 1 && distance <= max_distance) {
            result.emplace(static_cast<int64_t>(raw), length);
        }
    }
    return result;
}

size_t expected_cardinality(int length, int rounds) {
    if (length <= 0 || rounds <= 0) {
        return 0;
    }

    const int max_distance = std::min(length, rounds);
    size_t combinations = 1;
    size_t substitutions = 1;
    size_t total = 0;
    for (int distance = 1; distance <= max_distance; ++distance) {
        combinations =
            combinations * static_cast<size_t>(length - distance + 1) /
            static_cast<size_t>(distance);
        substitutions *= 3;
        total += combinations * substitutions;
    }
    return total;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void verify_case(const int64_seq& seq, int rounds, bool check_oracle) {
    const mutation_set legacy =
        legacy_generate_mutated_barcodes(seq, rounds);
    const mutation_set direct =
        mutation_tools::generate_mutated_barcodes(seq, rounds);

    require(
        direct == legacy,
        "legacy/direct set mismatch for length=" +
            std::to_string(seq.length) +
            ", rounds=" + std::to_string(rounds));

    if (seq.length <= 32) {
        require(
            direct.size() ==
                expected_cardinality(static_cast<int>(seq.length), rounds),
            "unexpected cardinality for length=" +
                std::to_string(seq.length) +
                ", rounds=" + std::to_string(rounds));
    }

    if (check_oracle) {
        const mutation_set oracle = exhaustive_oracle(seq, rounds);
        require(
            direct == oracle,
            "direct/oracle set mismatch for length=" +
                std::to_string(seq.length) +
                ", rounds=" + std::to_string(rounds));
    }

    if (!seq.bits.empty()) {
        require(
            direct.find(seq) == direct.end(),
            "original barcode was included");
        for (const int64_seq& candidate : direct) {
            require(
                candidate.length == seq.length &&
                    candidate.bits.size() == 1,
                "mutated barcode shape changed");
        }
    }
}

void run_equivalence_tests() {
    verify_case(int64_seq{}, -1, false);
    verify_case(int64_seq{}, 0, false);
    verify_case(int64_seq{}, 2, false);
    verify_case(int64_seq(0, 0), -1, true);
    verify_case(int64_seq(0, 0), 0, true);
    verify_case(int64_seq(0, 0), 2, true);

    // Exhaust every origin through length four and include rounds beyond the
    // sequence length, where the legacy frontier naturally becomes empty.
    for (int length = 1; length <= 4; ++length) {
        const uint64_t sequence_count = uint64_t{1} << (2 * length);
        for (uint64_t raw = 0; raw < sequence_count; ++raw) {
            for (int rounds = -1; rounds <= length + 1; ++rounds) {
                verify_case(
                    int64_seq(static_cast<int64_t>(raw), length),
                    rounds, true);
            }
        }
    }

    // Preserve unused high bits in non-canonical raw inputs.
    verify_case(
        int64_seq(static_cast<int64_t>(0xfedcba9876540000ULL), 8),
        3, true);

    const std::vector<uint64_t> patterns = {
        0x0000000000000000ULL,
        0xffffffffffffffffULL,
        0xaaaaaaaaaaaaaaaaULL,
        0x5555555555555555ULL,
        0x0123456789abcdefULL,
        0xfedcba9876543210ULL
    };
    for (const int length : {8, 14, 16, 31, 32}) {
        for (const uint64_t pattern : patterns) {
            verify_case(
                int64_seq(static_cast<int64_t>(pattern), length),
                1, false);
            verify_case(
                int64_seq(static_cast<int64_t>(pattern), length),
                2, false);
        }
    }

    // Distance three is materially larger, so cover representative production
    // and top-bit inputs without making every CTest invocation allocator-heavy.
    verify_case(int64_seq("ACGTACGTACGTACGT"), 3, false);
    verify_case(
        int64_seq(static_cast<int64_t>(0xffffffffffffffffULL), 32),
        3, false);

    std::mt19937_64 generator(0x5241445f4d555441ULL);
    for (int sample = 0; sample < 128; ++sample) {
        const int length = 1 + static_cast<int>(generator() % 32);
        const int rounds = static_cast<int>(generator() % 3);
        verify_case(
            int64_seq(static_cast<int64_t>(generator()), length),
            rounds, length <= 8);
    }

    require(
        mutation_tools::generate_mutated_barcodes(
            int64_seq("ACGTACGTACGTACGT"), 1).size() == 48,
        "16bp distance-one cardinality must be 48");
    require(
        mutation_tools::generate_mutated_barcodes(
            int64_seq("ACGTACGTACGTACGT"), 2).size() == 1128,
        "16bp distance-two cardinality must be 1128");
    require(
        mutation_tools::generate_mutated_barcodes(
            int64_seq("ACGTACGTACGTACGT"), 3).size() == 16248,
        "16bp distance-three cardinality must be 16248");
}

enum class implementation {
    legacy,
    direct
};

uint64_t benchmark(
    implementation selected, size_t calls,
    const std::vector<int64_seq>& inputs) {
    uint64_t checksum = 0;
    for (size_t index = 0; index < calls; ++index) {
        mutation_set result =
            selected == implementation::legacy
                ? legacy_generate_mutated_barcodes(
                      inputs[index % inputs.size()], 2)
                : mutation_tools::generate_mutated_barcodes(
                      inputs[index % inputs.size()], 2);
        // Consume an order-independent property. Equivalent unordered sets
        // need not have the same begin() element or iteration order.
        checksum =
            checksum * 0x9e3779b97f4a7c15ULL +
            static_cast<uint64_t>(result.size()) + index;
    }
    return checksum;
}

void run_benchmark(const std::string& selected_name, size_t calls) {
    const implementation selected =
        selected_name == "legacy" ? implementation::legacy
                                  : implementation::direct;
    require(
        selected_name == "legacy" || selected_name == "direct",
        "benchmark implementation must be legacy or direct");

    std::mt19937_64 generator(0x5241445f42454e43ULL);
    std::vector<int64_seq> inputs;
    inputs.reserve(256);
    for (int index = 0; index < 256; ++index) {
        // A canonical 16-base packed sequence uses only its low 32 bits.
        inputs.emplace_back(
            static_cast<int64_t>(generator() & 0xffffffffULL), 16);
    }

    // Keep correctness checks and allocator warm-up outside the timed region.
    for (const int64_seq& input : inputs) {
        require(
            legacy_generate_mutated_barcodes(input, 2) ==
                mutation_tools::generate_mutated_barcodes(input, 2),
            "benchmark input failed set-equivalence precondition");
    }
    benchmark(selected, 64, inputs);

    const auto start = std::chrono::steady_clock::now();
    const uint64_t checksum = benchmark(selected, calls, inputs);
    const auto finish = std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration<double>(finish - start).count();

    std::cout << "implementation=" << selected_name
              << " calls=" << calls
              << " candidates_per_call=1128"
              << " seconds=" << seconds
              << " calls_per_second=" << calls / seconds
              << " candidates_per_second="
              << (calls * 1128.0) / seconds
              << " checksum=" << checksum << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--benchmark") {
            run_benchmark(argv[2], std::stoull(argv[3]));
            return 0;
        }
        require(argc == 1, "usage: test_generate_mutated_barcodes "
                           "[--benchmark legacy|direct CALLS]");
        run_equivalence_tests();
        std::cout << "mutation generator equivalence tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "failure: " << error.what() << '\n';
        return 1;
    }
}
