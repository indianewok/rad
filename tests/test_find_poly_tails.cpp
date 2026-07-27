#include "rad/sigstring.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// Frozen copy of the pre-optimization implementation. Keep this deliberately
// straightforward: its purpose is to detect any semantic drift in faster
// implementations of aligner_tools::find_poly_tails.
static_alignments legacy_find_poly_tails(const std::string& query,
                                         const std::string& sequence,
                                         int window_size) {
    static_alignments result;
    result.success = false;
    result.edit_distance = 1;
    if (sequence.length() < static_cast<size_t>(window_size)) {
        return result;
    }
    char poly_base = std::toupper(query[0]);
    int min_count = static_cast<int>(window_size * 0.9);
    int min_gap = 3;
    int i = 0;
    while (i <= static_cast<int>(sequence.length()) - window_size) {
        int count = std::count_if(
            sequence.begin() + i,
            sequence.begin() + i + window_size,
            [poly_base](char c) { return std::toupper(c) == poly_base; }
        );
        if (count >= min_count) {
            int current_start = i;
            int current_end = i + window_size - 1;
            int non_poly_count = 0;
            int last_poly_pos = current_end;
            while (current_end + 1 < static_cast<int>(sequence.length())) {
                if (std::toupper(sequence[current_end + 1]) == poly_base) {
                    if (non_poly_count <= min_gap) {
                        current_end++;
                        last_poly_pos = current_end;
                        non_poly_count = 0;
                    } else {
                        break;
                    }
                } else {
                    current_end++;
                    non_poly_count++;
                    if (non_poly_count > min_gap) {
                        current_end = last_poly_pos;
                        break;
                    }
                }
            }
            result.success = true;
            result.edit_distance = window_size - count;
            result.positions.emplace_back(current_start + 1, current_end + 1);
            i = current_end + 1;
            i += min_gap;
        } else {
            i++;
        }
    }
    return result;
}

std::string describe(const static_alignments& result) {
    std::ostringstream out;
    out << "success=" << std::boolalpha << result.success
        << ", edit_distance=" << result.edit_distance << ", positions=[";
    for (std::size_t i = 0; i < result.positions.size(); ++i) {
        if (i != 0) {
            out << ", ";
        }
        out << '(' << result.positions[i].first << ','
            << result.positions[i].second << ')';
    }
    return out.str() + ']';
}

bool same_result(const static_alignments& lhs,
                 const static_alignments& rhs) {
    return lhs.success == rhs.success &&
           lhs.edit_distance == rhs.edit_distance &&
           lhs.positions == rhs.positions;
}

bool check_equal(const std::string& label,
                 const std::string& query,
                 const std::string& sequence,
                 int window_size,
                 const static_alignments& expected,
                 aligner_tools& aligner) {
    const static_alignments actual =
        aligner.find_poly_tails(query, sequence, window_size);
    if (same_result(actual, expected)) {
        return true;
    }

    std::cerr << "FAIL: " << label << '\n'
              << "  query: " << query << '\n'
              << "  sequence: " << sequence << '\n'
              << "  window_size: " << window_size << '\n'
              << "  expected: " << describe(expected) << '\n'
              << "  actual:   " << describe(actual) << '\n';
    return false;
}

static_alignments expected_result(
    bool success,
    int edit_distance,
    std::vector<std::pair<int, int>> positions = {}) {
    static_alignments result;
    result.success = success;
    result.edit_distance = edit_distance;
    result.positions = std::move(positions);
    return result;
}

bool check_explicit_cases(aligner_tools& aligner) {
    struct TestCase {
        const char* label;
        std::string query;
        std::string sequence;
        int window_size;
        static_alignments expected;
    };

    const std::vector<TestCase> cases = {
        {
            "sequence shorter than window",
            "A{12,}+",
            std::string(13, 'A'),
            14,
            expected_result(false, 1),
        },
        {
            "no matching bases",
            "A",
            std::string(14, 'C'),
            14,
            expected_result(false, 1),
        },
        {
            "eleven of fourteen does not meet truncated threshold",
            "A",
            std::string(11, 'A') + "CCC",
            14,
            expected_result(false, 1),
        },
        {
            "twelve of fourteen meets truncated threshold",
            "A",
            std::string(12, 'A') + "CC",
            14,
            expected_result(true, 2, {{1, 14}}),
        },
        {
            "lowercase bases remain case insensitive",
            "a{12,}+",
            std::string(12, 'a') + "cc",
            14,
            expected_result(true, 2, {{1, 14}}),
        },
        {
            "threshold window can begin inside a non-poly prefix",
            "A",
            "CCC" + std::string(14, 'A'),
            14,
            expected_result(true, 2, {{2, 17}}),
        },
        {
            "three-base internal gap is bridged",
            "A",
            std::string(14, 'A') + "CCC" + std::string(4, 'A'),
            14,
            expected_result(true, 0, {{1, 21}}),
        },
        {
            "four-base internal gap ends at last poly base",
            "A",
            std::string(14, 'A') + "CCCC" + std::string(4, 'A'),
            14,
            expected_result(true, 0, {{1, 14}}),
        },
        {
            "short trailing non-poly run at EOF remains included",
            "A",
            std::string(14, 'A') + "CCC",
            14,
            expected_result(true, 0, {{1, 17}}),
        },
        {
            "multiple tails preserve legacy skip and coordinates",
            "A",
            std::string(14, 'A') + std::string(8, 'C') +
                std::string(14, 'A'),
            14,
            expected_result(true, 2, {{1, 14}, {21, 36}}),
        },
        {
            "poly-T query uses first query character",
            "T{12,}+",
            "GG" + std::string(12, 't'),
            14,
            expected_result(true, 2, {{1, 14}}),
        },
    };

    bool ok = true;
    for (const auto& test : cases) {
        ok = check_equal(test.label, test.query, test.sequence,
                         test.window_size, test.expected, aligner) && ok;
    }
    return ok;
}

bool check_exhaustive_binary_sequences(aligner_tools& aligner) {
    constexpr int window_size = 14;
    constexpr int max_length = 18;

    for (int length = 0; length <= max_length; ++length) {
        const std::uint64_t sequence_count = std::uint64_t{1} << length;
        for (std::uint64_t mask = 0; mask < sequence_count; ++mask) {
            std::string sequence(static_cast<std::size_t>(length), 'C');
            for (int position = 0; position < length; ++position) {
                if ((mask & (std::uint64_t{1} << position)) != 0) {
                    sequence[static_cast<std::size_t>(position)] = 'A';
                }
            }

            const static_alignments expected =
                legacy_find_poly_tails("A{12,}+", sequence, window_size);
            if (!check_equal("exhaustive A/C sequence", "A{12,}+",
                             sequence, window_size, expected, aligner)) {
                std::cerr << "  exhaustive mask: " << mask
                          << ", length: " << length << '\n';
                return false;
            }
        }
    }
    return true;
}

bool check_post_hit_reseed_cases(aligner_tools& aligner) {
    const std::vector<int> window_sizes = {2, 3, 5, 8, 14, 17, 32};
    const std::vector<int> gap_lengths = {3, 4, 5, 6, 9};

    for (const int window_size : window_sizes) {
        for (int prefix_length = 0; prefix_length <= 5; ++prefix_length) {
            for (const int gap_length : gap_lengths) {
                for (int tail_extension = 0; tail_extension <= 3;
                     ++tail_extension) {
                    std::string sequence(
                        static_cast<std::size_t>(prefix_length), 'C');
                    sequence += std::string(
                        static_cast<std::size_t>(
                            window_size + tail_extension),
                        'A');
                    sequence += std::string(
                        static_cast<std::size_t>(gap_length), 'C');
                    sequence += std::string(
                        static_cast<std::size_t>(window_size + 3), 'A');
                    sequence += "CCCC";
                    sequence += std::string(
                        static_cast<std::size_t>(window_size + 1), 'A');

                    // Exercise the same re-seeding positions with mixed-case
                    // input as well as the usual uppercase FASTQ alphabet.
                    for (std::size_t i = 0; i < sequence.size(); ++i) {
                        if ((i % 2) != 0) {
                            sequence[i] = static_cast<char>(
                                std::tolower(
                                    static_cast<unsigned char>(sequence[i]))
                            );
                        }
                    }

                    const static_alignments expected =
                        legacy_find_poly_tails(
                            "a{12,}+", sequence, window_size);
                    if (!check_equal(
                            "post-hit rolling-window re-seed",
                            "a{12,}+",
                            sequence,
                            window_size,
                            expected,
                            aligner)) {
                        std::cerr
                            << "  prefix length: " << prefix_length
                            << ", gap length: " << gap_length
                            << ", tail extension: " << tail_extension
                            << '\n';
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool check_deterministic_mixed_corpus(aligner_tools& aligner) {
    constexpr std::uint64_t seed = 0x524144504f4c5954ULL;
    constexpr int cases = 5000;
    const std::string alphabet = "ACGTNacgtn";
    const std::vector<std::string> queries = {
        "A", "C", "G", "T", "a{12,}+", "t{12,}+",
    };

    std::mt19937_64 random(seed);
    std::uniform_int_distribution<int> length_distribution(0, 512);
    std::uniform_int_distribution<int> window_distribution(1, 32);
    std::uniform_int_distribution<std::size_t> base_distribution(
        0, alphabet.size() - 1);
    std::uniform_int_distribution<std::size_t> query_distribution(
        0, queries.size() - 1);

    for (int case_index = 0; case_index < cases; ++case_index) {
        const int length = length_distribution(random);
        const int window_size = window_distribution(random);
        const std::string& query = queries[query_distribution(random)];
        std::string sequence(static_cast<std::size_t>(length), 'N');
        for (char& base : sequence) {
            base = alphabet[base_distribution(random)];
        }

        const static_alignments expected =
            legacy_find_poly_tails(query, sequence, window_size);
        if (!check_equal("fixed-seed mixed corpus", query, sequence,
                         window_size, expected, aligner)) {
            std::cerr << "  corpus seed: " << seed
                      << ", case index: " << case_index << '\n';
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    aligner_tools aligner;

    if (!check_explicit_cases(aligner)) {
        return 1;
    }
    if (!check_exhaustive_binary_sequences(aligner)) {
        return 1;
    }
    if (!check_post_hit_reseed_cases(aligner)) {
        return 1;
    }
    if (!check_deterministic_mixed_corpus(aligner)) {
        return 1;
    }

    std::cout << "find_poly_tails regression test passed\n";
    return 0;
}
