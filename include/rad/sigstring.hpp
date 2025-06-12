#pragma once
#include "rad_headers.h"

/**
 * @brief Represents a processed element with alignment information
 */
struct seq_element {
    std::string class_id;               // Unique identifier for the element
    std::string global_class;           // Global classification
    std::optional<int> edit_distance;    //Edit distance from alignment
    std::pair<int, int> position;       // Start and stop positions
    std::string type;                   // Element type
    int order;                          // Position in layout
    std::string direction;              // Orientation
    std::string flags;
    std::optional<bool> element_pass;  // Whether element passed validation
    std::optional<std::string> seq;    // Element sequence if available
    std::optional<std::string> qual;   // Element quality scores if available
    std::optional<std::string> original_seq;  // corrected sequence if available


    seq_element(
        std::string class_id,
        std::string global_class,
        std::optional<int> edit_distance,
        std::pair<int, int> position,
        std::string type,
        int order,
        std::string flags,
        std::string direction,
        std::optional<bool> element_pass = std::nullopt,
        std::optional<std::string> seq = std::nullopt,
        std::optional<std::string> qual = std::nullopt,
        std::optional<std::string> original_seq = std::nullopt
    ) : class_id(std::move(class_id)),
        global_class(std::move(global_class)),
        edit_distance(edit_distance),
        position(position),
        type(std::move(type)),
        order(order),
        flags(std::move(flags)),
        direction(std::move(direction)),
        element_pass(element_pass),
        seq(std::move(seq)),
        qual(std::move(qual)),
        original_seq(std::move(original_seq)) {}
};

/**
 * @brief Represents a static processed element with alignment information
 */
struct static_alignments {
    std::vector<std::pair<int, int>> positions;
    int edit_distance;
    bool success;
    std::string seq;
    std::string cigar; 
    int score;
    std::pair<int,int> pos; 
};

// Tags for the multi-index container
struct sig_id_tag {};
struct sig_global_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};
struct sig_pass_tag {};
struct sig_read_tag {};

// Define the multi-index container
typedef boost::multi_index::multi_index_container<
    seq_element,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_id_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::class_id>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_global_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::global_class>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_ed_tag>,
            boost::multi_index::member<seq_element, std::optional<int>, &seq_element::edit_distance>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_order_tag>,
            boost::multi_index::member<seq_element, int, &seq_element::order>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_dir_tag>,
            boost::multi_index::member<seq_element, std::string, &seq_element::direction>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_pass_tag>,
            boost::multi_index::member<seq_element, std::optional<bool>, &seq_element::element_pass>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<sig_read_tag>,
            boost::multi_index::composite_key<
                seq_element,
                boost::multi_index::member<seq_element, std::string, &seq_element::direction>,
                boost::multi_index::member<seq_element, std::optional<bool>, &seq_element::element_pass>,
                boost::multi_index::member<seq_element, std::pair<int,int>, &seq_element::position>
            >
        >
    >
> SigElement;

class aligner_tools {
private:
    int get_max_consecutive_matches(const char* cigar) {
        int max_matches = 0;
        int current_matches = 0;
        int number = 0;
        for (const char* p = cigar; *p; p++) {
            if (std::isdigit(*p)) {
                number = number * 10 + (*p - '0');
            } else {
                // For simplicity, treat 'M' and '=' as matches.
                if (*p == 'M' || *p == '=') {
                    current_matches = number;
                } else {
                    current_matches = 0;
                }
                if (current_matches > max_matches)
                    max_matches = current_matches;
                number = 0;
            }
        }
        return max_matches;
    }
    
    int get_max_consecutive_matches_with_indels(const char* cigar) {
        // First, parse the cigar string into segments.
        std::vector<std::pair<int, char>> segments;
        while (*cigar) {
            int count = 0;
            while (isdigit(*cigar)) {
                count = count * 10 + (*cigar - '0');
                ++cigar;
            }
            char op = *cigar;
            ++cigar; // skip the op
            segments.push_back({count, op});
        }
        
        int best = 0;
        // For each segment that is a match, consider it as a window start.
        for (size_t i = 0; i < segments.size(); i++) {
            // Only start at a match segment.
            if (!(segments[i].second == '=' || segments[i].second == 'M'))
                continue;
            int current = 0;
            bool used_indel = false;
            // Now extend the window from i forward.
            for (size_t j = i; j < segments.size(); j++) {
                char op = segments[j].second;
                int count = segments[j].first;
                if (op == '=' || op == 'M') {
                    current += count;
                } else if ((op == 'I' || op == 'D' || op == 'X') && count <= 2) {
                    if (!used_indel) {
                        used_indel = true;
                    } else {
                        // Already used an indel; break out.
                        break;
                    }
                } else {
                    // For any clipping (S/H) or other operation, break.
                    break;
                }
                best = std::max(best, current);
            }
        }
        return best;
    }

    int compute_edit_distance(const char* cigar) {
        int edit_distance = 0;
        int number = 0;
        while (*cigar) {
            if (std::isdigit(*cigar)) {
                number = number * 10 + (*cigar - '0');
            } else {
                // In an extended CIGAR:
                // '=' represents a match (ignored),
                // 'X' represents a mismatch,
                // 'I' represents an insertion,
                // 'D' represents a deletion.
                if (*cigar == 'X' || *cigar == 'I' || *cigar == 'D') {
                    edit_distance += number;
                }
                number = 0;
            }
            ++cigar;
        }
        return edit_distance;
    }

public:   
    static_alignments align_static_elements(const std::string& query, const std::string& target, bool verbose,
                                        int max_edit_distance = -1, const std::string& masked_query = "", 
                                        bool primary = true, int expected_start = 1, int expected_end = -1) {

        static_alignments alignment;
        alignment.success = false;
        alignment.edit_distance = -1;
        // Our threshold for a “good” alignment
        const int min_match_bases = 5;
        
        // If expected_end is not provided, use the target's length
        if (expected_end < 0)
            expected_end = target.size();

        // Custom equality rules:
        // Uppercase A, C, T, G match their lowercase forms and also the pad character 'x'.
        // The masked letter N matches uppercase A, C, T, G and 'x', but not lowercase.
        const int numEq = 13;
        EdlibEqualityPair customEqualities[numEq] = {
            {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
            {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
            {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'},
            {'N','x'}
        };

        std::string query_to_use = query;

        // --- Phase 1: Run Edlib to obtain candidate intervals ---
        EdlibAlignConfig config = edlibNewAlignConfig(max_edit_distance,
                                                    EDLIB_MODE_HW,
                                                    EDLIB_TASK_LOC,
                                                    customEqualities, numEq);

        EdlibAlignResult edlibResult = edlibAlign(query_to_use.c_str(), query_to_use.size(),
                                                target.c_str(), target.size(),
                                                config);

        std::vector<std::pair<int,int>> candidates;
        if (edlibResult.status == EDLIB_STATUS_OK &&
            edlibResult.numLocations > 0 &&
            edlibResult.editDistance > -1) {
            std::vector<std::pair<int,int>> intervals;
            for (int i = 0; i < edlibResult.numLocations; ++i) {
                int start = edlibResult.startLocations[i] + 1;  // Convert to 1-indexed
                int end   = edlibResult.endLocations[i] + 1;
                intervals.push_back({start, end});
            }
            std::sort(intervals.begin(), intervals.end(),
                    [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
                        return a.first < b.first;
                    });
            // Collapse overlapping intervals.
            if (!intervals.empty()) {
                std::pair<int,int> current = intervals[0];
                for (size_t i = 1; i < intervals.size(); ++i) {
                    if (intervals[i].first <= current.second) {
                        current.second = std::max(current.second, intervals[i].second);
                    } else {
                        candidates.push_back(current);
                        current = intervals[i];
                    }
                }
                candidates.push_back(current);
            }
        }

        // --- Early Out: If Edlib returned exactly one candidate and it is within the expected region, use it.
        if (candidates.size() == 1) {
            auto cand = candidates.front();
            if (cand.first >= expected_start && cand.second <= expected_end) {
                alignment.positions.push_back(cand);
                alignment.success = true;
                alignment.edit_distance = edlibResult.editDistance;
                // Optionally set alignment.seq, alignment.cigar, etc.
                alignment.seq = target.substr(cand.first - 1, cand.second - cand.first + 1);
                // Return immediately without SSW
                edlibFreeAlignResult(edlibResult);
                return alignment;
            }
        }
        edlibFreeAlignResult(edlibResult);
        // --- Phase 2: Run SSW on the entire target (for no candidate or multiple candidates) ---
        // Prepare uppercase strings for SSW.
        std::string query_ssw = query;
        std::string target_ssw = target;
        std::transform(query_ssw.begin(), query_ssw.end(), query_ssw.begin(), ::toupper);
        std::transform(target_ssw.begin(), target_ssw.end(), target_ssw.begin(), ::toupper);

        int32_t maskLen = static_cast<int32_t>(query_ssw.size() / 2);
        if (maskLen < 15)
            maskLen = 15;
        //defaults of match/mismatch/gap_open/gap_extend is 2/2/3/1
        //thoughts for a shorter sequence--i'm okay with a gap being opened, but the longer the gap goes the 
        //more it should cost. so if we score gap open as a mismatch, but then gap extend as a mismatch as well?
        //playing with 2,2,3,2 and 2,2,3,1
        StripedSmithWaterman::Aligner ssw_aligner(2, 2, 3, 2);
        StripedSmithWaterman::Filter ssw_filter;
        ssw_filter.report_cigar = true;
        StripedSmithWaterman::Alignment sswAlign;
        ssw_aligner.Align(target_ssw.c_str(), query_ssw.c_str(), query_ssw.size(),
                        ssw_filter, &sswAlign, maskLen);

        if (sswAlign.query_begin >= 0 && sswAlign.query_end >= sswAlign.query_begin) {
            int ssw_start = sswAlign.query_begin + 1; // Convert to 1-indexed
            int ssw_end = sswAlign.query_end + 1;
            int ssw_length = ssw_end - ssw_start + 1;
            int ssw_max_matches_ind = get_max_consecutive_matches_with_indels(sswAlign.cigar_string.c_str());
            int ssw_max_matches = get_max_consecutive_matches(sswAlign.cigar_string.c_str());
            // Calculate deviation from expected region.
            int deviation = 0;
            if (ssw_start < expected_start)
                deviation += (expected_start - ssw_start);
            if (ssw_end > expected_end)
                deviation += (ssw_end - expected_end);

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "SSW alignment: region ["
                            << ssw_start << ", " << ssw_end << "], length = "
                            << ssw_length << ", max_matches with indels = " 
                            << ssw_max_matches_ind
                            << " , max_matches = " << ssw_max_matches
                            << ", deviation = " << deviation
                            << ", cigar: " << sswAlign.cigar_string << "\n";
                    std::cout << oss.str();
                }
            }
            // Accept candidate if it meets threshold. first case is within the deviation allowed. second case is for concats.
            if (
            (ssw_max_matches > min_match_bases && deviation <= 100 && ssw_length >= 10 && ssw_max_matches_ind >= 8 && ssw_max_matches >= 8)||
            (deviation > 100 && ssw_length >= 10 && ssw_max_matches >= 10)) {
                alignment.positions.push_back({ssw_start, ssw_end});
                alignment.edit_distance = compute_edit_distance(sswAlign.cigar_string.c_str());
                alignment.success = alignment.edit_distance < max_edit_distance + 1 ? true : false;
                alignment.cigar = sswAlign.cigar_string;
                alignment.seq = target.substr(ssw_start - 1, ssw_end - ssw_start + 1);
            }
        }
    return alignment;
}

    static_alignments find_poly_tails(const std::string& query, const std::string& sequence, int window_size) {
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
                // Move i to just after the end of current tail
                i = current_end + 1;
                // Look for next potential poly-tail after a gap
                i += min_gap;
            } else {
                i++;
            }
        }
            return result;
    }
};

namespace barcode_correction {
    // barcode correction result(s)
    // Count-based quality check for barcodes, looking at total summated counts to curb against false positives and spurious corrections
    bool passes_quality_check(const int64_seq& candidate, const whitelist::wl_entry* wl, 
                              const std::string& whitelist_source, bool verbose, bool raw_match = false) {
        
            // Get the appropriate whitelist reference
            const auto& whitelist = (whitelist_source == "global") ? wl->global_bcs : wl->true_bcs;
            
            int raw_count = whitelist.get_bc_count(candidate, barcode_counts::raw);
            int total_count = whitelist.get_bc_count(candidate, barcode_counts::total);
            int filtered_count = whitelist.get_bc_count(candidate, barcode_counts::filtered);
            int corrected_count = whitelist.get_bc_count(candidate, barcode_counts::corrected);
        
            // KEY FIX: If this barcode has no history (total_count = 0), always pass
            // This handles the case where we're correcting to a barcode that doesn't exist yet
            if (total_count == 0) {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Quality check for " << candidate.bits_to_sequence()
                            << " [" << whitelist_source << "]: total=0 (new barcode) -> PASS" << std::endl;
                        std::cout << oss.str();
                    }
                }
                return true;
            }

            // For barcodes with sufficient history, check correction ratio
            double correction_ratio = static_cast<double>(corrected_count) / total_count;
            // Quality checks
            bool low_filtered = filtered_count <= raw_count;
            bool correction_limit = true;
            if(total_count > 10){
                correction_limit = correction_ratio <= 0.5;  // Max 50% corrected
            }
            bool overall_pass = correction_limit;
        return overall_pass;
    }
    
    // Resolves multiple hits for a barcode candidate against a set of candidates, fails if there is no true winner
    std::pair<std::optional<int64_seq>, std::optional<int>> resolve_multiple_hits(
                    const int64_seq& query, const std::unordered_set<int64_seq>& candidates, 
                    int max_dist, bool verbose,const std::string& wl_type, const whitelist::wl_entry* wl = nullptr) {

            auto sorted = mutation_tools::int64_lvdist(query, candidates, max_dist);
            const auto& target_whitelist = (wl_type == "global") ? wl->global_bcs : wl->true_bcs;

            if (sorted.empty()) {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "NO_MATCH\n";
                        std::cout << oss.str();
                    }
                }
                return {std::nullopt, std::nullopt};
            }
            // Get the lowest edit distance

            auto it = sorted.begin();
            int best_dist = it->first;
            const auto& best_set = it->second;

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Best edit distance = " << best_dist << ", # of candidates = " << best_set.size() << "\n";
                    std::cout << oss.str();
                }
            }
            if (best_set.size() == 1) {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Best candidate found: " << best_set.begin()->bits_to_sequence() << "\n";
                        std::cout << oss.str();
                    }
                }
                return {*best_set.begin(), best_dist};
            }
            
            // Multiple candidates with same edit distance
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Ambiguous best hits found: " << best_set.size() 
                        << " candidates with edit distance " << best_dist;
                    if (wl != nullptr) {
                        oss << ". Attempting count-based tie breaking...";
                    }
                    oss << "\n";
                    std::cout << oss.str();
                }
            }
            
            // If whitelist is provided, attempt count-based tie breaking
            if (wl != nullptr) {
                struct candidate_info {
                    int64_seq barcode;
                    int raw_count;
                    int corrected_count;
                    int concat_count;
                    double score;
                };

                std::vector<candidate_info> candidate_infos;
                candidate_infos.reserve(best_set.size());
                // Collect count information for each candidate
                for (const auto& candidate : best_set) {
                    candidate_info info;
                    info.barcode = candidate;
                    info.raw_count = target_whitelist.get_bc_count(candidate, barcode_counts::raw);
                    info.corrected_count = target_whitelist.get_bc_count(candidate, barcode_counts::corrected);
                    int forw_concat = target_whitelist.get_bc_count(candidate, barcode_counts::forw_concat);
                    int rev_concat = target_whitelist.get_bc_count(candidate, barcode_counts::rev_concat);
                    info.concat_count = forw_concat + rev_concat;
                    // Calculate ranking score: Higher raw counts = better, lower corrected/concat = better
                    info.score = (info.raw_count * 2.0) - (info.corrected_count * 1.0) - (info.concat_count * 0.5);
                    candidate_infos.push_back(info);
                    
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "  Candidate: " << candidate.bits_to_sequence() 
                                << " [" << wl_type << "] raw = " << info.raw_count 
                                << " corrected = " << info.corrected_count 
                                << " concat = " << info.concat_count 
                                << " score = " << info.score << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }

                // Sort by score (highest first), then by raw count as tiebreaker
                std::sort(candidate_infos.begin(), candidate_infos.end(), 
                        [](const candidate_info& a, const candidate_info& b) {
                            if (a.score != b.score) {
                                return a.score > b.score;
                            }
                            if (a.raw_count != b.raw_count) {
                                return a.raw_count > b.raw_count;
                            }
                            if (a.corrected_count != b.corrected_count) {
                                return a.corrected_count < b.corrected_count;
                            }
                            return a.concat_count < b.concat_count;
                        });

                // Check if we have a clear winner
                const auto& winner = candidate_infos[0];
                bool has_clear_winner = true;
                if (candidate_infos.size() > 1 && candidate_infos[1].score == winner.score) {
                    has_clear_winner = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Multiple candidates with same top score (" << winner.score 
                                << "), cannot resolve tie" << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
                if (has_clear_winner) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Count-based tie breaking resolved. Winner: " 
                                << winner.barcode.bits_to_sequence() 
                                << " [score=" << winner.score 
                                << " raw=" << winner.raw_count 
                                << " corrected=" << winner.corrected_count 
                                << " concat=" << winner.concat_count << "]" << std::endl;
                            std::cout << oss.str();
                        }
                    }
                    return {winner.barcode, best_dist};
                }

                // Count-based tie breaking failed
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Could not resolve ambiguous hits using count information" << std::endl;
                        std::cout << oss.str();
                    }
                }
            } else {
                // No whitelist provided - show candidates as before
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        for (const auto& s : best_set) {
                            oss << "  " << s.bits_to_sequence() << "\n";
                        }
                        std::cout << oss.str();
                    }
                }
            }
    
        // Return the minimum distance even if we couldn't resolve the ambiguity
        return {std::nullopt, best_dist};
    }

// Enhanced correct_barcode function with quality checking
    std::optional<int64_seq> correct_barcode(const seq_element& elem, const ReadLayout& layout, bool verbose, int mut_dist, int shift_dist) {
        // pick the right whitelist
        auto key = seq_utils::remove_rc(elem.class_id);
        auto &wl = layout.wl_map.maps.at(key).get();
        int max_dist = 4;
        
        // extract and reverse‐complement the raw string
        // encoded barcode and reverse complement
        std::string raw = elem.seq.value();
        if (elem.direction == "reverse") {
            raw = seq_utils::revcomp(raw);
        }
        int64_seq bc, rc_bc;;
        bc.sequence_to_bits(raw);
        rc_bc.sequence_to_bits(seq_utils::revcomp(raw));
        int bc_len = static_cast<int>(bc.length);
        
        // === filtering for messy barcodes ===
        bool filtered_hit = !wl.filter_bcs.empty() && (wl.filter_bcs.check_wl_for(bc) || wl.filter_bcs.check_wl_for(rc_bc));
        if (filtered_hit || seq_utils::int_kmerize(raw, 2) < 4) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "FILTER_CHECK_FOUND\nNO_CHECK_WORKED\n";
                    std::cout << oss.str();
                }
            }
            return std::nullopt;
        }

        // === exact match in global whitelist ===
        if (!wl.global_bcs.empty() && wl.global_bcs.check_wl_for(bc)) {
            auto matched = wl.global_bcs.return_putative_correct_bcs(bc);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "GLOBAL_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq candidate = *matched.begin();
                if (passes_quality_check(candidate, &wl, "global", verbose) || bc == candidate || rc_bc == candidate) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_WORKED\n";
                            std::cout << oss.str();
                        }
                    }
                    return candidate;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_QUALITY_FAILED - falling through to true whitelist\n";
                            std::cout << oss.str();
                        }
                    }
                    // Fall through to true whitelist check
                }
            } else {
                // Use enhanced resolution with count-based tie breaking for global
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose, "global", &wl);
                if (resolved.has_value()) {
                    // Quality check the resolved candidate
                    if (passes_quality_check(resolved.value(), &wl, "global", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_CHECK_MULTIPLE_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_CHECK_MULTIPLE_RESOLVED_QUALITY_FAILED - falling through to true whitelist\n";
                                std::cout << oss.str();
                            }
                        }
                        // Fall through to true whitelist check
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "GLOBAL_CHECK_MULTIPLE_UNRESOLVED";
                            if (min_dist.has_value()) {
                                oss << " (min_dist=" << min_dist.value() << ")";
                            }
                            oss << " - falling through to true whitelist\n";
                            std::cout << oss.str();
                        }
                    }
                    // Fall through to true whitelist check
                }
            }
        }
        
        // === exact match in true barcodes ===
        if (!wl.true_bcs.empty() && wl.true_bcs.check_wl_for(bc)) {
            auto matched = wl.true_bcs.return_putative_correct_bcs(bc);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "ORIGINAL_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq candidate = *matched.begin();
                // Quality check: if this fails on true whitelist, we reject the barcode entirely
                if (passes_quality_check(candidate, &wl, "true", verbose) || bc == candidate || rc_bc == candidate) {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_CHECK_WORKED\n";
                            std::cout << oss.str();
                        }
                    }
                    return candidate;
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_CHECK_QUALITY_FAILED - rejecting barcode\n";
                            std::cout << oss.str();
                        }
                    }
                    return std::nullopt; // Fail the barcode
                }
            } else {
                // Use enhanced resolution with count-based tie breaking for true
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose, "true", &wl);
                if (resolved.has_value()) {
                    // Quality check the resolved candidate
                    if (passes_quality_check(resolved.value(), &wl, "true", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "ORIGINAL_MATCH_COLLISION_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "ORIGINAL_MATCH_COLLISION_RESOLVED_QUALITY_FAILED - rejecting barcode\n";
                                std::cout << oss.str();
                            }
                        }
                        return std::nullopt; // Fail the barcode
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "ORIGINAL_MATCH_COLLISION_UNRESOLVED";
                            if (min_dist.has_value()) {
                                oss << " (min_dist=" << min_dist.value() << ")";
                            }
                            oss << "\n";
                            std::cout << oss.str();
                        }
                    }
                    return std::nullopt; // Fail the barcode
                }
            }
        }

        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "ORIGINAL_MATCH_NOT_FOUND\n";
                std::cout << oss.str();
            }
        }
        
        // === mutation match - check global first, then true ===
        auto muts = mutation_tools::generate_mutated_barcodes(bc, mut_dist);
        // Check mutations against global whitelist first
        if (!wl.global_bcs.empty() && wl.global_bcs.check_wl_for(muts)) {
            auto matched = wl.global_bcs.return_putative_correct_bcs(muts);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "GLOBAL_MUTATION_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq putative_candidate = *matched.begin();
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "GLOBAL_MUTATION_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\n";
                        std::cout << oss.str();
                    }
                }
                int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_dist);
                if (res >= 0) {
                    // Quality check
                    if (passes_quality_check(putative_candidate, &wl, "global", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "LVDIST::" << res << "\nGLOBAL_MUTATION_CHECK_MATCHED\n";
                                std::cout << oss.str();
                            }
                        }
                        return putative_candidate;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_MUTATION_CHECK_QUALITY_FAILED - falling through to true whitelist\n";
                                std::cout << oss.str();
                            }
                        }
                        // Fall through to true whitelist
                    }
                }
            } else {
                // Use enhanced resolution for global mutations
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose, "global", &wl);
                if (resolved.has_value()) {
                    if (passes_quality_check(resolved.value(), &wl, "global", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_MUTATION_MULTIPLE_MATCHED_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - falling through to true whitelist\n";
                                std::cout << oss.str();
                            }
                        }
                        // Fall through to true whitelist
                    }
                }
            }
        }
        // Check mutations against true barcodes
        if (!wl.true_bcs.empty() && wl.true_bcs.check_wl_for(muts)) {
            auto matched = wl.true_bcs.return_putative_correct_bcs(muts);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "MUTATION_CHECK_FOUND (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq putative_candidate = *matched.begin();
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "MUTATION_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\n";
                        std::cout << oss.str();
                    }
                }
                int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_dist);
                if (res >= 0) {
                    // Quality check - if fails on true whitelist, reject entirely
                    if (passes_quality_check(putative_candidate, &wl, "true", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "LVDIST::" << res << "\nMUTATION_CHECK_MATCHED\n";
                                std::cout << oss.str();
                            }
                        }
                        return putative_candidate;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "MUTATION_CHECK_QUALITY_FAILED - rejecting barcode\n";
                                std::cout << oss.str();
                            }
                        }
                        return std::nullopt;
                    }
                }
            } else {
                // Use enhanced resolution for true mutations
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose, "true", &wl);
                if (resolved.has_value()) {
                    if (passes_quality_check(resolved.value(), &wl, "true", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "MUTATION_MULTIPLE_MATCHED_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "MUTATION_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode\n";
                                std::cout << oss.str();
                            }
                        }
                        return std::nullopt;
                    }
                }
            }
        }
        
        // === shift fuzzy match - check global first, then true ===
        auto shifts = mutation_tools::generate_shifted_barcodes(bc, shift_dist);
        // Check shifts against global whitelist first
        if (!wl.global_bcs.empty() && wl.global_bcs.check_wl_for(shifts)) {
            auto matched = wl.global_bcs.return_putative_correct_bcs(shifts);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "GLOBAL_SHIFT_CHECK_WORKED (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq putative_candidate = *matched.begin();
                int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_dist);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "GLOBAL_SHIFT_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\nLVDIST::" << res << "\n";
                        std::cout << oss.str();
                    }
                }
                if (res >= 0) {
                    // Quality check
                    if (passes_quality_check(putative_candidate, &wl, "global", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_SHIFT_CHECK_MATCHED\n";
                                std::cout << oss.str();
                            }
                        }
                        return putative_candidate;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_SHIFT_CHECK_QUALITY_FAILED - falling through to true whitelist\n";
                                std::cout << oss.str();
                            }
                        }
                        // Fall through to true whitelist
                    }
                }
            } else {
                // Use enhanced resolution for global shifts
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose, "global", &wl);
                if (resolved.has_value()) {
                    if (passes_quality_check(resolved.value(), &wl, "global", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_SHIFT_MULTIPLE_MATCHED_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "GLOBAL_SHIFT_MULTIPLE_RESOLVED_QUALITY_FAILED - falling through to true whitelist\n";
                                std::cout << oss.str();
                            }
                        }
                        // Fall through to true whitelist
                    }
                }
            }
        }
        // Check shifts against true barcodes
        if (!wl.true_bcs.empty() && wl.true_bcs.check_wl_for(shifts)) {
            auto matched = wl.true_bcs.return_putative_correct_bcs(shifts);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "SHIFT_CHECK_WORKED (" << matched.size() << " candidates)\n";
                    std::cout << oss.str();
                }
            }
            
            if (matched.size() == 1) {
                int64_seq putative_candidate = *matched.begin();
                int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_dist);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "SHIFT_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\n";
                        oss << "LVDIST::" << res << "\n";
                        std::cout << oss.str();
                    }
                }
                if (res >= 0) {
                    // Quality check - if fails on true whitelist, reject entirely
                    if (passes_quality_check(putative_candidate, &wl, "true", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "SHIFT_CHECK_MATCHED\n";
                                std::cout << oss.str();
                            }
                        }
                        return putative_candidate;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "SHIFT_CHECK_QUALITY_FAILED - rejecting barcode\n";
                                std::cout << oss.str();
                            }
                        }
                        return std::nullopt;
                    }
                }
            } else {
                // Use enhanced resolution for true shifts
                auto [resolved, min_dist] = resolve_multiple_hits(bc, matched, max_dist, verbose,"true", &wl);
                if (resolved.has_value()) {
                    if (passes_quality_check(resolved.value(), &wl, "true", verbose)) {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "SHIFT_MULTIPLE_MATCHED_RESOLVED\n";
                                std::cout << oss.str();
                            }
                        }
                        return resolved;
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "SHIFT_MULTIPLE_RESOLVED_QUALITY_FAILED - rejecting barcode\n";
                                std::cout << oss.str();
                            }
                        }
                        return std::nullopt;
                    }
                }
            }
        }
        // === no match ===
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "NO_MATCH_FOUND\n";
                std::cout << oss.str();
            }
        }
        return std::nullopt;
    }

};

/**
 * @brief a sequence with its associated metadata and elements
 */
class SigString {
    SigElement sig_elements;
    std::string sequence_id;
    int sequence_length;
    std::string read_type;

public:
    // Constructor
    SigString(std::string id = "", int length = 0, std::string type = "undefined")
        : sequence_id(std::move(id)),
          sequence_length(length),
          read_type(std::move(type)) {}

    // Metadata accessors
    const std::string& id() const { return sequence_id; }
    int length() const { return sequence_length; }
    const std::string& type() const { return read_type; }

    // Container access methods
    const SigElement& elements() const { return sig_elements; }
    SigElement& elements() { return sig_elements; }

    // Index accessors
    auto& by_id() { return sig_elements.get<sig_id_tag>(); }
    auto& by_global() { return sig_elements.get<sig_global_tag>(); }
    auto& by_edit_distance() { return sig_elements.get<sig_ed_tag>(); }
    auto& by_order() { return sig_elements.get<sig_order_tag>(); }
    auto& by_direction() { return sig_elements.get<sig_dir_tag>(); }
    auto& by_pass() { return sig_elements.get<sig_pass_tag>(); }
    auto& by_read() { return sig_elements.get<sig_read_tag>(); }

    // Const versions of index accessors
    const auto& by_id() const { return sig_elements.get<sig_id_tag>(); }
    const auto& by_global() const { return sig_elements.get<sig_global_tag>(); }
    const auto& by_edit_distance() const { return sig_elements.get<sig_ed_tag>(); }
    const auto& by_order() const { return sig_elements.get<sig_order_tag>(); }
    const auto& by_direction() const { return sig_elements.get<sig_dir_tag>(); }
    const auto& by_pass() const { return sig_elements.get<sig_pass_tag>(); }
    const auto& by_read() const { return sig_elements.get<sig_read_tag>(); }

    // Element manipulation
    void add_element(const seq_element& element) {
        sig_elements.insert(element);
    }
    
    void add_element(seq_element&& element) {
        sig_elements.insert(std::move(element));
    }

    // Metadata setters
    void set_id(const std::string& id) { sequence_id = id; }
    void set_length(int length) { sequence_length = length; }
    void set_type(const std::string& type) { read_type = type; }

    // Container operations
    size_t size() const { return sig_elements.size(); }
    bool empty() const { return sig_elements.empty(); }
    void clear() { 
        sig_elements.clear(); 
    }

    // Iterators
    auto begin() { return sig_elements.begin(); }
    auto end() { return sig_elements.end(); }
    auto begin() const { return sig_elements.begin(); }
    auto end() const { return sig_elements.end(); }

    template<typename mod> bool edit_elem(const std::string &class_id, mod m) {
      auto &idx = sig_elements.get<sig_id_tag>();
      auto it = idx.find(class_id);
      if (it == idx.end()){
        return false;
      }
      idx.modify(it, m);
      return true;
    }

private:
    //switched from map to multimap to be able to reference multiple positions    
    bool generate_variable_elements(const ReadElement* layout_elem,  const std::multimap<std::string, const seq_element*>& static_refs,
                                  const std::string& read_seq, SigElement& sig_elements, bool verbose) {
        auto& id_index = sig_elements.get<sig_id_tag>();
        auto it = id_index.find(layout_elem->class_id);
        if (it == id_index.end() || !layout_elem->ref_pos){
            return false;
        }
        const auto& ref_pos = *layout_elem->ref_pos;
        std::pair<int, int> var_positions = {-1, -1};

        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Processing variable element: " << layout_elem->class_id << std::endl;
                std::cout << oss.str();
            }
        }
        return map_positions(ref_pos, static_refs, read_seq, it, id_index, var_positions, verbose);
    }

    bool validate_var_positions(const std::pair<int, int>& positions, size_t read_length) {
        return positions.first > 0 && 
               positions.second > 0 && 
               positions.first <= static_cast<int>(read_length) &&
               positions.second <= static_cast<int>(read_length) &&
               positions.second > positions.first;
    }

    // master function to map variable element positions--super bulky and unweldy, but works. needs to be broken up into constituents
    bool map_positions(const ReferencePositions& ref_pos, const std::multimap<std::string, const seq_element*>& static_refs,
                   const std::string& read_seq,  SigElement::index<sig_id_tag>::type::iterator var_it,
                   SigElement::index<sig_id_tag>::type& id_index,  std::pair<int, int>& var_positions,
                   bool verbose) {
    // Build an ordered list of static elements with the same direction as the variable element.
    std::vector<const seq_element*> ordered_elements;
    for (const auto& [key, elem] : static_refs) {
        if (elem->direction == var_it->direction) {
            ordered_elements.push_back(elem);
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Added static reference " << elem->class_id 
                              << " (order: " << elem->order 
                              << ", pos: " << elem->position.first << "-" << elem->position.second 
                              << ") to ordered list." << std::endl;
                    std::cout << oss.str();
                }
            }
        }
    }
    std::sort(ordered_elements.begin(), ordered_elements.end(),
    [](const seq_element* a, const seq_element* b) {
         if (a->position.first != b->position.first)
              return a->position.first < b->position.first;
         return a->order < b->order;
    });
    
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Ordered static references for variable element " 
                      << var_it->class_id << ": ";
            std::cout << oss.str();
            for (const auto& elem : ordered_elements) {
                std::ostringstream oss_elem;
                oss_elem << elem->class_id << "(" << elem->order << ") ";
                std::cout << oss_elem.str();
            }
            std::cout << std::endl;
        }
    }
    
    int start_pos = -1;
    // ----- START POSITION LOOKUP -----
    { // Primary start lookup
        auto range = static_refs.equal_range(ref_pos.primary_start.ref_id);
        const seq_element* primary_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            // Use candidate if it's present in the ordered list.
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                primary_candidate = it->second;
                break;
            }
        }
        if (primary_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.primary_start.ref_id 
                              << " as the primary start reference candidate.\n"
                              << "Primary start details: is_start=" << ref_pos.primary_start.is_start 
                              << ", offset=" << ref_pos.primary_start.offset 
                              << ", static pos=(" << primary_candidate->position.first << ","
                              << primary_candidate->position.second << ")" << std::endl;
                    std::cout << oss.str();
                }
            }
            auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
                                              [&](const seq_element* elem) {
                                                  return elem->class_id == ref_pos.primary_start.ref_id;
                                              });
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ", global_class: " << (*prev)->global_class << ")." << std::endl;
                        std::cout << oss.str();
                    }
                }
                if ((*prev)->global_class != "poly_tail" &&
                    (*prev)->position.second > primary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Primary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_candidate->position.first << ")." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary start candidate is the first element in the ordered list." << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
            if (is_ordered) {
                start_pos = ref_pos.primary_start.is_start ?
                    primary_candidate->position.first + ref_pos.primary_start.offset :
                    primary_candidate->position.second + ref_pos.primary_start.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using primary start mapping: calculated start_pos = " << start_pos << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Primary start reference " << ref_pos.primary_start.ref_id 
                              << " not found in static_refs." << std::endl;
                    std::cout << oss.str();
                }
            }
        }
    }
    
    // If primary failed, try secondary start.
    if (start_pos <= 0) {
        auto range = static_refs.equal_range(ref_pos.secondary_start.ref_id);
        const seq_element* secondary_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                secondary_candidate = it->second;
                break;
            }
        }
        auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
        [&](const seq_element* elem) {
            return elem->class_id == ref_pos.secondary_start.ref_id;
        });
        if (secondary_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.secondary_start.ref_id 
                              << " as the secondary start reference candidate."
                              << "Secondary start details: is_start=" << ref_pos.secondary_start.is_start 
                              << ", offset=" << ref_pos.secondary_start.offset 
                              << ", static pos=(" << secondary_candidate->position.first << ","
                              << secondary_candidate->position.second << ")" << std::endl;
                    std::cout << oss.str();
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                        std::cout << oss.str();
                    }
                }
                if ((*prev)->position.second > secondary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Secondary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << secondary_candidate->position.first << ")." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary start candidate is the first element in the ordered list." << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
            if (is_ordered) {
                start_pos = ref_pos.secondary_start.is_start ?
                    secondary_candidate->position.first + ref_pos.secondary_start.offset :
                    secondary_candidate->position.second + ref_pos.secondary_start.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using secondary start mapping: calculated start_pos = " << start_pos << std::endl;
                        std::cout << oss.str();
                    }
                }
            } 
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary start reference " << ref_pos.secondary_start.ref_id 
                              << " not found in static_refs." << std::endl;
                    std::cout << oss.str();
                }
            }
        if (!ref_pos.secondary_start.add_flags.empty() && ref_pos.secondary_start.add_flags == "left_truncated") {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary start flagged as left_truncated; attempting fallback." << std::endl;
                    std::cout << oss.str();
                }
            }
                // Look for the nearest static element to the right.
            auto nearest_static = std::find_if(ref_order_pos, ordered_elements.end(),
                                                [](const seq_element* elem) { return elem->type == "static"; });
                if (nearest_static != ordered_elements.end()) {
                    auto fallback_range = static_refs.equal_range((*nearest_static)->class_id);
                    const seq_element* fallback_candidate = nullptr;
                    for (auto it = fallback_range.first; it != fallback_range.second; ++it) {
                        if (it->second == *nearest_static) {
                            fallback_candidate = it->second;
                            break;
                        }
                    }
                    if (fallback_candidate != nullptr) {
                        start_pos = fallback_candidate->position.first + ref_pos.secondary_start.offset;
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Using left terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with start_pos = " << start_pos << std::endl;
                                std::cout << oss.str();
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Fallback static reference for secondary start not found." << std::endl;
                                std::cout << oss.str();
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "No suitable fallback static element found for secondary start." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            }
        }
    }
    
    // ----- STOP POSITION LOOKUP -----
    int stop_pos = -1;
    { // Primary stop lookup
        auto range = static_refs.equal_range(ref_pos.primary_stop.ref_id);
        const seq_element* primary_stop_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                primary_stop_candidate = it->second;
                break;
            }
        }
        if (primary_stop_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.primary_stop.ref_id 
                              << " as the primary stop reference candidate.\n"
                              << "Primary stop details: is_start=" << ref_pos.primary_stop.is_start 
                              << ", offset=" << ref_pos.primary_stop.offset 
                              << ", static pos=(" << primary_stop_candidate->position.first << ","
                              << primary_stop_candidate->position.second << ")" << std::endl;
                    std::cout << oss.str();
                }
            }
            auto ref_order_pos = std::find_if(ordered_elements.begin(), ordered_elements.end(),
                                              [&](const seq_element* elem) {
                                                  return elem->class_id == ref_pos.primary_stop.ref_id;
                                              });
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary stop: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                        std::cout << oss.str();
                    }
                }
                if ((*prev)->position.second > primary_stop_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Primary stop out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_stop_candidate->position.first << ")." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Primary stop candidate is the first element in the ordered list." << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
            if (is_ordered) {
                stop_pos = ref_pos.primary_stop.is_start ?
                    primary_stop_candidate->position.first + ref_pos.primary_stop.offset :
                    primary_stop_candidate->position.second + ref_pos.primary_stop.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using primary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Primary stop reference " << ref_pos.primary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                    std::cout << oss.str();
                }
            }
        }
    }
    
    // If primary stop mapping failed, try secondary.
    if (stop_pos <= 0) {
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Primary stop mapping failed; attempting secondary stop mapping.\n"
                          << "Secondary stop reference ID: " << ref_pos.secondary_stop.ref_id
                          << "\nSecondary flags: " << ref_pos.secondary_stop.add_flags << std::endl;
                std::cout << oss.str();
            }
        }
        auto range = static_refs.equal_range(ref_pos.secondary_stop.ref_id);
        const seq_element* secondary_stop_candidate = nullptr;
        for (auto it = range.first; it != range.second; ++it) {
            if (std::find(ordered_elements.begin(), ordered_elements.end(), it->second) != ordered_elements.end()) {
                secondary_stop_candidate = it->second;
                break;
            }
        }
        //moved this out of the verbosity loop
        auto ref_order_pos = std::find_if(
            ordered_elements.begin(), ordered_elements.end(),[&](const seq_element* elem) {
               return elem->class_id == ref_pos.secondary_stop.ref_id;
                }
        );
        if (secondary_stop_candidate != nullptr) {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found " << ref_pos.secondary_stop.ref_id 
                              << " as the secondary stop reference candidate."
                              << "Secondary stop details: is_start=" << ref_pos.secondary_stop.is_start 
                              << ", offset=" << ref_pos.secondary_stop.offset 
                              << ", static pos=(" << secondary_stop_candidate->position.first << ","
                              << secondary_stop_candidate->position.second << ")" << std::endl;
                    std::cout << oss.str();
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop: Previous element in ordered list is " 
                                << (*prev)->class_id << " (pos: " 
                                << (*prev)->position.first << "-" << (*prev)->position.second 
                                << ")." << std::endl;
                        std::cout << oss.str();
                    }
                }
                if ((*prev)->position.second > secondary_stop_candidate->position.first)
                    is_ordered = false;
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop candidate is the first element in the ordered list." << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
            if (is_ordered) {
                stop_pos = ref_pos.secondary_stop.is_start ?
                    secondary_stop_candidate->position.first + ref_pos.secondary_stop.offset :
                    secondary_stop_candidate->position.second + ref_pos.secondary_stop.offset;
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Using secondary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Secondary stop reference " << ref_pos.secondary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                    std::cout << oss.str();
                }
            }
            if (!ref_pos.secondary_stop.add_flags.empty() && ref_pos.secondary_stop.add_flags == "right_truncated") {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Secondary stop flagged as right_truncated; attempting fallback." << std::endl;
                        std::cout << oss.str();
                    }
                }
                // Look for the nearest static element to the left.
                auto nearest_static = std::find_if(std::make_reverse_iterator(ref_order_pos),
                                                    ordered_elements.rend(),
                                                    [](const seq_element* elem) { return elem->type == "static"; });
                if (nearest_static != ordered_elements.rend()) {
                    auto fallback_range = static_refs.equal_range((*nearest_static)->class_id);
                    const seq_element* fallback_candidate = nullptr;
                    for (auto it = fallback_range.first; it != fallback_range.second; ++it) {
                        if (it->second == *nearest_static) {
                            fallback_candidate = it->second;
                            break;
                        }
                    }
                    if (fallback_candidate != nullptr) {
                        stop_pos = fallback_candidate->position.second + ref_pos.secondary_stop.offset;
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Using right terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with stop_pos = " << stop_pos << std::endl;
                                std::cout << oss.str();
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Fallback static reference for secondary stop not found." << std::endl;
                                std::cout << oss.str();
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "No suitable fallback static element found for secondary stop." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            }
        }
    }    
    
    var_positions = {start_pos, stop_pos};
    
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Final calculated positions: " << start_pos << ":" << stop_pos 
                      << " for variable element " << var_it->class_id << std::endl;
            std::cout << oss.str();
        }
    }
    
    // Validate positions relative to read length.
    if (validate_var_positions(var_positions, read_seq.length())) {
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Positions validated for variable element " << var_it->class_id << std::endl;
                std::cout << oss.str();
            }
        }
        std::string var_seq = read_seq.substr(var_positions.first - 1, 
                                              var_positions.second - var_positions.first + 1);
        id_index.modify(var_it, [&](seq_element& elem) {
            elem.position = var_positions;
            elem.seq = var_seq;
            elem.element_pass = true;
        });
        return true;
    }
    
    id_index.modify(var_it, [](seq_element& elem) {
        elem.element_pass = false;
    });
    if (verbose) {
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "Position validation failed for variable element " << var_it->class_id << std::endl;
            std::cout << oss.str();
        }
    }
    return false;
}

    // check if the read is long enough to contain the expected length of static adapters
    bool read_adapter_sum_comp(const std::string& read_seq, int total_expected_length) const {
        return (read_seq.length() >= 100 && read_seq.length() >= static_cast<size_t>(total_expected_length));
    }
    
    // calculate the total expected length of static adapters
    int calc_total_static_len(const ReadLayout& layout) const {
        int total = 0;
        for (const auto& elem : sig_elements) {
            auto layout_elem = layout.by_id().find(elem.class_id);
            if (layout_elem != layout.by_id().end() && layout_elem->expected_length) {
                total += *layout_elem->expected_length;
            }
        }
        return total;
    }

    // filter reads that are shorter than adapter length summatively
    bool filter_short_reads(const std::vector<std::reference_wrapper<const seq_element>>& elems, const ReadLayout& layout) const {
        size_t read_len = 0;
        // find the single "read" element
        for (auto& e_ref : elems) {
            auto const& e = e_ref.get();
            if (e.global_class == "read" && e.seq) {
                read_len = e.seq->size();
                break;
            }
        }
        // if no read at all -> drop
        if (read_len == 0){
            return false;
        }
        // sum up all the expected static‐adapter lengths from the layout
        int adapters = calc_total_static_len(layout);
        // if the read is long enough, say yes; if the read is too short, say no
        return read_len >= static_cast<size_t>(adapters);
    }

    // validate the positions of the elements
    bool validate_sig_positions(seq_element& e, const std::string& direction, bool verbose){
        if (e.position.first <= 0 || e.position.second <= 0 || e.position.second <= e.position.first) {
            auto& idx = sig_elements.get<sig_id_tag>();
            idx.modify(idx.find(e.class_id), [](seq_element& x){ x.element_pass = false; });
            if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "  ["<< direction << "] invalid positions on ""<< e.class_id << ""\n";
                        std::cout << oss.str();
                    }
                }
                return false;
            }
        return true;
    }

    // filter out overlapping elements
    void filter_overlaps(const std::vector<std::reference_wrapper<seq_element>>& elems, bool verbose) {
        for (size_t i = 0; i < elems.size(); ++i) {
            auto &e1 = elems[i].get();
            if (e1.global_class=="start"||e1.global_class=="stop"
             || e1.type!="variable") continue;
            for (size_t j = i+1; j < elems.size(); ++j) {
                auto &e2 = elems[j].get();
                if (e2.global_class=="start"||e2.global_class=="stop"
                 || e2.type!="variable") continue;
                if (e1.position.second >= e2.position.first + 3) {
                    edit_elem(e1.class_id, [](seq_element &x){ x.element_pass = false; });
                    edit_elem(e2.class_id, [](seq_element &x){ x.element_pass = false; });
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "  Overlap detected between elements: "
                                  << e1.class_id << " and " << e2.class_id
                                  << " at positions (" << e1.position.first << "-" 
                                  << e1.position.second << ") and ("
                                  << e2.position.first << "-" 
                                  << e2.position.second << ")." << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
            }
        }
    }

    // Group elements by direction
    auto group_directionally(){
        std::map<std::string, std::vector<std::reference_wrapper<const seq_element>>> direction_elements;
        for (auto& elem : sig_elements)
        direction_elements[elem.direction].push_back(std::ref(elem));
        return direction_elements;
    }

    //  update_bc_counts with whitelist selection and global filtering
    void update_bc_counts(SigString &sig, const ReadLayout &layout, bool verbose, bool skip_global_writes = true) {
        auto type = sig.read_type; // "forward", "reverse", "concatenate", or "filtered"
        bool has_global_only_barcodes = false;
    
        for (auto const &elem : sig.elements()) {
            if (elem.global_class != "barcode" || !elem.seq.has_value()) {
                continue;
            }
            
            // Get current sequence and original sequence
            std::string bc_seq = elem.seq.value();
            std::string original_seq = elem.original_seq.value_or("");
            int64_seq bc(bc_seq);
            
            // Look up the right whitelist
            auto key = seq_utils::remove_rc(elem.class_id);
            auto &wl = layout.wl_map.maps.at(key).get();
            
            // Determine if this was corrected or raw
            bool is_corrected = false;
            if (!original_seq.empty() && original_seq != "") {
                std::string original_for_comparison = original_seq;
                std::string original_revcomp = seq_utils::revcomp(original_seq);
                if (bc_seq != original_for_comparison && bc_seq != original_revcomp) {
                    is_corrected = true;
                }
            }

            // === select the right wl ===
            bool found_in_true = wl.true_bcs.check_wl_for(bc);
            bool found_in_global = wl.global_bcs.check_wl_for(bc);
            // Simple selection: true or global
            auto* target_whitelist = &wl.filter_bcs;
            std::string whitelist_used = "filtered";
            if (found_in_true) {
                target_whitelist = &wl.true_bcs;
                whitelist_used = "true";
            } else if (found_in_global) {
                target_whitelist = &wl.global_bcs;
                whitelist_used = "global";
                has_global_only_barcodes = true;
            }

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Updating " << bc.bits_to_sequence() << " in " << whitelist_used 
                        << " whitelist (" << (is_corrected ? "CORRECTED" : "RAW") << ")" << std::endl;
                    std::cout << oss.str();
                }
            }
            // === update counts in selected wl ===
            if (type == "filtered" || type != elem.direction && type != "concatenate") {
                target_whitelist->update_bc_count(bc, barcode_counts::filtered);
            }

            if (type == "forward" && elem.direction == "forward") {
                target_whitelist->update_bc_count(bc, barcode_counts::total);
                target_whitelist->update_bc_count(bc, barcode_counts::forw);
                if (is_corrected) {
                    target_whitelist->update_bc_count(bc, barcode_counts::corrected);
                } else {
                    target_whitelist->update_bc_count(bc, barcode_counts::raw);
                }
            }

            if (type == "reverse" && elem.direction == "reverse") {
                target_whitelist->update_bc_count(bc, barcode_counts::total);
                target_whitelist->update_bc_count(bc, barcode_counts::rev);
                if (is_corrected) {
                    target_whitelist->update_bc_count(bc, barcode_counts::corrected);
                } else {
                    target_whitelist->update_bc_count(bc, barcode_counts::raw);
                }
            }

            if (type == "concatenate") {
                if (elem.direction == "forward") {
                    target_whitelist->update_bc_count(bc, barcode_counts::total);
                    target_whitelist->update_bc_count(bc, barcode_counts::forw_concat);
                    if (is_corrected) {
                        target_whitelist->update_bc_count(bc, barcode_counts::corrected);
                    } else {
                        target_whitelist->update_bc_count(bc, barcode_counts::raw);
                    }
                }

                if (elem.direction == "reverse") {
                    target_whitelist->update_bc_count(bc, barcode_counts::total);
                    target_whitelist->update_bc_count(bc, barcode_counts::rev_concat);
                    if (is_corrected) {
                        target_whitelist->update_bc_count(bc, barcode_counts::corrected);
                    } else {
                        target_whitelist->update_bc_count(bc, barcode_counts::raw);
                    }
                }
            }
        }
        // === handle global filters ===
        if (skip_global_writes && has_global_only_barcodes && sig.read_type != "filtered") {
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Read " << sig.sequence_id << " has global-only barcodes."
                        << "Setting read_type to 'skipped' due to skip_global_writes = true" << std::endl;
                    std::cout << oss.str();
                }
            }
            sig.set_type("skipped");
        }
    }
public:
    //sigalign_static
   void sigalign_static(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
    aligner_tools aligner;
    auto& type_index = layout.by_type();
    auto static_range = type_index.equal_range("static");

    if(verbose){
        #pragma omp critical
        {
            std::ostringstream oss;
            oss << "\n=== Starting static alignment for " << sequence_id << " ===" << std::endl;
            std::cout << oss.str();
        }
    }
    int max_distance = -1;
    // Make a mutable copy of the read so we can mask aligned regions.
    std::string mutable_seq = read.seq;
    size_t read_length = read.seq.length();

    for (auto it = static_range.first; it != static_range.second; ++it) {
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Processing static element: " << it->class_id << std::endl;
                std::cout << oss.str();
                if(it->aligned_positions){
                    std::ostringstream aligned;
                    aligned << "  Aligned positions: " 
                              << it->aligned_positions->start_stats.first << ", "
                              << it->aligned_positions->start_stats.second << std::endl;
                    std::cout << aligned.str();
                } else {
                    std::ostringstream fail;
                    fail << "  No aligned positions available." << std::endl;
                    std::cout << fail.str();
                }
                if(it->misaligned_positions){
                    std::ostringstream misaligned;
                    misaligned << "  Misaligned positions: " 
                              << it->misaligned_positions->start_stats.first << ", "
                              << it->misaligned_positions->start_stats.second << std::endl;
                    std::cout << misaligned.str();
                } else {
                    std::ostringstream fail_again;
                    fail_again << "  No misaligned positions available." << std::endl;
                    std::cout << fail_again.str();
                }
            }
        }
        // For 'start' or 'stop' types, add an element without alignment.
        if (it->global_class == "start" || it->global_class == "stop") {
            add_element(seq_element(
                it->class_id,
                it->global_class,
                std::nullopt,
                (it->global_class == "start")
                    ? std::make_pair(1, 1)
                    : std::make_pair(static_cast<int>(read_length), static_cast<int>(read_length)),
                "static",
                it->order,
                "",
                it->direction,
                std::nullopt,
                std::nullopt
            ));
            continue;
        }

        // Update max_distance from misalignment_threshold if available. If not, set to ~20% of the adapter seq length.
        // this will be the case for R1/R2 reads
        if (it->misalignment_threshold) {
            max_distance = std::get<0>(*it->misalignment_threshold);
        } else {
            max_distance = static_cast<int>(it->seq.length() * 0.2);
        }

        // poly-tail solution--importantly, only will look for poly-tails if you tell it to
        if (it->global_class == "poly_tail") {
            auto result = aligner.find_poly_tails(it->seq, read.seq, 14);
            if (result.success) {
                for (const auto& positions : result.positions) {
                    add_element(seq_element(
                        it->class_id,
                        it->global_class,
                        result.edit_distance,
                        positions,
                        "static",
                        it->order,
                        "",
                        it->direction,
                        true,
                        mutable_seq.substr(positions.first - 1,
                                             positions.second - positions.first + 1)
                    ));
                    // Mask the found region so that it is not re-aligned.
                    std::fill(mutable_seq.begin() + positions.first - 1,
                              mutable_seq.begin() + positions.second, 'X');
                }
            }
            continue;
        }
        // For other static elements (non-poly_tail), use the alignment statistics
        // to define the expected region for the adapter
        int expected_start, expected_end;
        // Check if the element has aligned and misaligned positions
        if (it->aligned_positions && it->misaligned_positions) {
            const auto& [start_mean, start_var] = it->aligned_positions->start_stats;
            const auto& [mis_start, start_mvar] = it->misaligned_positions->start_stats;
            // Compute expected region boundaries based on the aligned stats
            expected_start = (start_mean < 50.0) ? 1 
            : static_cast<int>(std::max(0.0, ((start_mean - (start_var * 1.25)) / 100.0) * read_length));

            expected_end = (start_mean < 50.0)
                ? static_cast<int>(std::min(((start_mean + (start_var * 1.25)) / 100.0) * read_length, static_cast<double>(read_length)))
                : static_cast<int>(read_length);

        } else {
            // Default to using the entire read length if no statistics are available.
            expected_start = 1;
            expected_end = static_cast<int>(read_length);
        }
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Expected region for " << it->class_id << ": " 
                          << expected_start << " - " << expected_end << std::endl;
                std::cout << oss.str();
            }
        }

        // Use the entire mutable sequence as target.
        auto result = aligner.align_static_elements(it->seq,
                                                    mutable_seq,
                                                    verbose,
                                                    max_distance,
                                                    it->masked_seq,
                                                    verbose,
                                                    expected_start,
                                                    expected_end);
            if (result.success) {
                // Positions returned are relative to the full read.
                for (const auto& pos : result.positions) {
                    int adj_start = pos.first;
                    int adj_end = pos.second;
                    if (adj_start < 1) {
                        adj_start = 1;
                    }
                    if (adj_end > static_cast<int>(read_length)) {
                        adj_end = read_length;
                    }
                    std::pair<int, int> adjusted_positions = {adj_start, adj_end};

                    add_element(seq_element(
                        it->class_id,
                        it->global_class,
                        result.edit_distance,
                        adjusted_positions,
                        "static",
                        it->order,
                        "",
                        it->direction,
                        true,
                        mutable_seq.substr(adj_start - 1,
                                             adj_end - adj_start + 1)
                    ));

                    // Mask the aligned region so that it is not re-aligned.
                    std::fill(mutable_seq.begin() + adj_start - 1,
                              mutable_seq.begin() + adj_end, 'X');
                }
                // Primary region yielded a valid match; skip further processing for this element.
                continue;
            }
        }
}
    //this version iterates across reverse elements in reverse order, which fixes issues with rc umi/rc barcode dependencies
    void sigalign_variable(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "\n=== Starting variable alignment for " << sequence_id << " ===" << std::endl;
                std::cout << oss.str();
            }
        }

        // Split variables into read and non-read
        std::vector<const ReadElement*> non_read_vars;
        std::vector<const ReadElement*> read_vars;
        
        auto& type_index = layout.by_type();
        auto variable_range = type_index.equal_range("variable");
        
        for (auto it = variable_range.first; it != variable_range.second; ++it) {
            // Add a placeholder element into our sig_elements container.
            add_element(seq_element(
                it->class_id,
                it->global_class,
                std::nullopt,
                {-1, -1},
                "variable",
                it->order,
                "",
                it->direction,
                std::nullopt,
                std::nullopt
            ));
        
            if (it->global_class == "read") {
                read_vars.push_back(&(*it));
            } else {
                non_read_vars.push_back(&(*it));
            }
        }
        
        // Partition non-read variables by direction.
        std::vector<const ReadElement*> non_read_forward;
        std::vector<const ReadElement*> non_read_reverse;
        for (const auto* var : non_read_vars) {
            if (var->direction == "reverse") {
                non_read_reverse.push_back(var);
            } else {
                non_read_forward.push_back(var);
            }
        }
        
        int positioned_count = 0;
        
        // Process forward non-read variables in natural order.
        std::multimap<std::string, const seq_element*> static_refs;
        for (const auto& elem : sig_elements) {
            if (elem.type == "static") {
                static_refs.insert({elem.class_id, &elem});
                if (verbose) {
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "Found static reference: " << elem.class_id 
                                  << " at position " << elem.position.first 
                                  << ":" << elem.position.second 
                                  << " with edit distance " << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "none")
                                  << " and sequence " << (elem.seq ? elem.seq.value() : "none") 
                                  << std::endl;
                        std::cout << oss.str();
                    }
                }
            }
        }
        
        // Process forward non-read variables.
        for (const auto* var : non_read_forward) {
            if (generate_variable_elements(var, static_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
                // Add the newly generated variable element(s) to the static_refs multimap so that secondary positions can be calculated.
                auto found = sig_elements.get<sig_id_tag>().find(var->class_id);
                if (found != sig_elements.get<sig_id_tag>().end()) {
                    static_refs.insert({var->class_id, &(*found)});
                }
            }
        }
        
        // Process reverse non-read variables in reverse order.
        for (auto it = non_read_reverse.rbegin(); it != non_read_reverse.rend(); ++it) {
            const auto* var = *it;
            if (generate_variable_elements(var, static_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
                auto found = sig_elements.get<sig_id_tag>().find(var->class_id);
                if (found != sig_elements.get<sig_id_tag>().end()) {
                    static_refs.insert({var->class_id, &(*found)});
                }
            }
        }
        
        // Build a container for total references.
        std::multimap<std::string, const seq_element*> total_refs;
        for (const auto& elem : sig_elements) {
            total_refs.insert({elem.class_id, &elem});
            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "Found reference: " << elem.class_id 
                              << " at position " << elem.position.first 
                              << ":" << elem.position.second << std::endl;
                    std::cout << oss.str();
                }
            }
        }
        // Process read variables last.
        for (const auto* var : read_vars) {
            if (generate_variable_elements(var, total_refs, read.seq, sig_elements, verbose)) {
                positioned_count++;
            }
        }
        if (verbose) {
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Successfully positioned " << positioned_count << " variable elements\n"
                 << "Final sigstring elements: " << sig_elements.size() << std::endl;
                std::cout << oss.str();
            }
        }
    }
    // contains per-unit barcode correction
    void sigalign_filter(const read_streaming::sequence &read, const ReadLayout& layout, 
                         int gen_mut, int gen_shift, bool verbose) {
        // make a resever for the qual that doesn't interfere w/ASCII characters
        constexpr char qual_mask = '\x7F';
        auto direction_elements = group_directionally();
        // For recording validity and count per direction.
        std::map<std::string, bool> direction_valid;
        std::map<std::string, int> pass_counts, static_counts;
        std::unordered_map<std::string, std::unordered_set<int64_seq>> seen_bcs;
        bool skip_global_writes = false;
        aligner_tools aligner;

        // Process each direction separately.
        for (auto& [direction, elements] : direction_elements) {
            int barcode_count = 0;
            bool multiple_barcodes = false;
            std::string masked_read = read.seq;
            std::string masked_qual = read.is_fastq ? read.qual : "";
            //step one--mask elements within the read from overlapping adapters/poly-tails
            for (auto const& ref : elements) {
                const auto& e = ref.get();
                if (e.global_class == "read" || e.position.first <= 0 || e.position.second <= e.position.first){
                    continue;
                }
                // for this direction, if there's more than one barcode, we can set the multiple barcodes flag
                // to true
                if (e.global_class == "barcode"){
                    barcode_count++;
                    if (barcode_count > 1) {
                        multiple_barcodes = true;
                    }
                }
                size_t s = e.position.first - 1;
                size_t len = e.position.second - e.position.first + 1;
                //this populates the masked read with 'N's
                //if the read is shorter than the expected length, skip
                if (s + len <= masked_read.size()) {
                    std::fill(masked_read.begin() + s, masked_read.begin() + s + len, 'N');
                    if (read.is_fastq && masked_qual.size() == masked_read.size()) {
                        std::fill(masked_qual.begin() + s, masked_qual.begin() + s + len, qual_mask);
                    }
                }
            }
            
            //step two--trim the read to remove 'N's
            for (auto const& ref : elements) {
                const auto& e = ref.get();
                if (e.global_class != "read" || !e.seq){
                    continue;
                }

                size_t start = e.position.first - 1;
                size_t length = e.position.second - e.position.first + 1;
                if (start + length >= masked_read.size() || start < 1 || length <= 1 || start + length <= 1) {
                    if(verbose){
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Skipping a read element due to out-of-bounds parameters.\n"
                                    << "Start: " << start << ", Length: " << length 
                                    << ", Read length: " << masked_read.size() << "\n"
                                    << "Masked read: " << masked_read << "\n"
                                    << "Direction: " << direction << std::endl;
                            std::cout << oss.str();
                        }
                    }
                    break;  // out of bounds—skip
                }

                // extract, then strip all 'N's
                std::string window = masked_read.substr(start, length);
                std::string window_qual = read.is_fastq ? masked_qual.substr(start, length) : "";

                std::string cleaned_seq;
                std::string cleaned_qual;
                // reserve space for cleaned sequences
                cleaned_seq.reserve(window.size());
                // if the read is fastq, reserve space for qualscores
                if (read.is_fastq){
                    cleaned_qual.reserve(window.size());
                }

                for (size_t i = 0; i < window.size(); ++i) {
                    if (window[i] != 'N') {
                        cleaned_seq.push_back(window[i]);
                        if (read.is_fastq) {
                            // use the corresponding quality character
                            cleaned_qual.push_back(window_qual[i]);
                        }
                    }
                }
                // small lambda for dealing with editing elements
                edit_elem(e.class_id, [cleaned_seq, cleaned_qual, &read, verbose](seq_element &el) {
                    *el.seq = cleaned_seq;
                    if (read.is_fastq) {
                        el.qual = cleaned_qual;
                        if(verbose){
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Masked qual: " << cleaned_qual << "\n";
                                std::cout << oss.str();
                            }
                        }
                    } else {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Masked qual: " << el.qual.value() << "\n";
                            std::cout << oss.str();
                        }
                    }
                });

                // verbose logging
                if (verbose) {
                    #pragma omp critical
                    {
                    std::ostringstream oss;
                    oss << "Trimmed read for " << e.class_id 
                              << " -> " << cleaned_seq << "\n"
                              << "Masked read: " << masked_read << std::endl;
                    std::cout << oss.str();

                    }
                    if(read.is_fastq) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Masked qual: " << cleaned_qual << std::endl;
                            std::cout << oss.str();
                        }
                    }
                }
                break;  // only one read element per direction
            }

            bool valid_direction = true;
            int count = 0;
            
            //per-direction check and length filter
            if (!filter_short_reads(elements, layout)) {
                direction_valid[direction] = false;
                if(verbose){
                    #pragma omp critical
                    {
                        std::ostringstream oss;
                        oss << "FILTERED_READ_LENGTH" << std::endl;
                        std::cout << oss.str();
                    }
                }
                continue;
            }

            // Sort by order.
            std::sort(elements.begin(), elements.end(), [](const auto& a, const auto& b) {
                return a.get().order < b.get().order;
            });

            // per element checks in this direction, and filters.
            for (auto &elem_ref : elements) {
                int statics = 0;
                const auto& elem = elem_ref.get();
                // Skip "start" or "stop" elements for filtering/counting.
                if (elem.global_class == "start" || elem.global_class == "stop" || elem.global_class == "poly_tail"){
                    continue;
                }
                // Count static elements in this direction.
                if(elem.type == "static"){
                    ++statics;
                    static_counts[direction] = statics;
                }
               
                // If any variable element has invalid positions, mark the entire direction as invalid and flag the element.
                if (elem.position.first <= 0 || elem.position.second <= 0 || elem.position.second <= elem.position.first) {
                    valid_direction = false;
                    auto& id_index = sig_elements.get<sig_id_tag>();
                    id_index.modify(
                        id_index.find(elem.class_id), [](seq_element &e) { 
                                        e.element_pass = false; 
                                    });
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << 
                            "Filtered element " << 
                            elem.class_id << 
                            " (invalid positions) in direction " << direction << std::endl;
                            std::cout << oss.str();
                        }
                    }
                    // Exit processing for this direction.
                    break;
                }
                //barcode correction here!!!
                if (elem.global_class == "barcode") {
                    if(verbose){
                        #pragma omp critical
                        {
                            std::ostringstream oss;
                            oss << "Barcode correction for " << elem.seq.value() << std::endl;
                            std::cout << oss.str();
                        }
                    }
                    auto& wl = layout.wl_map.maps.at(seq_utils::remove_rc(elem.class_id)).get();
                    //had originally included a check for global barcode class not empty (&& !wl.true_bcs.empty())
                    //but this is not necessary, as the barcode correction function will handle it
                    if(!wl.true_bcs.empty()){
                        auto out = barcode_correction::correct_barcode(elem, layout, verbose, gen_mut, gen_shift);
                        bool matched = out.has_value();
                        if(matched) {
                            int64_seq original_bc;
                            original_bc.sequence_to_bits(elem.seq.value());
                            int64_seq correct_bc = out.value();
                            seen_bcs[seq_utils::remove_rc(elem.class_id)].insert(out.value());
                            //fixing this because corrected is only in the forward direction anyways so need to revcomp
                            std::string final_bc = correct_bc.bits_to_sequence();
                            // Check if the corrected barcode exists in true_bcs or global_bcs
                            bool found_in_true = wl.true_bcs.check_wl_for(correct_bc);
                            bool found_in_global = wl.global_bcs.check_wl_for(correct_bc);
                            //first check if the corrected barcode is the same as an original sequence
                            bool already_keyed = wl.true_bcs.check_wl_for(original_bc) || wl.global_bcs.check_wl_for(original_bc);         
                            if (!already_keyed) {
                                if (found_in_true &! found_in_global) {
                                    auto true_range = wl.true_bcs.equal_range(correct_bc);
                                    if (!wl.true_bcs.check_wl_for(original_bc) && true_range.first != true_range.second) {
                                        const barcode_entry& correct_entry = (*true_range.first);
                                        #pragma omp critical
                                        {
                                            wl.true_bcs.insert_bc_entry(original_bc, correct_entry);
                                        }
                                    }
                                }
                                else if (found_in_global &! found_in_true) {
                                    // Correct barcode found in global_bcs
                                    auto global_range = wl.global_bcs.equal_range(correct_bc);
                                    if (!wl.true_bcs.check_wl_for(original_bc) && global_range.first != global_range.second) {
                                        const barcode_entry& correct_entry = (*global_range.first);
                                        #pragma omp critical
                                        {
                                            wl.global_bcs.insert_bc_entry(original_bc, correct_entry);
                                        }
                                    }
                                }
                            }
                            // update original_seq with the original sequence, set the sequence to the final correccted barcode
                            // and set element pass to be true
                            auto& id_index = sig_elements.get<sig_id_tag>();
                            std::string final_wl;
                            if(found_in_true &! found_in_global){
                                final_wl = "true";
                            } else {
                                final_wl = "global";
                            }

                            id_index.modify(
                                id_index.find(elem.class_id), [&](seq_element &e){
                                    e.original_seq = e.seq;
                                    e.seq = final_bc;
                                    e.element_pass = true;
                                    e.flags = final_wl;
                                }
                            );
                        } else {
                            // else mark the element as failed
                            auto& id_index = sig_elements.get<sig_id_tag>();
                            id_index.modify(
                                id_index.find(elem.class_id),
                                [](seq_element &e) {
                                        e.element_pass = false; 
                                        e.flags = "filter";
                                    }
                                );
                            if(elem.seq.has_value()){
                                int64_seq failed_bc, rc_failed_bc;
                                failed_bc.sequence_to_bits(elem.seq.value());
                                rc_failed_bc.sequence_to_bits(seq_utils::revcomp(elem.seq.value()));
                                
                                #pragma omp critical
                                {
                                    // Check if neither the barcode nor its reverse complement are already keys
                                    if (!wl.filter_bcs.check_wl_for(failed_bc) && !wl.filter_bcs.check_wl_for(rc_failed_bc)) {
                                        if(verbose){
                                            std::ostringstream oss;
                                            oss << "Barcode failed for " << failed_bc.bits_to_sequence() 
                                                << "\nno direct match found for this key, adding as a new barcode\n";
                                            std::cout << oss.str();
                                        }                       
                                        barcode_entry failed_entry;
                                        failed_entry.barcode = failed_bc;
                                        failed_entry.edit_dist = 10;
                                        failed_entry.filtered = true;
                                        failed_entry.flags = "failed";  
                                        // Double-check that it's still not present
                                        auto failed_range = wl.filter_bcs.equal_range(failed_bc);
                                        auto rc_failed_range = wl.filter_bcs.equal_range(rc_failed_bc);
                                        if((failed_range.first == failed_range.second) && (rc_failed_range.first == rc_failed_range.second)) {
                                            // Entry doesn't exist - insert it
                                            wl.filter_bcs.insert_bc_entry(failed_bc, failed_entry);
                                        }
                                    }
                                }
                            }
                            // adding this in to try to manage if there are multiple barcodes, whether we should set this entire direction to false. right now, we *only* return completely correct barcodes (all multiples are also correct) 
                            if(!multiple_barcodes){
                                // so if multiple barcodes is false, then we set the direction value to false and set the pass counts to 0.
                                valid_direction = false;
                                pass_counts[direction] = 0;    
                            }
                            continue;
                        }
                    } /*else {
                        //this branch is for generating a purely standalone whitelist that will occur IFF there is no default whitelist
                        //--default whitelist generation has yet to happen considering just dumping them into the true_bcs. 
                        //will probably pull this out and write a standalone sorting function internally here.
                        if(elem.seq.has_value()){
                            std::string seq = elem.seq.value();
                            int64_seq putative_bc(seq);
                            if(!wl.filter_bcs.check_wl_for(putative_bc)){
                                #pragma omp critical
                                {
                                    wl.filter_bcs.insert_bc_entry(be.barcode, be);
                               }
                            }
                        }
                    }*/
                }
            }

            // If no static elements were found, mark the direction as invalid.
            if(static_counts[direction] == 0){
                valid_direction = false;
                pass_counts[direction] = 0;
            }

            // If the direction failed the per-element check, record zero passing variable elements.
            //removed pass_counts [direction] = 0 because we could blind ourselves to concatenates
            if (!valid_direction) {
                direction_valid[direction] = false;
                //pass_counts[direction] = 0;
                continue;
            }

            // Second pass: Check for overlaps among variable elements (skip start/stop elements).
            for (size_t i = 0; i < elements.size(); i++) {
                const auto& e1 = elements[i].get();
                if (e1.global_class == "start" || e1.global_class == "stop" || e1.type != "variable")
                    continue;
                for (size_t j = i + 1; j < elements.size(); j++) {
                    const auto& e2 = elements[j].get();
                    if (e2.global_class == "start" || e2.global_class == "stop" || e2.type != "variable")
                        continue;
                    // Check for overlap beyond a tolerance (here e1.position.second >= e2.position.first + 3)
                    if (e1.position.second >= e2.position.first + 3) {
                        auto& id_index = sig_elements.get<sig_id_tag>();
                        id_index.modify(id_index.find(e1.class_id),
                                        [](seq_element &e) { e.element_pass = false; });
                        id_index.modify(id_index.find(e2.class_id),
                                        [](seq_element &e) { e.element_pass = false; });
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::ostringstream oss;
                                oss << "Overlap detected: Filtering variable elements " 
                                          << e1.class_id << " and " << e2.class_id 
                                          << " in direction " << direction << std::endl;
                                std::cout << oss.str();
                            }
                        }
                    }
                }
            }
            
            // Count passing variable elements after all modifications.
            for (auto &elem_ref : elements) {
                const auto& elem = elem_ref.get();
                if (elem.type == "variable" && elem.element_pass && 
                    (elem.global_class == "read" || elem.global_class == "barcode"))
                    count++;
            }

            // if we've gotten to this point, then we have a valid direction and a valid pass count
            direction_valid[direction] = true;
            pass_counts[direction] = count;

            if (verbose) {
                #pragma omp critical
                {
                    std::ostringstream oss;
                    oss << "\nDirection " << direction << " valid: " << direction_valid[direction] 
                              << " with " << count << " passing variable elements." << std::endl;
                    std::cout << oss.str();
                }
            }
        }
        
        // Decide the overall read type based on the pass counts.
        int forward_count = (direction_valid.count("forward") && direction_valid.at("forward")) ? pass_counts["forward"] : 0;
        int reverse_count = (direction_valid.count("reverse") && direction_valid.at("reverse")) ? pass_counts["reverse"] : 0;
        
        bool forward_valid = (forward_count > 0 && forward_count >= reverse_count);
        bool reverse_valid = (reverse_count > 0 && reverse_count >= forward_count);

        if (forward_valid && reverse_valid) {
            set_type("concatenate"); 
        } else if (forward_valid) {
            set_type("forward");
        } else if (reverse_valid) {
            set_type("reverse");
        } else {
            set_type("filtered");
        }

        // Update the read barcode counts.
        update_bc_counts(*this, layout, verbose);

        if(verbose){
            #pragma omp critical
            {
                std::ostringstream oss;
                oss << "Final read type: "
                          << read_type
                          << ", forward element(s) count: "
                          << forward_count
                          << ", reverse element(s) count: "
                          << reverse_count
                          << "\n"
                          << std::endl;
                std::cout << oss.str();
            }
        }
    }
     
    //metrics parallelized over reads
    static void sigalign(const std::string& fastq_path, const ReadLayout& layout, const std::string& output_prefix, 
                         std::optional<int> gen_mut, std::optional<int> gen_shift, bool verbose, 
                         int num_threads, size_t chunk_size, size_t max_reads, bool write_debug) {

    // Output file paths
    std::string sig_path, csv_path, metrics_path;
    std::string fastq_output_path = output_prefix + ".fq";
    bool compress_fastq = true;
    // Initialize writers
    std::unique_ptr<std::ofstream> metrics_file_ptr;
    std::unique_ptr<sigstring_writing> sig_writer, csv_writer, metrics_file;
    sigstring_writing fastqa_writer(fastq_output_path, sigstring_writing::format::FASTQA, compress_fastq, /*append=*/false);
    // Initialize the parallel writer
    parallel_writer writer;
    
    if(write_debug){
        sig_path = output_prefix + ".sig";
        csv_path = output_prefix + ".csv";
        metrics_path = output_prefix + ".metrics.tsv";
        sig_writer = std::make_unique<sigstring_writing>(sig_path, sigstring_writing::format::SIGSTRING, /*compress=*/false, /*append=*/false);
        csv_writer = std::make_unique<sigstring_writing>(csv_path, sigstring_writing::format::CSV, /*compress=*/false, /*append=*/false);
        // Write CSV header
        {
            SigString header("", 0);
            (*csv_writer)(std::vector<SigString>{header});
        }
        // Open metrics file
        metrics_file_ptr = std::make_unique<std::ofstream>(metrics_path);
        *metrics_file_ptr << "chunk_id\tseqs_in_chunk\tseqs_passed\tread_time_ms\tprocess_time_ms\tqueue_time_ms\ttotal_time_ms\n";    
    }
    // Initialize file and reader
    file_streaming files(fastq_path);
    read_streaming reader(files);
    
    // Counters
    size_t total_reads = 0;
    size_t total_passed = 0;
    size_t chunk_id = 0;
    
    // Timing totals
    double total_read_time = 0;
    double total_process_time = 0;
    double total_queue_time = 0;
    
    // Process chunks sequentially
    while (true) {
        chunk_id++;

        // Start timing
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Read a chunk of sequences
        std::vector<read_streaming::sequence> chunk;
        chunk.reserve(chunk_size);
        
        size_t local_count = 0;
        while (local_count < chunk_size) {
            if (max_reads > 0 && total_reads >= max_reads) break;
            
            auto seq = reader.next_sequence();
            if (!seq) break;
            
            chunk.push_back(*seq);
            local_count++;
            total_reads++;
        }
        
        // If no more sequences, break
        if (chunk.empty()){
            break;
        }
        // Read time ends here
        auto read_end_time = std::chrono::high_resolution_clock::now();
        
        std::vector<SigString> results;
        results.resize(chunk.size()); 
        // Pre-allocate space for results
        std::atomic<size_t> passed_count{0};
        
        // Process the sequences in parallel
        #pragma omp parallel num_threads(num_threads)
        {
            // Thread-local vector to collect passed SigStrings
            std::vector<SigString> thread_results;
            thread_results.reserve(chunk.size() / num_threads);
            
            #pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < chunk.size(); i++) {
                const auto& read = chunk[i];
                SigString sig(read.id, read.seq.length());
                sig.sigalign_static(read, layout, verbose);
                sig.sigalign_variable(read, layout, verbose);
                sig.sigalign_filter(read, layout, gen_mut.value_or(2), gen_shift.value_or(3), verbose);
                
                if (sig.read_type != "filtered") {
                    //counting everything that's not filtered, but could be skipped
                    thread_results.push_back(std::move(sig));
                }
            }

            // Merge thread results into global results vector
            #pragma omp critical
            {
                for (auto& sig : thread_results) {
                       results[passed_count++] = std::move(sig);
                }
            }
        }
        
        // Resize results to actual number of passed reads
        results.resize(passed_count);
        total_passed += passed_count;
        
        // Process time ends here
        auto process_end_time = std::chrono::high_resolution_clock::now();
        
        // Queue results for writing (non-blocking)
        if (!results.empty()) {
            if(write_debug && sig_writer && csv_writer){
                writer.write_all(*sig_writer, *csv_writer, fastqa_writer, results);
            } else {
                writer.write(fastqa_writer, results);
            }
        }
        
        // Queue time ends here
        auto queue_end_time = std::chrono::high_resolution_clock::now();
        
        // Calculate timing information
        auto read_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(read_end_time - start_time).count();
        auto process_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(process_end_time - read_end_time).count();
        auto queue_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(queue_end_time - process_end_time).count();
        auto total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(queue_end_time - start_time).count();
        
        // Update timing totals
        total_read_time += read_time_ms;
        total_process_time += process_time_ms;
        total_queue_time += queue_time_ms;
        
        // Log metrics
        if(write_debug && metrics_file_ptr){
            *metrics_file_ptr << chunk_id << "\t" << chunk.size() << "\t" << passed_count << "\t"
            << read_time_ms << "\t" << process_time_ms << "\t" << queue_time_ms << "\t" << total_time_ms << "\n";
        }

        
        // Print progress
        std::cout << "[chunk_stats] " << chunk_id << ": processed " << chunk.size() << " reads in " 
                  << total_time_ms / 1000.0 << " seconds ("
                  << read_time_ms / 1000.0 << "s read, "
                  << process_time_ms / 1000.0 << "s process, "
                  << queue_time_ms / 1000.0 << "s queue), "
                  << passed_count << " passed (" << (double)passed_count / chunk.size() * 100.0 << "%)" << std::endl;
        
        // Free memory explicitly
        std::vector<read_streaming::sequence>().swap(chunk);
        //std::vector<read_streaming::sequence>().clear();
        std::vector<SigString>().swap(results);
    }
    
    // Wait for writer to finish all queued writes
    writer.stop();
    // Close metrics file
    if (write_debug && metrics_file_ptr) {
        metrics_file_ptr->close();
    }

   // Calculate total time
    double total_time = total_read_time + total_process_time + total_queue_time;
    
    // Print summary
    std::cout << "\n[sigalign] Performance Summary:" << std::endl;
    std::cout << "────────────────────────────────────────────────────" << std::endl;
    std::cout << "[sigalign] Total runtime: " << total_time / 1000.0 << " seconds" << std::endl;
    std::cout << "[sigalign] Total chunks processed: " << chunk_id << std::endl;
    std::cout << "[sigalign] Total reads processed: " << total_reads << std::endl;
    std::cout << "[sigalign] Reads passing filter: " << total_passed << " (" 
              << (total_reads > 0 ? (total_passed * 100.0 / total_reads) : 0) << "%)" << std::endl;
    std::cout << "[sigalign] Timing breakdown:" << std::endl;
    std::cout << "  - Read time: " << total_read_time / 1000.0 << " seconds (" 
              << (total_read_time * 100.0 / total_time) << "%)" << std::endl;
    std::cout << "  - Process time: " << total_process_time / 1000.0 << " seconds (" 
              << (total_process_time * 100.0 / total_time) << "%)" << std::endl;
    std::cout << "  - Queue time: " << total_queue_time / 1000.0 << " seconds (" 
              << (total_queue_time * 100.0 / total_time) << "%)" << std::endl;
    
    if(write_debug){
            std::cout << "\n[sigalign] Output written to:\n"
            << "[sigalign] [sigstring]: " << sig_path << "\n"
            << "[sigalign]       [csv]: " << csv_path << "\n"
            << "[sigalign]     [fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n"
            << "[sigalign]   [metrics]: " << metrics_path << std::endl;
    } else {
        std::cout << "\n[sigalign] Output written to:\n"
                  << "[sigalign]     [fastq]: " << fastq_output_path << (compress_fastq ? ".gz" : "") << "\n";
    }
}
    
    // sigstring format
    std::string to_sigstring() const {
        std::stringstream ss;
        std::string overall = read_type; // e.g., "forward", "reverse", "concatenate", "filtered"

        // Group elements by direction
        std::map<std::string, std::vector<std::reference_wrapper<const seq_element>>> direction_elements;
        for (const auto& elem : sig_elements) {
            direction_elements[elem.direction].push_back(std::ref(elem));
        }
        
        bool first_direction = true;
        // Process each directional group
        for (const auto& [direction, elements] : direction_elements) {
            // If overall read type restricts the output, skip the other direction
            if ((overall == "forward" && direction != "forward") || (overall == "reverse" && direction != "reverse")) {
                 continue;
            }
            if (elements.empty()) {
                continue; // Skip empty groups
            }
            if (!first_direction) {
                ss << "\n";
            }
            first_direction = false;
            
            // Sort elements by order
            std::vector<std::reference_wrapper<const seq_element>> sorted_elements = elements;
            std::sort(sorted_elements.begin(), sorted_elements.end(),
                    [](const auto& a, const auto& b) {
                        return a.get().position.first < b.get().position.first;
                    });
            
            bool first_elem = true;
            for (const auto& elem_ref : sorted_elements) {
                const auto& elem = elem_ref.get();

                if (!first_elem) {
                    ss << "|";
                }
                if(elem.position.first == -1 && elem.position.second == -1){
                    continue; // Skip elements with invalid positions
                }
                first_elem = false;
                ss << elem.class_id << ":"
                << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "0") << ":"
                << elem.position.first << ":"
                << elem.position.second;
            }
            // Append final tag: direction abbreviated as F or R
            std::string dir = (direction == "forward") ? "F" : "R";
            std::string concat = read_type == "concatenate" ? ":C" : "";
            ss << "<" << sequence_length << ":" << sequence_id << ":" << dir << concat << ">";
        }
        // Return empty string if nothing was added
        return ss.str().empty() ? "" : ss.str();
    }
   
    // FASTQ format
    std::string to_fastqa() const {
    std::vector<std::string> dirs;
    if(read_type == "skipped"){
        return ""; // Skip if read type is "skipped"
    }
    if (read_type == "concatenate") {
        dirs = { 
            "forward", 
            "reverse" 
        };
    } else {
        dirs = { 
            read_type 
        };
    }

    std::string all_records;
    for (const auto& dir : dirs) {
        std::vector<std::string> bc_keys;
        std::unordered_map<std::string, std::string> bc_map;
        std::unordered_map<std::string, std::string> bc_dir;
        std::unordered_map<std::string, std::string> cr_map;

        std::string umi, read_seq, read_qual;
        for (const auto& elem : sig_elements) {
            if (!elem.seq.has_value()) continue;
            if (elem.direction != dir) continue;

            if (elem.global_class == "barcode") {
                auto key = seq_utils::remove_rc(elem.class_id);
                bool is_fwd = (dir == "forward");

                if (bc_map.find(key) == bc_map.end()) {
                    bc_keys.push_back(key);
                    bc_map[key] = elem.seq.value();
                    bc_dir[key] = dir;
                    if (elem.original_seq.has_value()) {
                        cr_map[key] = elem.original_seq.value();
                    }
                } else if (is_fwd && bc_dir[key] == "reverse") {
                    bc_map[key] = elem.seq.value();
                    bc_dir[key] = dir;
                    if (elem.original_seq.has_value()) {
                        cr_map[key] = elem.original_seq.value();
                    }
                }
                continue;
            }

            if (elem.global_class == "umi") {
                if (elem.seq.has_value()) {
                    umi = elem.seq.value();
                }
                continue;
            }

            if (elem.global_class == "read") {
                read_seq = elem.seq.value();
                if (elem.qual.has_value()) {
                    read_qual = elem.qual.value();
                }
                continue;
            }

            if (elem.global_class == "poly_tail" || elem.global_class == "start" || elem.global_class == "stop") {
                continue;
            }
        }

        std::string cb_tag, cr_tag;
        for (size_t i = 0; i < bc_keys.size(); ++i) {
            const auto& key = bc_keys[i];
            if (i) {
                cb_tag += '-';
                if (!cr_tag.empty()) cr_tag += '-';
            }
            cb_tag += bc_map[key];
            if (cr_map.find(key) != cr_map.end()) {
                cr_tag += cr_map[key];
            }
        }

        bool is_fastq = !read_qual.empty();
        bool is_concatenate = (read_type == "concatenate");
        bool is_forward = (dir == "forward");
        std::stringstream ss;
        //generating modified sequence id for rad
        ss << (is_fastq ? '@' : '>') << sequence_id << (is_forward ? "-F" : "-R") << (is_concatenate ? "-CT" : "");
        //adding barcode tag
        if (!cb_tag.empty()) ss << "\tCB:Z:" << cb_tag;
        //adding corrected read tag for SAM
        if (!cr_tag.empty() && (cb_tag != cr_tag)){
            ss << "\tCR:Z:" << cr_tag;
        } else {
            // Ensure CR tag is present even if empty
            ss << "\tCR:Z:"; 
        }
        //adding transcript tag for SAM
        if (!umi.empty()) ss << "\tUB:Z:" << umi;

        if(is_forward){
            ss << "\tTS:A:+";
        } else {
            ss << "\tTS:A:-";
        }

        ss << "\n" << read_seq << "\n";
        if (is_fastq) {
            ss << "+\n" << read_qual << "\n";
        }
        all_records += ss.str();
    }

    return all_records;
}

    // CSV conversion
    std::string to_csv(bool write_header = false) const {
        std::stringstream ss;
        if (write_header) {
            ss << "id,elem,read_type,seq\n";
        }
        for (const auto& elem : sig_elements) {
            // Process only variable elements that have a non-empty sequence.
            if (elem.type == "variable" && elem.seq.has_value() && !elem.seq->empty()) {
                // Determine the element's direction.
                std::string read_mintype = "forward";
                std::string seq = elem.seq.value();
                std::string final_id = elem.class_id;
                if (elem.class_id.substr(0, 3) == "rc_") {
                    seq = seq_utils::revcomp(seq);
                    final_id = elem.class_id.substr(3);
                    read_mintype = "reverse";
                }
                // If the overall read type is "forward", skip reverse elements.
                if (read_type == "forward" && read_mintype != "forward")
                    continue;
                // If the overall read type is "reverse", skip forward elements.
                if (read_type == "reverse" && read_mintype != "reverse")
                    continue;
                if(read_type == "concatenate"){
                    read_mintype = read_mintype + "_concatenate";
                }
                // If read_type is "concatenate", include both.
                ss << sequence_id << "," 
                << final_id << "," 
                << read_mintype << "," 
                << seq << "\n";
            }
        }
        return ss.str();
    }

    //add split per barcode on demux
    //igblast->airr.tsv (That's a little downstream)
    //immcantation framework is separable by pipes (made into data.table, col.name is in the info)
    //use equals for immcantation (depends on how it gets parsed, but the colon has other meanings)
    //can you encode/decode stuff in a .fasta/.fastq file 

};