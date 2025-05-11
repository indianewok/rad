#pragma once
#include "rad_headers.h"

/**
 * @brief Represents a processed element with alignment information
 */
struct seq_element {
    std::string class_id;          // Unique identifier for the element
    std::string global_class;      // Global classification
    std::optional<int> edit_distance;  //Edit distance from alignment
    std::pair<int, int> position;  // Start and stop positions
    std::string type;              // Element type
    int order;                     // Position in layout
    std::string direction;         // Orientation
    std::optional<bool> element_pass;  // Whether element passed validation
    std::optional<std::string> seq;    // Element sequence if available
    std::optional<std::string> qual;   // Element quality scores if available

    seq_element(
        std::string class_id,
        std::string global_class,
        std::optional<int> edit_distance,
        std::pair<int, int> position,
        std::string type,
        int order,
        std::string direction,
        std::optional<bool> element_pass = std::nullopt,
        std::optional<std::string> seq = std::nullopt,
        std::optional<std::string> qual = std::nullopt
    ) : class_id(std::move(class_id)),
        global_class(std::move(global_class)),
        edit_distance(edit_distance),
        position(position),
        type(std::move(type)),
        order(order),
        direction(std::move(direction)),
        element_pass(element_pass),
        seq(std::move(seq)),
        qual(std::move(qual)) {}
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
static_alignments align_static_elements(const std::string& query, const std::string& target,
    bool verbose,
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
                // Return immediately without SSW.
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
                    std::cout << "SSW alignment: region ["
                            << ssw_start << ", " << ssw_end << "], length = "
                            << ssw_length << ", max_matches with indels = " 
                            << ssw_max_matches_ind
                            << " , max_matches = " << ssw_max_matches
                            << ", deviation = " << deviation
                            << ", cigar: " << sswAlign.cigar_string << "\n";
                }
            }
            // Accept candidate if it meets threshold. first case is within the deviation allowed. second case is for concats.
            if ((ssw_max_matches > min_match_bases && deviation <= 100 && ssw_length >= 10 && ssw_max_matches_ind > 8)||
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

    static_alignments find_poly_tails(const std::string& query, const std::string& sequence, int window_size = 12) {
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
    // returns the best candidate compared to the query candidate
    std::optional<int64_seq> resolve_multiple_hits(const int64_seq& query,
         const std::unordered_set<int64_seq>& candidates, int max_dist, bool verbose){
        auto sorted = mutation_tools::int64_lvdist(query, candidates, max_dist);
        if (sorted.empty()){
            if(verbose){
                #pragma omp critical
                std::cout << "NO_MATCH\n";
            }
            return std::nullopt;
        }
         // lowest edit distance
        auto it = sorted.begin();
        int best_dist = it->first;
        const auto& best_set = it->second;

        if(verbose){
            #pragma omp critical
            std::cout << "Best edit distance = " << best_dist << ", # of candidates = " << best_set.size() << "\n";
        }
    
        if (best_set.size() == 1) {
            if(verbose){
                #pragma omp critical
                std::cout << "Best candidate = " << best_set.begin()->bits_to_sequence() << "\n";
            }
            return *best_set.begin();
        }
        if(verbose){
            #pragma omp critical
            std::cout << "Ambiguous best hits:\n";
            for (auto &s : best_set) {
              std::cout << "  " << s.bits_to_sequence() << "\n";
            }
        }
        return std::nullopt;
    }

    // Returns true if we can correct this one barcode element into exactly one “true” whitelist entry.
    std::optional<int64_seq> correct_barcode(const seq_element& elem, const ReadLayout& layout, bool verbose){
        // pick the right whitelist
        auto key = seq_utils::remove_rc(elem.class_id);
        auto &wl = layout.wl_map.maps.at(key).get();
        // extract and reverse‐complement the raw string
        std::string raw = elem.seq.value();
        if (elem.direction == "reverse"){
            raw = seq_utils::revcomp(raw);
        }
        // encoded barcode
        int64_seq bc;
        bc.sequence_to_bits(raw);
        int max_shift_dist = static_cast<int>(bc.length*0.25);
        int max_resolve_dist = static_cast<int>(bc.length*0.2);
        int max_mutation_dist = static_cast<int>(bc.length*0.6);
        // === filter match ===
        if(!wl.filter_bcs.empty() && wl.filter_bcs.check_wl_for(bc)){
            if(verbose){
                #pragma omp critical
                std::cout << "FILTER_CHECK_FOUND\n";
                std::cout << "NO_CHECK_WORKED\n";
            }
            return std::nullopt;
        }
        // === exact match ===
        if (!wl.true_bcs.empty() && wl.true_bcs.check_wl_for(bc)) {
            auto matched = wl.true_bcs.return_putative_correct_bcs(bc);
            if(verbose){
                #pragma omp critical
                std::cout << "ORIGINAL_CHECK_FOUND (" << matched.size()
                        << " candidates)\n";
            }
            if (matched.size() == 1) {
                if(verbose){
                    #pragma omp critical
                    std::cout << "ORIGINAL_CHECK_WORKED\n";    
                }
              return std::optional<int64_seq>(*matched.begin()); 
            }
            if (auto one = resolve_multiple_hits(bc, matched, max_resolve_dist, verbose)) {
                if(verbose){
                    #pragma omp critical
                    std::cout << "ORIGINAL_MATCH_COLLISION_RESOLVED\n";    
                }
                
                return std::optional<int64_seq>(*one);
            }
            if(verbose){
                #pragma omp critical
                std::cout << "ORIGINAL_MATCH_COLLISION_UNRESOLVED\n";
            }
        }
        if(verbose){
            #pragma omp critical
            std::cout << "ORIGINAL_MATCH_NOT_FOUND\n";
        }
        // === global whitelist ===
        {
            if(!wl.global_bcs.empty() && wl.global_bcs.check_wl_for(bc)){
                auto matched = wl.global_bcs.return_putative_correct_bcs(bc);
                if(verbose){
                    #pragma omp critical
                    std::cout << "GLOBAL_CHECK_FOUND (" << matched.size() << " candidates)\n";
                }
                if (matched.size() == 1) {
                    if(verbose){
                        #pragma omp critical
                        std::cout << "GLOBAL_CHECK_WORKED\n";    
                    }
                    return std::optional<int64_seq>(*matched.begin());
                }
                if(verbose){
                    #pragma omp critical
                    std::cout << "GLOBAL_MATCH_COLLISION_UNRESOLVED\n";    
                }
            }
        }
        // === mutation match ===
        {
        auto muts = mutation_tools::generate_mutated_barcodes(bc, 2);
            if (wl.true_bcs.check_wl_for(muts)) {
                auto matched = wl.true_bcs.return_putative_correct_bcs(muts);
                if(verbose){
                    #pragma omp critical
                    std::cout << "MUTATION_CHECK_FOUND (" << matched.size() << " candidates)\n";
                }
                if (matched.size() == 1) {
                    int64_seq putative_candidate = *matched.begin();
                    putative_candidate.bits_to_sequence();
                    if(verbose){
                        #pragma omp critical
                        std::cout << "MUTATION_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\n";
                    }
                    int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_mutation_dist);
                    if(res >= 0 && res <= max_resolve_dist){
                        if(verbose){
                            #pragma omp critical
                            std::cout << "LVDIST::" << res << "\n";
                            std::cout << "MUTATION_CHECK_MATCHED\n";    
                        }
                        return std::optional<int64_seq>(*matched.begin());
                    }
                }
                if (auto one = resolve_multiple_hits(bc, matched, max_mutation_dist, verbose)) {
                    if(verbose){
                        #pragma omp critical
                        std::cout << "WL_MULTIPLE_MATCHED_RESOLVED\n";
                    }
                    return std::optional<int64_seq>(*one);
                }
                if(verbose){
                    #pragma omp critical
                    std::cout << "WL_MULTIPLE_MATCHED_FAIL\n";
                }
            }
        }
        // === shift fuzzy match ===
        {
        auto shifts = mutation_tools::generate_shifted_barcodes(bc, max_shift_dist);
            if (wl.true_bcs.check_wl_for(shifts)) {
                auto matched = wl.true_bcs.return_putative_correct_bcs(shifts);
                if(verbose){
                    #pragma omp critical
                    std::cout << "SHIFT_CHECK_WORKED (" << matched.size()
                            << " candidates)\n";    
                }
                if (matched.size() == 1) {
                    int64_seq putative_candidate = *matched.begin();
                    int res = mutation_tools::int64_lvdist(bc, putative_candidate, max_mutation_dist);
                    if(verbose){
                        #pragma omp critical
                        std::cout << "SHIFT_CHECK_CANDIDATE::" << putative_candidate.bits_to_sequence() << "\n";
                        std::cout << "LVDIST::" << res << "\n";    
                    }
                    if(res >= 0 && res <= max_resolve_dist){
                        if(verbose){
                            #pragma omp critical
                            std::cout << "SHIFT_CHECK_MATCHED\n";    
                        }
                        return std::optional<int64_seq>(*matched.begin());
                    }
                }

                if (auto one = resolve_multiple_hits(bc, matched, max_resolve_dist, verbose)) {  
                        if(verbose){
                            #pragma omp critical
                            std::cout << "SHIFT_MULTIPLE_MATCHED_RESOLVED\n";
                        }
                    return std::optional<int64_seq>(*one);
                }
                if(verbose){
                    #pragma omp critical
                    std::cout << "SHIFT_MULTIPLE_MATCHED_FAIL\n";
                }
            }
        }
        // === no match ===
        if(verbose){
            #pragma omp critical
            std::cout << "NO_CHECK_WORKED\n";
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
    void clear() { sig_elements.clear(); }

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
                std::cout << "Processing variable element: " << layout_elem->class_id << std::endl;
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
    bool map_positions(
        const ReferencePositions& ref_pos, const std::multimap<std::string, const seq_element*>& static_refs,
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
                    std::cout << "Added static reference " << elem->class_id 
                              << " (order: " << elem->order 
                              << ", pos: " << elem->position.first << "-" << elem->position.second 
                              << ") to ordered list." << std::endl;
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
            std::cout << "Ordered static references for variable element " 
                      << var_it->class_id << ": ";
            for (const auto& elem : ordered_elements) {
                std::cout << elem->class_id << "(" << elem->order << ") ";
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
                    std::cout << "Found " << ref_pos.primary_start.ref_id 
                              << " as the primary start reference candidate." << std::endl;
                    std::cout << "Primary start details: is_start=" << ref_pos.primary_start.is_start 
                              << ", offset=" << ref_pos.primary_start.offset 
                              << ", static pos=(" << primary_candidate->position.first << ","
                              << primary_candidate->position.second << ")" << std::endl;
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
                        std::cout << "Primary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ", global_class: " << (*prev)->global_class << ")." << std::endl;
                    }
                }
                if ((*prev)->global_class != "poly_tail" &&
                    (*prev)->position.second > primary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << "Primary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_candidate->position.first << ")." << std::endl;
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Primary start candidate is the first element in the ordered list." << std::endl;
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
                        std::cout << "Using primary start mapping: calculated start_pos = " << start_pos << std::endl;
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::cout << "Primary start reference " << ref_pos.primary_start.ref_id 
                              << " not found in static_refs." << std::endl;
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
                    std::cout << "Found " << ref_pos.secondary_start.ref_id 
                              << " as the secondary start reference candidate." << std::endl;
                    std::cout << "Secondary start details: is_start=" << ref_pos.secondary_start.is_start 
                              << ", offset=" << ref_pos.secondary_start.offset 
                              << ", static pos=(" << secondary_candidate->position.first << ","
                              << secondary_candidate->position.second << ")" << std::endl;
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Secondary start: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                    }
                }
                if ((*prev)->position.second > secondary_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << "Secondary start out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << secondary_candidate->position.first << ")." << std::endl;
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Secondary start candidate is the first element in the ordered list." << std::endl;
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
                        std::cout << "Using secondary start mapping: calculated start_pos = " << start_pos << std::endl;
                    }
                }
            } 
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::cout << "Secondary start reference " << ref_pos.secondary_start.ref_id 
                              << " not found in static_refs." << std::endl;
                }
            }
        if (!ref_pos.secondary_start.add_flags.empty() && ref_pos.secondary_start.add_flags == "left_truncated") {
            if (verbose) {
                #pragma omp critical
                {
                    std::cout << "Secondary start flagged as left_truncated; attempting fallback." << std::endl;
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
                                std::cout << "Using left terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with start_pos = " << start_pos << std::endl;
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::cout << "Fallback static reference for secondary start not found." << std::endl;
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << "No suitable fallback static element found for secondary start." << std::endl;
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
                    std::cout << "Found " << ref_pos.primary_stop.ref_id 
                              << " as the primary stop reference candidate." << std::endl;
                    std::cout << "Primary stop details: is_start=" << ref_pos.primary_stop.is_start 
                              << ", offset=" << ref_pos.primary_stop.offset 
                              << ", static pos=(" << primary_stop_candidate->position.first << ","
                              << primary_stop_candidate->position.second << ")" << std::endl;
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
                        std::cout << "Primary stop: Previous element in ordered list is " 
                                  << (*prev)->class_id << " (pos: " 
                                  << (*prev)->position.first << "-" << (*prev)->position.second 
                                  << ")." << std::endl;
                    }
                }
                if ((*prev)->position.second > primary_stop_candidate->position.first) {
                    is_ordered = false;
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << "Primary stop out-of-order: previous element's end (" 
                                      << (*prev)->position.second 
                                      << ") > candidate's beginning (" 
                                      << primary_stop_candidate->position.first << ")." << std::endl;
                        }
                    }
                }
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Primary stop candidate is the first element in the ordered list." << std::endl;
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
                        std::cout << "Using primary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::cout << "Primary stop reference " << ref_pos.primary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                }
            }
        }
    }
    
    // If primary stop mapping failed, try secondary.
    if (stop_pos <= 0) {
        if(verbose){
            #pragma omp critical
            {
                std::cout << "Primary stop mapping failed; attempting secondary stop mapping." << std::endl;
                std::cout << "Secondary stop reference ID: " << ref_pos.secondary_stop.ref_id << std::endl;
                std::cout << "Secondary flags: " << ref_pos.secondary_stop.add_flags << std::endl;
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
                    std::cout << "Found " << ref_pos.secondary_stop.ref_id 
                              << " as the secondary stop reference candidate." << std::endl;
                    std::cout << "Secondary stop details: is_start=" << ref_pos.secondary_stop.is_start 
                              << ", offset=" << ref_pos.secondary_stop.offset 
                              << ", static pos=(" << secondary_stop_candidate->position.first << ","
                              << secondary_stop_candidate->position.second << ")" << std::endl;
                }
            }
            bool is_ordered = true;
            if (ref_order_pos != ordered_elements.begin()) {
                auto prev = std::prev(ref_order_pos);
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Secondary stop: Previous element in ordered list is " 
                                << (*prev)->class_id << " (pos: " 
                                << (*prev)->position.first << "-" << (*prev)->position.second 
                                << ")." << std::endl;
                    }
                }
                if ((*prev)->position.second > secondary_stop_candidate->position.first)
                    is_ordered = false;
            } else {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Secondary stop candidate is the first element in the ordered list." << std::endl;
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
                        std::cout << "Using secondary stop mapping: calculated stop_pos = " << stop_pos << std::endl;
                    }
                }
            }
        } else {
            if (verbose) {
                #pragma omp critical
                {
                    std::cout << "Secondary stop reference " << ref_pos.secondary_stop.ref_id 
                              << " not found in static_refs." << std::endl;
                }
            }
            if (!ref_pos.secondary_stop.add_flags.empty() && ref_pos.secondary_stop.add_flags == "right_truncated") {
                if (verbose) {
                    #pragma omp critical
                    {
                        std::cout << "Secondary stop flagged as right_truncated; attempting fallback." << std::endl;
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
                                std::cout << "Using right terminal linked fallback: " 
                                            << fallback_candidate->class_id 
                                            << " with stop_pos = " << stop_pos << std::endl;
                            }
                        }
                    } else {
                        if (verbose) {
                            #pragma omp critical
                            {
                                std::cout << "Fallback static reference for secondary stop not found." << std::endl;
                            }
                        }
                    }
                } else {
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << "No suitable fallback static element found for secondary stop." << std::endl;
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
            std::cout << "Final calculated positions: " << start_pos << ":" << stop_pos 
                      << " for variable element " << var_it->class_id << std::endl;
        }
    }
    
    // Validate positions relative to read length.
    if (validate_var_positions(var_positions, read_seq.length())) {
        if (verbose) {
            #pragma omp critical
            {
                std::cout << "Positions validated for variable element " << var_it->class_id << std::endl;
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
        std::cout << "Position validation failed for variable element " << var_it->class_id << std::endl;
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
                std::cout << "  ["<<direction<<"] invalid positions on “"
                        << e.class_id << "”\n";
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
                        std::cout << "Overlap: " 
                                  << e1.class_id << " vs " << e2.class_id << "\n";
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

    // update barcode counts
    void update_bc_counts(const SigString &sig, const ReadLayout &layout, bool verbose = false) {
        auto type = sig.type(); // “forward”, “reverse”, “concatenate”, or “filtered”
        for (auto const &elem : sig.elements()) {
            // ||  !elem.element_pass.value_or(false)
            if (elem.global_class != "barcode" || !elem.seq.has_value()) {
                continue;
            }
        // recover the raw sequence (reverse‐comp if needed)
        std::string bc_seq = elem.seq.value();
        if (elem.direction == "reverse" && type != "forward"){
            bc_seq = seq_utils::revcomp(bc_seq);
        }
        int64_seq bc(bc_seq);

        // look up the right whitelist
        auto key = seq_utils::remove_rc(elem.class_id);
        auto &wl = layout.wl_map.maps.at(key).get();

        // always bump the total counter
        if(verbose){
            #pragma omp critical
            {
                std::cout << "UPDATE_BC_COUNT_TYPE " << type << " FOR ELEM.DIR: " << elem.direction << std::endl;
            }
        }
        if(type == "filtered" || type != elem.direction){
            wl.true_bcs.update_bc_count(bc, barcode_counts::filtered);
        }
        if (type == "forward" && elem.direction == "forward") {
                wl.true_bcs.update_bc_count(bc, barcode_counts::total);
                wl.true_bcs.update_bc_count(bc, barcode_counts::forw);
            } else if (type == "reverse" && elem.direction == "reverse") {
                wl.true_bcs.update_bc_count(bc, barcode_counts::total);
                wl.true_bcs.update_bc_count(bc, barcode_counts::rev);
            } else if (type == "concatenate") {
                if (elem.direction == "forward" && elem.direction == "forward") {
                    wl.true_bcs.update_bc_count(bc, barcode_counts::total);
                    wl.true_bcs.update_bc_count(bc, barcode_counts::forw_concat);
                }
                if(elem.direction == "reverse" && elem.direction == "reverse"){
                    wl.true_bcs.update_bc_count(bc, barcode_counts::total);
                    wl.true_bcs.update_bc_count(bc, barcode_counts::rev_concat);
                }
            }
            int bc_count = wl.true_bcs.get_bc_count(bc);
            if(verbose){
                #pragma omp critical
                std::cout << "Barcode count for " << bc.bits_to_sequence() << ": " << bc_count << std::endl;
            }
        }
    }

public:
    //sigalign_static v2
   void sigalign_static(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
    aligner_tools aligner;
    auto& type_index = layout.by_type();
    auto static_range = type_index.equal_range("static");

    if(verbose){
        #pragma omp critical
        {
            std::cout << "\n=== Starting static alignment for " << sequence_id << " ===" << std::endl;
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
                std::cout << "Processing static element: " << it->class_id << std::endl;
                if(it->aligned_positions){
                    std::cout << "  Aligned positions: " 
                              << it->aligned_positions->start_stats.first << ", "
                              << it->aligned_positions->start_stats.second << std::endl;
                } else {
                    std::cout << "  No aligned positions available." << std::endl;
                }
                if(it->misaligned_positions){
                    std::cout << "  Misaligned positions: " 
                              << it->misaligned_positions->start_stats.first << ", "
                              << it->misaligned_positions->start_stats.second << std::endl;
                } else {
                    std::cout << "  No misaligned positions available." << std::endl;
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
            auto result = aligner.find_poly_tails(it->seq, read.seq, 12);
            if (result.success) {
                for (const auto& positions : result.positions) {
                    add_element(seq_element(
                        it->class_id,
                        it->global_class,
                        result.edit_distance,
                        positions,
                        "static",
                        it->order,
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

                std::cout << "Expected region for " << it->class_id << ": " 
                          << expected_start << " - " << expected_end << std::endl;
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
                std::cout << "\n=== Starting variable alignment for " << sequence_id << " ===" << std::endl;
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
                        std::cout << "Found static reference: " << elem.class_id 
                                  << " at position " << elem.position.first 
                                  << ":" << elem.position.second 
                                  << " with edit distance " << (elem.edit_distance ? std::to_string(elem.edit_distance.value()) : "none")
                                  << " and sequence " << (elem.seq ? elem.seq.value() : "none") 
                                  << std::endl;
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
                    std::cout << "Found reference: " << elem.class_id 
                              << " at position " << elem.position.first 
                              << ":" << elem.position.second << std::endl;
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
                std::cout << "Successfully positioned " << positioned_count << " variable elements" << std::endl;
                std::cout << "Final sigstring elements: " << sig_elements.size() << std::endl;
            }
        }
    }
    // contains per-unit barcode correction
    void sigalign_filter(const read_streaming::sequence &read, const ReadLayout& layout, bool verbose) {
        // make a resever for the qual that doesn't interfere w/ASCII characters
        constexpr char qual_mask = '\x7F';
        auto direction_elements = group_directionally();
        // For recording validity and count per direction.
        std::map<std::string, bool> direction_valid;
        std::map<std::string, int> pass_counts, static_counts;
        std::unordered_map<std::string, std::unordered_set<int64_seq>> seen_bcs;
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
                        std::cout << "Skipping a read element due to out-of-bounds parameters." << std::endl;
                        std::cout << "Start: " << start << ", Length: " << length 
                                  << ", Read length: " << masked_read.size() << std::endl;
                        std::cout << "Masked read: " << masked_read << std::endl;
                        std::cout << "Direction: " << direction << std::endl;
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
                            std::cout << "Masked qual: " << el.qual.value() << "\n";
                        }
                    } else {
                        //el.qual.reset();
                        std::cout << "Masked qual: " << el.qual.value() << "\n";
                    }
                });

                // verbose logging
                if (verbose) {
                    #pragma omp critical
                    std::cout << "Trimmed read for " << e.class_id 
                              << " -> " << cleaned_seq << "\n";
                    std::cout << "Masked read: " << masked_read << "\n";
                    if(read.is_fastq) {
                        std::cout << "Masked qual: " << cleaned_qual << "\n";
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
                        std::cout << "FILTERED_READ_LENGTH" << std::endl;
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
                        id_index.find(elem.class_id),
                                      [](seq_element &e) { 
                                        e.element_pass = false; 
                                    });
                    if (verbose) {
                        #pragma omp critical
                        {
                            std::cout << 
                            "Filtered element " << 
                            elem.class_id << 
                            " (invalid positions) in direction " << direction << ".\n";
                        }
                    }
                    // Exit processing for this direction.
                    break;
                }

                //barcode correction here
                if (elem.global_class == "barcode") {
                    if(verbose){
                        #pragma omp critical
                        {
                            std::cout << "Barcode correction for " << elem.seq.value() << std::endl;
                        }
                    }
                    auto& wl = layout.wl_map.maps.at(seq_utils::remove_rc(elem.class_id)).get();
                    //had originally included a check for global barcode class not empty (&& !wl.true_bcs.empty())
                    //but this is not necessary, as the barcode correction function will handle it
                    if(!wl.true_bcs.empty()){
                        auto out = barcode_correction::correct_barcode(elem, layout, verbose);
                        bool corrected = out.has_value();
                        if(corrected) {
                            int64_seq correct_bc = out.value();
                            seen_bcs[seq_utils::remove_rc(elem.class_id)].insert(out.value());
                            std::string final_bc = elem.direction == "forward" ? correct_bc.bits_to_sequence() : 
                            seq_utils::revcomp(correct_bc.bits_to_sequence());
                            // add to corrected barcode count
                            if(final_bc !=elem.seq.value()) {
                                wl.true_bcs.update_bc_count(correct_bc, barcode_counts::corrected);
                            }
                            auto& id_index = sig_elements.get<sig_id_tag>();
                            id_index.modify(
                                id_index.find(elem.class_id),
                                [&](seq_element &e){
                                    e.seq = final_bc;
                                    e.element_pass = true;
                                }
                            );
                        } else {
                            auto& id_index = sig_elements.get<sig_id_tag>();
                            id_index.modify(
                                id_index.find(elem.class_id),
                                [](seq_element &e) {
                                    e.element_pass = false; 
                                    }
                                );
                                if(elem.seq.has_value()){
                                    std::string seq = elem.seq.value();
                                    int64_seq putative_bc(seq);
                                    barcode_entry be;
                                    be.barcode = putative_bc;
                                    be.flags = "";
                                    be.edit_dist = 0;
                                    be.filtered = true;
                                   wl.filter_bcs.insert_bc_entry(be.barcode, be);
                                }
                            // adding this in to try to manage if there are multiple barcodes, 
                            // whether we should set this entire direction to false. 
                            // right now, we *only* return completely correct barcodes (all multiples are also correct) 
                            if(!multiple_barcodes){
                                // so if multiple barcodes is false, then we set the direction value to false
                                // and set the pass counts to 0.
                                valid_direction = false;
                                pass_counts[direction] = 0;    
                            }
                            continue;
                        }
                    } else {
                        //this branch is for generating a purely standalone whitelist that will occur IFF there is no default whitelist
                        //--default whitelist generation has yet to happen considering just dumping them into the true_bcs. 
                        //will probably pull this out and write a standalone sorting function internally here.
                        if(elem.seq.has_value()){
                            std::string seq = elem.seq.value();
                            int64_seq putative_bc(seq);
                            barcode_entry be;
                            be.barcode = putative_bc;
                            be.flags = "";
                            be.edit_dist = 0;
                            be.filtered = true;
                           wl.filter_bcs.insert_bc_entry(be.barcode, be);
                        }
                    }
                }
            }

            // If no static elements were found, mark the direction as invalid.
            if(static_counts[direction] == 0){
                valid_direction = false;
                pass_counts[direction] = 0;
            }

            // If the direction failed the per-element check, record zero passing variable elements.
            if (!valid_direction) {
                direction_valid[direction] = false;
                pass_counts[direction] = 0;
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
                                std::cout << "Overlap detected: Filtering variable elements " 
                                          << e1.class_id << " and " << e2.class_id 
                                          << " in direction " << direction << ".\n";
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
                    std::cout << "\nDirection " << direction << " valid: " << direction_valid[direction] 
                              << " with " << count << " passing variable elements.\n";
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
                std::cout << "Final read type: "
                          << read_type
                          << ", forward count: "
                          << forward_count
                          << ", reverse count: "
                          << reverse_count
                          << "\n"
                          << std::endl;
            }
        }
    }
    // full sigalign implementation (took out const on read layout in sigalign filter->whitelist population)
    static void sigalign(const std::string& fastq_path, const ReadLayout& layout, const std::string& output_prefix, bool verbose, int num_threads = 1, 
        size_t max_reads = -1, bool split_bc = false) {
        //std::string home = std::getenv("HOME");
        //std::string desktop_path = home + "/Desktop/";
        std::string sig_path = output_prefix + ".sig";
        std::string csv_path = output_prefix + ".csv";
        std::string fastq_output_path = output_prefix + ".fq";

        sigstring_writing sig_writer(sig_path, sigstring_writing::format::SIGSTRING, /*compress=*/false, /*append=*/false);
        sigstring_writing csv_writer(csv_path, sigstring_writing::format::CSV, /*compress=*/true, /*append=*/false);
        sigstring_writing fastqa_writer(fastq_output_path, sigstring_writing::format::FASTQA, /*compress=*/true, /*append=*/false);

        // csv writer will emit exactly what to_csv(true) produces
        {
            SigString header("",0);
            csv_writer(std::vector<SigString>{header});
        }

        // Define the chunk processing function
        using ChunkFunc = std::function<void(const std::vector<read_streaming::sequence>&, const std::string&)>;
        chunk_streaming<read_streaming::sequence,ChunkFunc> streamer;
    
        // Set up the chunk processing function to process each chunk of reads
        ChunkFunc process_func = [&](auto const &chunk, auto const & /*unused_file*/) {
            // collect SigString objects
            std::vector<SigString> to_write;
            to_write.reserve(chunk.size());
            for (auto const& read : chunk) {
                SigString sig(read.id, read.seq.length());
                sig.sigalign_static(read, layout, verbose);
                sig.sigalign_variable(read, layout, verbose);
                sig.sigalign_filter(read, layout, verbose);

                if (sig.read_type != "filtered") {
                    to_write.push_back(std::move(sig));
                }
            }
    
            if (!to_write.empty()) {
                #pragma omp critical
                {
                    sig_writer(to_write);
                    csv_writer(to_write);
                    fastqa_writer(to_write);
                }
            }
        };

        // Process the chunks
        streamer.process_chunks(fastq_path, process_func, num_threads, static_cast<int64_t>(max_reads));

        std::cout <<
                   "[sigalign] Output written to:\n"
                << "[sigalign] [sigstring]: " << sig_path << "\n"
                << "[sigalign]       [csv]: " << csv_path << "\n"
                << "[sigalign]     [fastq]: " << fastq_output_path << std::endl;
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
   
    //to fastqa
    std::string to_fastqa() const {
        // 1) decide which directions to emit
        std::vector<std::string> dirs;
        if (read_type == "concatenate") {
            dirs = { "forward", "reverse" };
        } else {
            dirs = { 
                read_type 
            };
        }
    
        std::string all_records;
    
        // 2) for each requested direction, build one record
        for (auto const &dir : dirs) {
            // collect barcodes in insertion order, collapsing RC if needed
            std::vector<std::string> bc_keys;
            std::unordered_map<std::string, std::string> bc_map;
            std::unordered_map<std::string, std::string> bc_dir;
    
            std::string umi, read_seq, read_qual;
    
            for (auto const &elem : sig_elements) {
                if (!elem.seq.has_value()){
                    continue;
                }
    
                // only take elements in this direction
                if (elem.direction != dir) 
                    continue;
    
                // barcode
                if (elem.global_class == "barcode") {
                    auto key = seq_utils::remove_rc(elem.class_id);
                    bool is_fwd = (dir == "forward");
    
                    auto it = bc_map.find(key);
                    if (it == bc_map.end()) {
                        bc_keys.push_back(key);
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                    }
                    else if (is_fwd && bc_dir[key] == "reverse") {
                        bc_map[key] = elem.seq.value();
                        bc_dir[key] = dir;
                    }
                    continue;
                }
    
                // umi 
                if (elem.global_class == "umi") {
                    if(elem.seq.has_value()){
                        umi = elem.seq.value();
                    }
                    continue;
                }
    
                // read and qual
                if (elem.global_class == "read") {
                    read_seq = elem.seq.value();
                    if (elem.qual.has_value()) {
                        read_qual = elem.qual.value();
                    }
                    continue;
                }
    
                if (elem.global_class == "poly_tail" || elem.global_class == "start" || elem.global_class == "stop")
                    {
                        continue;
                }
            }
    
            //add CB tag
            std::string cb_tag;
            for (size_t i = 0; i < bc_keys.size(); ++i) {
                if (i) cb_tag += '-';
                cb_tag += bc_map[bc_keys[i]];
            }
    
            bool is_fastq = !read_qual.empty();
            std::stringstream ss;
            ss << (is_fastq ? '@' : '>') << sequence_id;
            if (!cb_tag.empty()) ss << " CB:Z:" << cb_tag;
            if (!umi.empty())    ss << " UB:Z:" << umi;
            ss << "\n" << read_seq << "\n";
            if (is_fastq){
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

    //SEAMLESSLY COMPATIBLE WITH MINIMAP2
    //add split per barcode on demux
    //igblast->airr.tsv (That's a little downstream)
    //immcantation framework is separable by pipes (made into data.table, col.name is in the info)
    //use equals for immcantation (depends on how it gets parsed, but the colon has other meanings)
    //can you encode/decode stuff in a .fasta/.fastq file 

};