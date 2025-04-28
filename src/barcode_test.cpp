#include "include/rad/barcode_correction_old.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include "edlib/edlib/include/edlib.h"
#include "include/ssw/ssw_cpp.h"
//#include "rapidfuzz/rapidfuzz/rapidfuzz_all.hpp"

// Trim whitespace from a string
std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string::npos) {
        return ""; // String is all whitespace
    }
    
    auto end = s.find_last_not_of(" \t\n\r\f\v");
    return s.substr(start, end - start + 1);
}

// Load barcodes from file
std::vector<std::string> load_barcodes(const std::string& filename) {
    std::vector<std::string> barcodes;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return barcodes;
    }
    
    std::string line;
    // Skip header if present
    if (std::getline(file, line) && (line.find("Barcode") == 0 || line.find("barcode") == 0)) {
        // Skipped header
    } else {
        // This wasn't a header, process it as a barcode
        size_t delimiter_pos = line.find_first_of(",\t\r");
        if (delimiter_pos != std::string::npos) {
            std::string barcode = trim(line.substr(0, delimiter_pos));
            if (!barcode.empty()) {
                barcodes.push_back(barcode);
            }
        } else if (!line.empty()) {
            barcodes.push_back(trim(line));
        }
    }
    
    // Read remaining lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        size_t delimiter_pos = line.find_first_of(",\t\r");
        if (delimiter_pos != std::string::npos) {
            std::string barcode = trim(line.substr(0, delimiter_pos));
            if (!barcode.empty()) {
                barcodes.push_back(barcode);
            }
        } else {
            barcodes.push_back(trim(line));
        }
    }
    
    file.close();
    return barcodes;
}

// Load barcodes and convert to bit representation
std::unordered_map<int64_t, std::string> load_bit_barcodes(const std::string& filename) {
    std::unordered_map<int64_t, std::string> barcode_map;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return barcode_map;
    }
    
    std::string line;
    // Skip header if present
    if (std::getline(file, line) && (line.find("Barcode") == 0 || line.find("barcode") == 0)) {
        // Skipped header
    } else {
        // Process this line as it's not a header
        size_t delimiter_pos = line.find_first_of(",\t\r");
        std::string barcode = delimiter_pos != std::string::npos ? 
                              trim(line.substr(0, delimiter_pos)) : 
                              trim(line);
        
        if (!barcode.empty()) {
            int64_seq seq(barcode);
            if (!seq.bits.empty()) {
                barcode_map[seq.bits[0]] = barcode;
            }
        }
    }
    
    // Read remaining lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        size_t delimiter_pos = line.find_first_of(",\t\r");
        std::string barcode = delimiter_pos != std::string::npos ? 
                              trim(line.substr(0, delimiter_pos)) : 
                              trim(line);
        
        if (!barcode.empty()) {
            int64_seq seq(barcode);
            if (!seq.bits.empty()) {
                barcode_map[seq.bits[0]] = barcode;
            }
        }
    }
    
    file.close();
    return barcode_map;
}

// Edlib comparison
void edlib_compare(const std::string& query, const std::vector<std::string>& barcodes, int k = 3) {
    std::cout << "\n=== Edlib Comparison for '" << query << "' ===" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Set up custom equality for DNA bases
    const int numEq = 13;
    EdlibEqualityPair customEqualities[numEq] = {
        {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
        {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
        {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'},
        {'N','x'}
    };

    
    EdlibAlignConfig config = edlibNewAlignConfig(k, // max edit distance
                                           EDLIB_MODE_HW, // global alignment
                                           EDLIB_TASK_DISTANCE, // just get edit distance
                                           customEqualities, numEq);
    
    // Store results: <barcode, edit_distance>
    std::vector<std::pair<std::string, int>> results;
    std::string padded_barcode;
    // Process each barcode
    for (const auto& barcode : barcodes) {
        padded_barcode = "NNN" + barcode + "NNN";
        EdlibAlignResult result = edlibAlign(query.c_str(), query.length(),
                                            padded_barcode.c_str(), padded_barcode.length(),
                                            config);
        
        if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
            results.push_back({padded_barcode, result.editDistance});
        }
        
        edlibFreeAlignResult(result);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Sort by edit distance
    std::sort(results.begin(), results.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Group by edit distance
    std::unordered_map<int, std::vector<std::string>> grouped;
    for (const auto& [barcode, dist] : results) {
        grouped[dist].push_back(barcode);
    }
    
    // Print results
    std::cout << "Found " << results.size() << " matches with edit distance <= " << k << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    std::cout << "Putative time for ~10K unique sequences: " 
              << (elapsed * 10000 / results.size()) << " ms" << std::endl;

}

// Edlib comparison with multiple queries
void edlib_compare_batch(const std::vector<std::string>& queries, 
                         const std::vector<std::string>& barcodes, 
                         int k = 3, bool verbose = false) {
    std::cout << "\n=== Edlib Batch Comparison (" << queries.size() << " queries) ===" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Set up custom equality for DNA bases
    const int numEq = 13;
    EdlibEqualityPair customEqualities[numEq] = {
        {'A','a'}, {'C','c'}, {'T','t'}, {'G','g'},
        {'A','x'}, {'C','x'}, {'T','x'}, {'G','x'},
        {'N','A'}, {'N','C'}, {'N','T'}, {'N','G'},
        {'N','x'}
    };
    
    EdlibAlignConfig config = edlibNewAlignConfig(k, // max edit distance
                                         EDLIB_MODE_HW, // global alignment
                                         EDLIB_TASK_DISTANCE, // just get edit distance
                                         customEqualities, numEq);
    
    // Store results for each query: <query, vector<barcode, edit_distance>>
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, int>>>> all_results;
    int total_matches = 0;
    
    // Process each query barcode
    for (const auto& query : queries) {
        std::vector<std::pair<std::string, int>> query_results;
        std::string padded_barcode;
        
        // Compare against each reference barcode
        for (const auto& barcode : barcodes) {
            padded_barcode = barcode;
            EdlibAlignResult result = edlibAlign(query.c_str(), query.length(),
                                              padded_barcode.c_str(), padded_barcode.length(),
                                              config);
            
            if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
                query_results.push_back({padded_barcode, result.editDistance});
                total_matches++;
            }
            
            edlibFreeAlignResult(result);
        }
        
        // Sort by edit distance
        std::sort(query_results.begin(), query_results.end(), 
                [](const auto& a, const auto& b) { return a.second < b.second; });
        
        all_results.push_back({query, query_results});
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Print overall results
    std::cout << "Processed " << queries.size() << " queries against " << barcodes.size() << " barcodes" << std::endl;
    std::cout << "Found " << total_matches << " total matches with edit distance <= " << k << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    std::cout << "Average time per query: " << (elapsed / queries.size()) << " ms" << std::endl;
    std::cout << "Putative time for ~10K unique sequences: " 
              << ((elapsed / queries.size())*10000) << " ms" << std::endl;
}

// RapidFuzz comparison
void rapidfuzz_compare(const std::string& query, const std::vector<std::string>& barcodes, double threshold = 70.0) {
    std::cout << "\n=== RapidFuzz Comparison for '" << query << "' ===" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Store results: <barcode, score>
    std::vector<std::pair<std::string, double>> results;
    
    // Process each barcode
    for (const auto& barcode : barcodes) {
        // Calculate the partial ratio score (0-100, higher is better match)
        auto score_alignment = rapidfuzz::fuzz::partial_ratio_alignment(query, barcode, 0.8);
    
        // Extract the score from the alignment result
        double score = score_alignment.score;
        results.push_back({barcode, score});
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Sort by score (descending)
    std::sort(results.begin(), results.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Group by score ranges (100, 95-99, 90-94, etc.)
    std::unordered_map<int, std::vector<std::string>> grouped;
    for (const auto& [barcode, score] : results) {
        int bucket = static_cast<int>(score) / 5 * 5; // Round down to nearest 5
        grouped[bucket].push_back(barcode);
    }
    
    // Print results
    std::cout << "Found " << results.size() << " matches with score >= " << threshold << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    
    // Show scores by range (descending)
    std::vector<int> score_buckets;
    for (const auto& [bucket, _] : grouped) {
        score_buckets.push_back(bucket);
    }
    std::sort(score_buckets.begin(), score_buckets.end(), std::greater<int>());
    
    for (int bucket : score_buckets) {
        auto& matches = grouped[bucket];
        std::cout << "Score " << bucket << "-" << (bucket + 4) << " (" << matches.size() << " matches):" << std::endl;
        
        if (matches.size() <= 10) {
            // Print all matches if fewer than 10
            for (const auto& barcode : matches) {
                // Find the exact score for this barcode
                auto it = std::find_if(results.begin(), results.end(), 
                                     [&barcode](const auto& p) { return p.first == barcode; });
                double exact_score = it->second;
                std::cout << "  " << barcode << " (score: " << exact_score << ")" << std::endl;
            }
        } else {
            // Print first 5 and last 5 if more than 10
            for (int j = 0; j < 5; ++j) {
                auto it = std::find_if(results.begin(), results.end(), 
                                     [&](const auto& p) { return p.first == matches[j]; });
                double exact_score = it->second;
                std::cout << "  " << matches[j] << " (score: " << exact_score << ")" << std::endl;
            }
            std::cout << "  ... (" << (matches.size() - 10) << " more) ..." << std::endl;
            for (int j = matches.size() - 5; j < matches.size(); ++j) {
                auto it = std::find_if(results.begin(), results.end(), 
                                     [&](const auto& p) { return p.first == matches[j]; });
                double exact_score = it->second;
                std::cout << "  " << matches[j] << " (score: " << exact_score << ")" << std::endl;
            }
        }
    }
}

//in-house edit distance tool
// Code for fast edit distance calculation for short sequences modified from
// s2 is always assumed to be the shorter string (barcode)
uint edit_distance(const std::string& s1, const std::string& s2, uint &end, int max_editd){
    std::size_t len1 = s1.size()+1, len2 = s2.size()+1;
    const char * s1_c = s1.c_str(); const char * s2_c = s2.c_str();
    std::vector<uint> dist_holder(len1*len2);
    //initialise the edit distance matrix.
    //penalise for gaps at the start and end of the shorter sequence (j)
    //but not for shifting the start/end of the longer sequence (i,0)
    dist_holder[0]=0; //[0][0]
    for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
    for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0]; 
    int best=len2;
    end=len1-1;
    //loop over the distance matrix elements and calculate running distance
    for (unsigned int j = 1; j < len2; ++j) {
      bool any_below_threshold = false; // flag used for early exit
      for (unsigned int i = 1; i < len1; ++i) {
        int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1; // are the bases the same?
        // if yes, no need to increment distance
        if (sub == 0) {
          dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
        }
        // otherwise add insertion, deletion or substitution
        else {
          // clang-format off
          dist_holder[i*len2+j] = std::min({ //j+i*len2  //set[i][j]
            dist_holder[(i-1)*len2+j] + 1, //[i-1][j]
            dist_holder[i*len2+(j-1)] + 1, //[i][j-1]
            dist_holder[(i-1)*len2+(j-1)] + 1}); // ((s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1) });
          // clang-format on
        }
        if (dist_holder[i * len2 + j] <= max_editd) {
          any_below_threshold = true;
        }
        // if this is the last row in j
        if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
          // check if this is the best running score
          best = dist_holder[i * len2 + j];
          end = i; // update the end position of alignment
        }
      }
      // early exit to save time.
      if(!any_below_threshold) {
        return(100);
      }
    }
    return best; // return edit distance
}

std::vector<std::string> generate_all_shifted_barcodes(
    const std::string &barcode,
    int shift
) {
    std::vector<std::string> out;
    int L = (int)barcode.size();
    if (shift <= 0 || shift >= L) {
        // nothing to shift, or shift >= length → no valid variants
        return out;
    }

    // Core for left‑shifts: drop last `shift` bases
    std::string coreR = barcode.substr(0, L - shift);
    // Core for right‑shifts: drop first  `shift` bases
    std::string coreL = barcode.substr(shift);

    const char bases[4] = {'A','C','G','T'};
    std::string buf(shift, ' ');

    // 1) Build left‑shifted variants: all prefixes + coreR
    std::function<void(int)> dfs_left = [&](int pos) {
        if (pos == shift) {
            out.push_back(buf + coreR);
            return;
        }
        for (char b : bases) {
            buf[pos] = b;
            dfs_left(pos + 1);
        }
    };
    dfs_left(0);

    // 2) Build right‑shifted variants: coreL + all suffixes
    std::function<void(int)> dfs_right = [&](int pos) {
        if (pos == shift) {
            out.push_back(coreL + buf);
            return;
        }
        for (char b : bases) {
            buf[pos] = b;
            dfs_right(pos + 1);
        }
    };
    dfs_right(0);

    return out;
}

std::unordered_set<int64_seq> generate_int64_shifted_barcodes(
    const int64_seq &seq,
    int shift
) {
    std::unordered_set<int64_seq> result;
    // 1) Quick exits
    if (seq.bits.empty() || shift <= 0 || shift >= seq.length)
        return result;

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

//this one works the fastest 
/*
int bit_ld(int64_t p, int64_t t, int n = 16) {
    int64_t np = ~p;
    int64_t HMASK = (1LL << (n - 1));        // Mask to check the high-order bit of the current column
    int64_t VP = (1LL << n) - 1;               // Vertical positive (all ones initially)
    int64_t VN = 0;                          // Vertical negative (all zeros initially)
    int score = n;                           // Initialize the score to n

    for (int j = 0; j < n; ++j) {
        int64_t Bj;
        // Select Bj based on the bit in t at position j.
        if ((t & (1LL << j)) != 0)
            Bj = p;
        else
            Bj = np;

        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);

        if ((HP & HMASK) != 0)
            score += 1;
        else if ((HN & HMASK) != 0)
            score -= 1;

        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    }
    return score;
}
*/

// Myers’s bit‑parallel edit distance with an early cutoff.
// p:  bitmask where bit j=1 means pattern[j] matches '1'
// t:  bitmask for text only used to select p or np per column
// n:  number of columns (<=64)
// max_dist: as soon as even the best possible remaining score exceeds max_dist,
//           we return max_dist+1 as a sentinel.
int bit_ld_cutoff(int64_t p, int64_t t, int n , int max_dist = INT_MAX) {
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
        // early cutoff bound: even if we subtract 1 in every remaining column,
        // our bestPossible = score - (n - j - 1)
        int remaining = n - j - 1;
        if (score - remaining > max_dist)
            return max_dist + 1;
        // advance VP, VN
        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    }
    return score;
}

// Compute the Myers bit‑parallel edit distance between pattern P and text T.
// P may contain 'A','C','G','T','N' (where 'N' matches anything).
// T may also contain 'A','C','G','T','N' (treat 'N' in text as wildcard).
// Returns the full edit distance (insertions, deletions, substitutions).
// Requires P.size() <= 64.
//this one works great!
int bit_lev_wildcard_v1(const std::string &P, const std::string &T) {
    int m = (int)P.size();
    int n = (int)T.size();
    assert(m <= 64 && "Pattern too long for 64‑bit bit‑parallel algorithm");

    // 1) Build Peq: for each possible character (A,C,G,T,'N'), a mask of positions in P.
    //    For 'N' in P, we set that bit for *every* concrete base, so it always matches.
    //    Similarly, when we look up Peq[T[i]], if T[i]=='N' we'll treat it as all-ones (wildcard).
    uint64_t Peq[128] = {0};
    for (int i = 0; i < m; ++i) {
        char c = P[i];
        uint64_t bit = (1ULL << i);
        if (c == 'N') {
            // wildcard in pattern: set this bit for all four bases
            Peq['A'] |= bit;
            Peq['C'] |= bit;
            Peq['G'] |= bit;
            Peq['T'] |= bit;
        } else {
            Peq[(int)c] |= bit;
        }
    }
    // Make sure that if text has 'N', we treat it as matching everything:
    // we'll handle that below by substituting Eq = all_ones when T[i]=='N'.

    uint64_t VP = ~0ULL;              // VP = all 1s for m positions
    uint64_t VN = 0ULL;               // VN = all 0s
    int      currDist = m;            // initial distance = m (all deletions)

    uint64_t mask_msb = (1ULL << (m - 1));   // mask for the "leftmost" bit

    for (int j = 0; j < n; ++j) {
        char c = T[j];
        // 2) Lookup Eq mask for this text character:
        uint64_t Eq;
        if (c == 'N') {
            // wildcard in text matches any pattern position
            Eq = (m == 64 ? ~0ULL : ((1ULL << m) - 1));
        } else {
            Eq = Peq[(int)c];
        }

        // 3) Myers bit‑parallel update
        uint64_t X  = Eq | VN;
        uint64_t D0 = (((X & VP) + VP) ^ VP) | X;
        uint64_t HN = VP & D0;
        uint64_t HP = VN | ~(VP | D0);

        // 4) Update distance based on the msb
        if (HP & mask_msb) {
            currDist++;
        } else if (HN & mask_msb) {
            currDist--;
        }

        // 5) Advance VP, VN for next character
        uint64_t HP_shift = (HP << 1) | 1ULL;
        VN = HP_shift & D0;
        VP = (HN << 1) | ~(HP_shift | D0);
    }

    return currDist;
}

int bit_ld(int64_t p, int64_t t, int n = 16) {
    // 1) Build a one‑bit mask: bit j == 1 if and only if
    //    the jᵗʰ base in p equals the jᵗʰ base in t.
    int64_t Peq = 0;
    for (int j = 0; j < n; ++j) {
        // Extract the 2 bits for base j (LSB = position 0)
        int shift = 2 * j;
        int bp = (p >> shift) & 0b11;
        int bt = (t >> shift) & 0b11;
        if (bp == bt) {
            Peq |= (1LL << j);
        }
    }

    // 2) Run Myers’s bit‑parallel loop over that equality mask.
    int64_t p_eq  = Peq;
    int64_t np_eq = ~p_eq;                      // complement mask
    int64_t HMASK = 1LL << (n - 1);             // test the "leftmost" column
    int64_t VP    = (1LL << n) - 1;             // vertical positive = all 1s
    int64_t VN    = 0;                          // vertical negative = all 0s
    int     score = n;                          // start with distance = n

    for (int j = 0; j < n; ++j) {
        // Bj is always the equality mask
        int64_t Bj = p_eq;

        // Core Myers update
        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);

        // Update score based on the leftmost bit of HP/HN
        if (HP & HMASK)        score += 1;
        else if (HN & HMASK)   score -= 1;

        // Shift for next iteration
        int64_t X  = (HP << 1) | 1;
        VN         = X & D0;
        VP         = (HN << 1) | ~(X | D0);
    }

    // Now score is the edit distance in [0..n]
    return score;
}

// Reversed variant: iterate from left-most position to right (i.e. j goes from n-1 downto 0)
// and use a left mask that corresponds to the most–significant bit.
int bit_ld_rev(int n, int64_t p, int64_t t) {
    int64_t np = ~p;
    // LMASK set to check the leftmost column of the n–bit window.
    int64_t LMASK = (1LL << (n - 1));
    int64_t VP = (1LL << n) - 1;
    int64_t VN = 0;
    int score = n;
    
    for (int j = n - 1; j >= 0; --j) {
        int64_t Bj = ((t & (1LL << j)) != 0) ? p : np;
        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);
        
        if ((HP & LMASK) != 0)
            score += 1;
        else if ((HN & LMASK) != 0)
            score -= 1;
        
        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    }
    return score;
}

int bit_ld_rev_debug(int n, int64_t p, int64_t t) {
    int64_t np = ~p;
    int64_t LMASK = (1LL << (n - 1));  // left-most mask
    int64_t VP = (1LL << n) - 1;
    int64_t VN = 0;
    int score = n;
    
    std::cout << "n = " << n << ", LMASK = " << std::bitset<64>(LMASK) << std::endl;
    std::cout << "Initial VP = " << std::bitset<64>(VP)
              << ", VN = " << std::bitset<64>(VN)
              << ", score = " << score << std::endl;
              
    for (int j = n - 1; j >= 0; --j) {
        int64_t Bj = ((t & (1LL << j)) != 0) ? p : np;
        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);
        
        // Print intermediate state for this column.
        std::cout << "j = " << j << ": " 
                  << "Bj = " << std::bitset<64>(Bj) << ", "
                  << "D0 = " << std::bitset<64>(D0) << ", "
                  << "HP = " << std::bitset<64>(HP) << ", "
                  << "HN = " << std::bitset<64>(HN);
        
        if ((HP & LMASK) != 0) {
            score += 1;
            std::cout << " -> score +1";
        }
        else if ((HN & LMASK) != 0) {
            score -= 1;
            std::cout << " -> score -1";
        }
        std::cout << ", score now = " << score << std::endl;
        
        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
        
        std::cout << "After shift: VP = " << std::bitset<64>(VP)
                  << ", VN = " << std::bitset<64>(VN) << std::endl;
    }
    return score;
}

int bit_ld_window(int n, int64_t p, int64_t t, int offset, int L) {
    // Compute the starting bit (from the left): highest index = n - 1 minus offset.
    int start = 0;
    // The end index for the window:
    int end = 31;
    if (end < 0) {  // Adjust if the window exceeds available bits.
        end = 0;
        L = end - start + 1;
    }
    
    int64_t np = ~p;
    // Use a mask that checks the column corresponding to the end of the window.
    int64_t MASK = (1LL << end);
    // Setup DP vectors for the window length L.
    int64_t VP = (1LL << L) - 1;  // L ones.
    int64_t VN = 0;
    int score = L;  // Initial score equals the window length.
    
    // Process columns from j = start downto j = end.
    for (int j = start; j >= end; --j) {
        int64_t Bj = ((t & (1LL << j)) != 0) ? p : np;
        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);
        
        if ((HP & MASK) != 0)
            score += 1;
        else if ((HN & MASK) != 0)
            score -= 1;
        
        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    }
    return score;
}

int bit_ld_debug(int n, int64_t p, int64_t t) {
    int64_t np = ~p;
    int64_t HMASK = (1LL << (n - 1));        // Mask to check the most-significant bit (leftmost)
    int64_t VP = (1LL << n) - 1;               // VP: n ones (all bits set for positions 0 to n-1)
    int64_t VN = 0;                          // VN: starts at zero
    int score = n;                           // Start the score at n
    
    std::cout << "n = " << n 
              << ", HMASK = " << std::bitset<64>(HMASK) << "\n";
    std::cout << "Initial VP = " << std::bitset<64>(VP)
              << ", VN = " << std::bitset<64>(VN)
              << ", score = " << score << "\n";
    
    for (int j = 0; j < n; ++j) {
        int64_t Bj;
        // Select Bj based on the j-th bit of t (bit j, counting from LSB as 0)
        if ((t & (1LL << j)) != 0)
            Bj = p;
        else
            Bj = np;
        
        int64_t D0 = ((VP + (Bj & VP)) ^ VP) | Bj | VN;
        int64_t HN = VP & D0;
        int64_t HP = VN | ~(VP | D0);
        
        std::cout << "j = " << j << ": ";
        std::cout << "Bj = " << std::bitset<64>(Bj) << ", ";
        std::cout << "D0 = " << std::bitset<64>(D0) << ", ";
        std::cout << "HP = " << std::bitset<64>(HP) << ", ";
        std::cout << "HN = " << std::bitset<64>(HN);
        
        // Update score based on the leftmost bit (i.e., HMASK).
        if ((HP & HMASK) != 0) {
            score += 1;
            std::cout << " -> score +1";
        }
        else if ((HN & HMASK) != 0) {
            score -= 1;
            std::cout << " -> score -1";
        }
        std::cout << ", score now = " << score << "\n";
        
        int64_t X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
        
        std::cout << "    After shift: VP = " << std::bitset<64>(VP)
                  << ", VN = " << std::bitset<64>(VN) << "\n";
    }
    
    return score;
}

// Returns edit distance ≤ max_dist, or max_dist+1 if it exceeds max_dist early.
int bit_lev_wildcard(
    const std::string &P,
    const std::string &T,
    int max_dist = 3
) {
    int m = (int)P.size();
    int n = (int)T.size();
    assert(m <= 64 && "Pattern too long for 64‑bit bit‑parallel algorithm");

    // 1) Build Peq masks once
    uint64_t Peq[128] = {0};
    for (int i = 0; i < m; ++i) {
        char c = P[i];
        uint64_t bit = 1ULL << i;
        if (c == 'N') {
            Peq['A'] |= bit;
            Peq['C'] |= bit;
            Peq['G'] |= bit;
            Peq['T'] |= bit;
        } else {
            Peq[(int)c] |= bit;
        }
    }

    // 2) Initialize Myers state
    uint64_t VP   = ~0ULL;            // all 1s in low m bits
    uint64_t VN   = 0ULL;             // all 0s
    int      dist = m;               // start = m (all deletions)
    uint64_t MSK  = 1ULL << (m - 1);  // mask for the “leftmost” bit

    // 3) Scan the text with early cutoff
    for (int j = 0; j < n; ++j) {
        // lookup Eq mask, treating 'N' in text as wildcard
        uint64_t Eq = (T[j] == 'N')
          ? ((m == 64) ? ~0ULL : ((1ULL << m) - 1))
          : Peq[(int)T[j]];

        // Myers bit‑parallel update
        uint64_t X  = Eq | VN;
        uint64_t D0 = (((X & VP) + VP) ^ VP) | X;
        uint64_t HN = VP & D0;
        uint64_t HP = VN | ~(VP | D0);

        // update distance
        if      (HP & MSK) dist++;
        else if (HN & MSK) dist--;

        // **correct early cutoff**: even if we subtract 1 in every remaining column,
        // can we still get ≤ max_dist?
        int remaining = n - j - 1;
        if (dist - remaining > max_dist) {
            return -1;
        }

        // shift for next character
        uint64_t HPs = (HP << 1) | 1ULL;
        VN = HPs & D0;
        VP = (HN << 1) | ~(HPs | D0);
    }

    // final result
    return (dist > max_dist ? -1 : dist);
}


std::vector<std::pair<int64_seq, std::tuple<int,int,int>>> find_candidates(
    const int64_seq &query, 
    const std::unordered_set<int64_seq> &candidateSet, 
    int max_dist = 4)
{
    std::vector<std::pair<int64_seq, std::tuple<int,int,int>>> results;

    // If the query's bit vector is empty, return an empty result.
    if (query.bits.empty()) return results;
    // Use the query's length and its encoded value (assumed stored in bits[0]).
    int n = query.length;
    int64_t q_bits = query.bits[0];
    // Loop over each candidate.
    for (const auto &candidate : candidateSet) {
        std::string barcode = query.bits_to_sequence();
        //std::cout << "Query barcode: " << barcode << "\n";
        //query.print_int64_seq();
    
        std::string cand_barcode = candidate.bits_to_sequence();
        //std::cout << "Candidate barcode: " << cand_barcode << "\n";
        //candidate.print_int64_seq();
        // Candidate must have been encoded.
        if (candidate.bits.empty()) continue;

        // Optionally only compare candidates of equal length.
        if (candidate.length != n) continue;

        int64_t c_bits = candidate.bits[0];
        //bit_distance
        //int bit_dist = bit_ld(q_bits, c_bits);
        //string distance
        std::string pad_cand_bc = "NN" + cand_barcode + "NN";
        std::string pad_bc = "NN" + barcode + "NN";

        int bit_dist = bit_lev_wildcard(barcode, pad_cand_bc, 3);
        int bit_dist_rev = 0; // ensuring printing
        int bit_dist_window = 100;

        // If any one distance is within the threshold, store this candidate.
        if (bit_dist_rev <= max_dist || bit_dist_window <= max_dist) {
            results.push_back({ candidate, std::make_tuple(bit_dist, 0, 100) });
        }
    }
    return results;
}

//cite flexiplex--looks like the biggest issue is that the short-distance edlib 
uint edit_distance_wildcard(const std::string& s1, const std::string& s2, uint &end, int max_editd){
    std::size_t len1 = s1.size()+1, len2 = s2.size()+1;
    const char * s1_c = s1.c_str(); const char * s2_c = s2.c_str();
    std::vector<uint> dist_holder(len1*len2);
    //initialise the edit distance matrix.
    //penalise for gaps at the start and end of the shorter sequence (j)
    //but not for shifting the start/end of the longer sequence (i,0)
    dist_holder[0]=0; //[0][0]
    for(uint j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
    for(uint i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];
    int best=len2;
    end=len1-1;
    //loop over the distance matrix elements and calculate running distance
    for (uint j = 1; j < len2; ++j) {
        bool any_below_threshold = false; // flag used for early exit
        for (uint i = 1; i < len1; ++i) {
            // Modify the substitution cost to handle wildcards
            int sub = 1; // Default: different bases cost 1
            
            char c1 = s1_c[i - 1];
            char c2 = s2_c[j - 1];
            
            // wildcard match goes here
            if (c1 == c2 || c1 == 'N' || c2 == 'N') {
                sub = 0; // No cost
            }
            
            // If bases match (or wildcards), no cost
            if (sub == 0) {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
            }
            // Otherwise add insertion, deletion or substitution
            else {
                dist_holder[i*len2+j] = std::min({
                    dist_holder[(i-1)*len2+j] + 1,     // deletion
                    dist_holder[i*len2+(j-1)] + 1,     // insertion
                    dist_holder[(i-1)*len2+(j-1)] + 1  // substitution
                });
            }
            
            if (dist_holder[i * len2 + j] <= max_editd) {
                any_below_threshold = true;
            }
            
            // if this is the last row in j
            if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
                // check if this is the best running score
                best = dist_holder[i * len2 + j];
                end = i; // update the end position of alignment
            }
        }
        // early exit to save time.
        if(!any_below_threshold) {
            return(100);
        }
    }
    return best; // return edit distance
}

uint edit_distance_wildcard_2row(const std::string &s1, const std::string &s2, uint &end, int max_editd) {
    // s1: padded barcode; s2: query
    const std::size_t n = s1.size();
    const std::size_t m = s2.size();
    
    std::vector<uint> prev(m + 1), curr(m + 1);
    // initialize first row; note: for insertion penalties in s2, you might need to adjust these
    for (size_t j = 0; j <= m; ++j) {
        prev[j] = j;
    }
    uint best = m;
    end = 0;
    // Process rows: each row corresponds to a character in s1.
    for (size_t i = 1; i <= n; ++i) {
        // gap penalty for s1's position: same as before (0 cost for left edge)
        curr[0] = 0; 
        bool anyBelowThreshold = false;
        for (size_t j = 1; j <= m; ++j) {
            int sub_cost = 1;
            char c1 = s1[i - 1];
            char c2 = s2[j - 1];
            if(c1 == c2 || c1 == 'N' || c2 == 'N') {
                sub_cost = 0;
            }
            // Calculate using deletion (from previous row), insertion (current row, previous column),
            // and substitution.
            curr[j] = std::min({ prev[j] + 1,          // deletion
                                 curr[j - 1] + 1,      // insertion
                                 prev[j - 1] + sub_cost // substitution
                               });
            if (curr[j] <= (uint)max_editd) {
                anyBelowThreshold = true;
            }
            // If this is the last column, update the best score and record the alignment end position.
            if (j == m && curr[j] < best) {
                best = curr[j];
                end = i;
            }
        }
        
        // Early exit: if no cell in the current row is within threshold, there’s no need to continue.
        if (!anyBelowThreshold) {
            return 100;
        }
        prev.swap(curr);
    }
    return best;
}

// Test function for find_candidates
void test_find_candidates() {
    // Define a test barcode (the query)
    std::string test_barcode = "AACAGTTATTACTTCT";
    int64_seq query(test_barcode);
    

    // Define candidate barcode strings.
    std::vector<std::string> candidateStrings = {
        "AACAGTTATTACTTCT",   // exact match
        "AACAGTTATTACTTCA",   // one substitution at end
        "AACAGCTATTACTTCT",   // one substitution in middle
        "AACAGTTATTACTTCG",   // substitution at last base
        "TTCAGTTATTACTTCT",   // substitution at beginning
        "AACAGTTATTGCTTCT",   // one substitution inside
        "AACAGTTATTACTTAT",   // substitution at last base
        "AACAGTTATTACTTCC",   // substitution at end
        "AACAGTTATTACTTGT",   // substitution at end
        "AACAGTTATTACTTC",    // deletion: missing last base
        "AACAGTTATTACTTCTA",   // insertion: extra base at end
        "TTAGGTTATTACTTCT",    // multiple errors
        "AAAACAGTTATTACTTCT"  // two extra bases at beginning
    };

    // Expected best edit distances for each candidate.
    std::vector<int> candidate_dist = {
        0,  // "AACAGTTATTACTTCT" exact match: 0 edits.
        1,  // "AACAGTTATTACTTCA" expected edit distance 1.
        1,  // "AACAGCTATTACTTCT" expected edit distance 1.
        1,  // "AACAGTTATTACTTCG" expected edit distance 1.
        2,  // "TTCAGTTATTACTTCT" expected edit distance 2.
        1,  // "AACAGTTATTGCTTCT" expected edit distance 1.
        1,  // "AACAGTTATTACTTAT" expected edit distance 1.
        1,  // "AACAGTTATTACTTCC" expected edit distance 1.
        1,  // "AACAGTTATTACTTGT" expected edit distance 1.
        1,  // "AACAGTTATTACTTC"  expected edit distance 1 (due to deletion).
        1,   // "AACAGTTATTACTTCTA" expected edit distance 1 (due to insertion).
        4,   // "TTAGGTTATTACTTCT" expected edit distance 4.
        2    // "AAAACAGTTATTACTTCT" edit distance 2
    };
    
    // Build an unordered_set of int64_seq from candidateStrings.
    std::unordered_set<int64_seq> candidateSet;
    for (const auto &s : candidateStrings) {
        int64_seq candidate(s);
        if (candidate.length == query.length) // enforce same length
            candidateSet.insert(candidate);
    }
    
    int max_dist = 3;  // Set maximum edit distance threshold.
    auto candidates = find_candidates(query, candidateSet, max_dist);
    
    // Print the test results.
    std::cout << "\n=== Test find_candidates() ===" << std::endl;
    std::cout << "Query barcode: " << test_barcode << std::endl;
    std::cout << "Candidates within edit distance " << max_dist << ":" << std::endl;
    
    // Instead of using pointer arithmetic, use a standard index-based loop
    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto &pair = candidates[i];
        auto [d, d_rev, d_win] = pair.second;
        // Compute the "best" (lowest) distance among the three.
        int best = std::min({d, d_rev});
        std::cout << "Candidate: " << pair.first.bits_to_sequence() 
                  << ", Distances: (" << d << ", " << d_rev << ")"
                  << ", Best: " << best << std::endl;
    }
}

void test_find_candidates_batch(
    const std::unordered_map<int64_t, std::string>& queryMap,
    const std::unordered_map<int64_t, std::string>& candidateMap, 
    int max_dist = 3)
{
    // Build the candidate set from candidateMap.
    std::unordered_set<int64_seq> candidateSet;
    for (const auto& pair : candidateMap) {
        int64_seq candidate(pair.second);
        if (!candidate.bits.empty())
            candidateSet.insert(candidate);
    }
    
    // Container for results for each query.
    std::vector<std::vector<std::pair<std::string, std::tuple<int,int,int>>>> all_results;
    all_results.reserve(queryMap.size());
    
    std::cout << "\n=== Batch find_candidates() Test (using unordered_map for queries and candidates) ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Process each query in the queryMap.
    for (const auto &qPair : queryMap) {
        int64_seq query(qPair.second);
        auto candidate_matches = find_candidates(query, candidateSet, 3);
        
        std::vector<std::pair<std::string, std::tuple<int,int,int>>> query_results;
        for (const auto &match : candidate_matches) {
            query_results.push_back({ match.first.bits_to_sequence(), match.second });
        }
        // Sort results by the “best” distance (lowest among the three in the tuple).
        std::sort(query_results.begin(), query_results.end(), 
                  [](const auto &a, const auto &b) {
                      int best_a = std::min({ std::get<0>(a.second), std::get<1>(a.second), std::get<2>(a.second) });
                      int best_b = std::min({ std::get<0>(b.second), std::get<1>(b.second), std::get<2>(b.second) });
                      return best_a < best_b;
                  });
        
        all_results.push_back(query_results);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Print summary.
    size_t total_matches = 0;
    for (const auto& results : all_results) {
        total_matches += results.size();
    }
    std::cout << "Processed " << queryMap.size() << " queries against candidate data." << std::endl;
    std::cout << "Found a total of " << total_matches 
              << " candidate matches with edit distance <= " << 3 
              << " in " << elapsed << " ms." << std::endl;
    std::cout << "Average time per query: " << elapsed / queryMap.size() << " ms." << std::endl;
    
    // For each query, print only the best candidates.
    const size_t max_to_print = 5;
    size_t query_index = 0;
    for (const auto &qPair : queryMap) {
        const auto &results = all_results[query_index];
        if (results.empty()) {
            std::cout  << "";
        } else {
            // Use the best (lowest) edit distance among the first candidate's tuple.
            std::cout << "\nQuery (" << query_index << "): " << qPair.second << std::endl;

            int best_distance = std::min({ std::get<0>(results.front().second), 
                                           std::get<1>(results.front().second), 
                                           std::get<2>(results.front().second) });
            // Filter the best candidates that achieve this distance.
            std::vector<std::pair<std::string, std::tuple<int,int,int>>> best_candidates;
            std::copy_if(results.begin(), results.end(), std::back_inserter(best_candidates),
                         [best_distance](const std::pair<std::string, std::tuple<int,int,int>> &r) {
                             int best_val = std::min({ std::get<0>(r.second), std::get<1>(r.second), std::get<2>(r.second) });
                             return best_val == best_distance;
                         });
            std::cout << "Best edit distance = " << best_distance << "; ";
            if (best_candidates.size() > max_to_print) {
                std::cout << "showing first " << max_to_print << " of " << best_candidates.size() << " candidates:" << std::endl;
                for (size_t i = 0; i < max_to_print; ++i) {
                    auto [d, d_rev, d_win] = best_candidates[i].second;
                    std::cout << "  Candidate: " << best_candidates[i].first 
                              << ", Distances: (" << d << ", " << d_rev << ", " << d_win << ")" << std::endl;
                }
                std::cout << "  ... (" << best_candidates.size() - max_to_print << " more with edit distance " 
                          << best_distance << ")" << std::endl;
            } else {
                for (const auto &cand : best_candidates) {
                    auto [d, d_rev, d_win] = cand.second;
                    std::cout << "  Candidate: " << cand.first 
                              << ", Distances: (" << d << ", " << d_rev << ", " << d_win << ")" << std::endl;
                }
            }
        }
        ++query_index;
    }
}

void levenshtein_compare(const std::string& query, const std::vector<std::string>& barcodes, int max_dist = 4){
    std::cout << "\n=== Levenshtein Edit Distance Comparison for '" << query << "' ===" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string, int>> results;
    std::vector<std::string> padded_barcodes;
    for (const auto& barcode : barcodes) {
        padded_barcodes.push_back("NNN" + barcode + "NNN");
    }
    for (size_t i = 0; i < padded_barcodes.size(); ++i) {
        uint end = 22;
        int dist = edit_distance_wildcard(padded_barcodes[i], query, end, 4);
        if (dist <= max_dist) {
            // Push the original barcode instead of the padded version.
            results.push_back({barcodes[i], dist});
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Sort by edit distance
    std::sort(results.begin(), results.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Group by edit distance
    std::unordered_map<int, std::vector<std::string>> grouped;
    for (const auto& [barcode, dist] : results) {
        grouped[dist].push_back(barcode);
    }

    std::cout << "Found " << results.size() << " matches with edit distance <= " << max_dist << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    
}

void levenshtein_compare_opt(const std::string& query, const std::vector<std::string>& barcodes, int max_dist = 4){
    std::cout << "\n=== Levenshtein Edit Distance Comparison, TWO ROW for '" << query << "' ===" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<std::string, int>> results;
    std::vector<std::string> padded_barcodes;
    for (const auto& barcode : barcodes) {
        padded_barcodes.push_back("NNN" + barcode + "NNN");
    }
    for (const auto& barcode : padded_barcodes) {
        uint end = 22;
        int dist = edit_distance_wildcard_2row(barcode, query, end, 4);
        if (dist <= max_dist) {
            results.push_back({barcode, dist});
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Sort by edit distance
    std::sort(results.begin(), results.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Group by edit distance
    std::unordered_map<int, std::vector<std::string>> grouped;
    for (const auto& [barcode, dist] : results) {
        grouped[dist].push_back(barcode);
    }

    std::cout << "Found " << results.size() << " matches with edit distance <= " << max_dist << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    
}

//function that returns the unique best match if available.
std::string levenshtein_compare_sort_(const std::string& query, const std::vector<std::string>& barcodes, int max_dist = 4) {
    std::cout << "\n=== Levenshtein Edit Distance Comparison for '" << query << "' ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    // Store each candidate barcode (original, unpadded) with its computed edit distance.
    std::vector<std::pair<std::string, int>> results;
    std::vector<std::string> padded_barcodes;
    for (const auto& barcode : barcodes) {
        padded_barcodes.push_back("NNN" + barcode + "NNN");
    }
    uint end = 22;
    // Compute distances.
    for (size_t i = 0; i < padded_barcodes.size(); ++i) {
        int dist = edit_distance_wildcard(padded_barcodes[i], query, end, 4);
        if (dist <= max_dist) {
            results.push_back({barcodes[i], dist});
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    double distance_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    // Sort the results by edit distance (lowest first).
    std::sort(results.begin(), results.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Group the results by edit distance.
    std::unordered_map<int, std::vector<std::string>> grouped;
    for (const auto& [barcode, dist] : results) {
        grouped[dist].push_back(barcode);
    }
    
    // Determine the unique best result:
    // Find the smallest edit distance and check if its group has exactly one candidate.
    auto search_time = std::chrono::high_resolution_clock::now();
    int best_distance = -1;
    std::string best_result = "";
    for (const auto& kv : grouped) {
        if (best_distance == -1 || kv.first < best_distance) {
            best_distance = kv.first;
        }
    }
    if (best_distance != -1 && grouped[best_distance].size() == 1) {
        best_result = grouped[best_distance][0];
        std::cout << "Unique best result (edit distance " << best_distance << "): " << best_result << std::endl;
    } else if (best_distance != -1) {
        std::cout << "Ambiguous best results at edit distance " << best_distance 
                  << " (" << grouped[best_distance].size() << " candidates)." << std::endl;
    } else {
        std::cout << "No candidate barcodes found within the maximum edit distance." << std::endl;
    }

    auto end_find = std::chrono::high_resolution_clock::now();
    double find_time = std::chrono::duration<double, std::milli>(end_find - search_time).count();

    std::cout << "Found " << results.size() << " matches with edit distance <= " << max_dist << std::endl;
    std::cout << "Time taken for levenshtein distances: " << distance_time << " ms" << std::endl;
    std::cout << "Time taken for finding best match: " << find_time << " ms" << std::endl;
    std::cout << "Total run time: " << distance_time + find_time << " ms" << std::endl;

    return best_result;
}

void levenshtein_compare_batch(
    const std::vector<std::string>& queries, 
    const std::vector<std::string>& barcodes, 
    int max_dist = 4) {
    
    // Open the output CSV file
    std::ofstream output_file("/Users/cmv/Desktop/lv_output.csv");
    if(!output_file.is_open()) {
        std::cerr << "Error: Could not open output file for writing." << std::endl;
        return;
    }
    // Write CSV header
    output_file << "query,matched,edit_dist\n";
    
    // Store results for all queries: query_index -> [(barcode, distance)]
    std::vector<std::vector<std::pair<std::string, int>>> all_results(queries.size());
    size_t total_matches = 0;
    
    std::cout << "\n=== Batch Levenshtein Edit Distance Comparison for " << queries.size() << " queries ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::string> padded_barcodes;

    for (const auto& barcode : barcodes) {
        padded_barcodes.push_back(barcode);
    }

    for (size_t q = 0; q < queries.size(); ++q) {
        const std::string& query = queries[q];
        auto& results = all_results[q];
        for (size_t i = 0; i < barcodes.size(); ++i) {
            uint end = 16;
            int dist = edit_distance_wildcard(padded_barcodes[i], query, end, max_dist);
            if (dist <= max_dist) {
                results.push_back({barcodes[i], dist});
            }
        }

        std::sort(results.begin(), results.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
    }    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    // Write all results to the CSV file
    for (size_t q = 0; q < queries.size(); ++q) {
        const std::string& query = queries[q];
        const auto& results = all_results[q];
        
        total_matches += results.size();
        
        for (const auto& [barcode, dist] : results) {
            {
                output_file << query << "," << barcode << "," << dist << "\n";
            }
        }
    }
    
    output_file.close();
    
    std::cout << "Found " << total_matches << " total matches with edit distance <= " << max_dist << std::endl;
    std::cout << "Time taken: " << elapsed << " ms" << std::endl;
    std::cout << "Average time per query: " << elapsed / queries.size() << " ms" << std::endl;
    std::cout << "Results written to: " << "/Users/cmv/Desktop/lv_output.csv" << std::endl;
}

// Test file for barcode_correction functionalities
int main() {
    std::string barcode_file = "/Users/cmv/Desktop/la1_barcodes.csv";
    std::string barcode_to_correct_file = "/Users/cmv/Desktop/barcodes_to_correct.csv";
    // Load barcodes
    std::cout << "Loading barcodes from " << barcode_file << "..." << std::endl;
    auto barcodes = load_barcodes(barcode_file);
    auto putative_barcodes = load_barcodes(barcode_to_correct_file);
    std::cout << "Loaded " << putative_barcodes.size() << " barcodes to correct." << std::endl;
    auto bit_barcodes = load_bit_barcodes(barcode_file);
    std::cout << "Loaded " << barcodes.size() << " barcodes." << std::endl;
    std::cout << "Loaded " << bit_barcodes.size() << " barcodes and converted to bits." << std::endl;
    auto bit_barcodes_to_correct = load_bit_barcodes(barcode_to_correct_file);
    std::cout << "Loaded " << bit_barcodes_to_correct.size() << " barcodes to correct." << std::endl;
    if (barcodes.empty()) {
        std::cerr << "No barcodes loaded. Exiting." << std::endl;
        return 1;
    }    

    // --- Test 1: Encode/Decode a 16-base read ---
    std::string read = "AACAGTTATTACTTCT"; // 16 bases
    int64_seq encoded(read);
    std::string decoded = encoded.bits_to_sequence();
    
    std::cout << "Original read: " << read << "\n";
    std::cout << "Decoded read:  " << decoded << "\n";
    std::cout << "Encoded representation:\n";
    encoded.print_int64_seq();

    // --- Test 2: Barcode Mutation ---
    // Use a poly-A barcode of 16 bases and generate mutations with 2 rounds.
    std::string barcode = "AACAGTTATTACTTCT"; // 16 bases (all A's)
    int64_seq int64_barcode(barcode);
    
    std::cout << "\nSequence: " << barcode << "\n";
    int64_barcode.print_int64_seq();
    
    // Generate recursive mutations (2 rounds) for the first chunk of the poly-A barcode.
    // (Assuming your poly-A fits in one int64_t, since 16 < 32)

    auto start_mutations = std::chrono::high_resolution_clock::now();
    int rounds = 2;
    std::vector<int64_t> recursive_mutations = generate_mutated_barcodes(int64_barcode, rounds);
    auto stop_mutations = std::chrono::high_resolution_clock::now();
    double time_mutations = std::chrono::duration<double, std::milli>(stop_mutations - start_mutations).count();

    std::cout << "\n============= Recursive mutations for " << barcode << ": =============\n";
    std::cout << "\nGenerated " << recursive_mutations.size() << " barcodes after " << rounds << " rounds of mutation.";
    std::cout << "\nTime elapsed for generating mutations: " << time_mutations << " milliseconds.\n";
    auto seconds_per_million = time_mutations*1000;
    auto minutes_per_million = seconds_per_million/60;
    std::cout << "Assuming that this time holds, generating mutations for a dataset of ~1M unique barcodes single-threaded will take: " << 
    //time mutations by 1000 is 1 million barcodes divided by 1000 milliseconds to seconds conversion
    seconds_per_million <<
     " seconds (" <<
     minutes_per_million << " minutes)\n";

    auto seconds_per_barcode = time_mutations/1000;
    auto seconds_per_rageseq = seconds_per_barcode*7092;
    std::cout << "For a much more reasonable estimate, generating mutations for a dataset of 7092 barcodes single-threaded will take: " << seconds_per_rageseq << " seconds.\n\n";

    auto start_map_scan = std::chrono::high_resolution_clock::now();
    for (auto m : recursive_mutations) {
        int64_seq temp;
        temp.length = 16;
        temp.bits.push_back(m);
        auto it = bit_barcodes.find(m);
        if (it != bit_barcodes.end()) {
            temp.print_int64_seq();
        }
    }
    auto stop_map_scan = std::chrono::high_resolution_clock::now();
    double total_map_scan = std::chrono::duration<double, std::milli>(stop_map_scan - start_map_scan).count();
    std::cout << "Time elapsed for scanning " << recursive_mutations.size() << " sequences: " << total_map_scan << " milliseconds.\n";
    std::cout << "Average time scanning per sequence: " << total_map_scan/recursive_mutations.size() << " ms. \n";
    
    // --- Test 3: Barcode Shifting ---
    // Generate shifted barcode candidates for a selected barcode with shift of (?)
    //int shift = 6;
    auto start_shifting = std::chrono::high_resolution_clock::now();
    //std::vector<int64_t> shifted_candidates = generate_shifted_barcodes(int64_barcode, shift);
    int shift = 4;
    auto shifted_candidates = generate_int64_shifted_barcodes(int64_barcode, shift);
    
    auto stop_shifting = std::chrono::high_resolution_clock::now();
    double time_shifting = std::chrono::duration<double, std::milli>(stop_shifting - start_shifting).count();

    std::cout << "\n============= Shifting mutations for " << barcode << ": =============\n";
    std::cout << "\nGenerated " << shifted_candidates.size() << " barcodes via " << shift << " shifted bases.";
    std::cout << "\nTime elapsed for generating shifted_candidates: " << time_shifting << " milliseconds.\n\n";

    auto start_shift_scan = std::chrono::high_resolution_clock::now();
    for (auto s : shifted_candidates) {
        std::string shifted_barcode = s.bits_to_sequence();
        std::cout << "Shifted candidate: " <<  shifted_barcode << "\n";
    }
    auto stop_shift_scan = std::chrono::high_resolution_clock::now();
    double total_shift_scan = std::chrono::duration<double, std::milli>(stop_shift_scan - start_shift_scan).count();
    std::cout << "Time elapsed for scanning " << shifted_candidates.size() << " sequences: " << total_shift_scan << " milliseconds.\n";
    std::cout << "Average time scanning per sequence: " << total_shift_scan/shifted_candidates.size() << " ms. \n";

    // Parse optional parameters
    int max_edit_distance = 3;
    double min_similarity = 70.0;
    // Run comparisons
    edlib_compare(barcode, barcodes, 4);
   // rapidfuzz_compare(barcode, barcodes, min_similarity);
    levenshtein_compare(barcode, barcodes, 4);
    //levenshtein_compare_sort_(barcode, barcodes, 4);
    test_find_candidates();
    //test_find_candidates_batch(bit_barcodes_to_correct, bit_barcodes, 4);
    //levenshtein_compare_opt(barcode, barcodes, 4);
    //levenshtein_compare_batch(putative_barcodes, barcodes,  4);
    //edlib_compare_batch(putative_barcodes, barcodes, 4, false);
    return 0;
}
