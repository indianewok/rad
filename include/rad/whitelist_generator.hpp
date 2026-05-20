#pragma once
#include "rad_headers.h"

struct extracted_bc {
    std::string read_id;
    std::string sequence;
    int position;
    bool is_rc;
    std::string barcode_1;
    std::string barcode_2;
};

static const EdlibEqualityPair kWildcardEqualities[8] = {
                {'N', 'A'}, 
                {'N', 'C'}, 
                {'N', 'G'}, 
                {'N', 'T'},
                {'A', 'N'}, 
                {'C', 'N'}, 
                {'G', 'N'}, 
                {'T', 'N'}
};

// Compute allowed edit distance for adapter search given a ratio.
// - Ratio <= 0 -> exact match only (0 edits).
// - Positive ratios are rounded up so small adapters still allow >=1 edit.
// - Clamp to int range for safety.
static int compute_max_edit_distance(size_t adapter_len, double max_edit_distance_ratio) {
    if (adapter_len == 0 || !std::isfinite(max_edit_distance_ratio) || max_edit_distance_ratio <= 0.0) {
        return 0;
    }
    double raw = std::floor(static_cast<double>(adapter_len) * max_edit_distance_ratio);
    raw = std::min(raw, static_cast<double>(std::numeric_limits<int>::max()));
    int max_edits = static_cast<int>(raw);
    return std::max(1, max_edits);
}

// Extract barcode sequence given adapter_seq position
std::string extract_barcode(
    const std::string& sequence, int adapter_seq_pos, int adapter_seq_len, int bc_length, int m_left, int m_right, bool is_rc) {
    int start_pos, end_pos; 
    std::string final_bc;
    if (is_rc) {
        // For reverse complement, extract before the adapter_seq
        start_pos = std::max(0, adapter_seq_pos - m_right - bc_length);
        end_pos = std::min((int)sequence.length(), adapter_seq_pos - m_right);
    } else {
        // For forward, extract after the adapter_seq
        start_pos = std::max(0, adapter_seq_pos + adapter_seq_len + m_left);
        end_pos = std::min((int)sequence.length(), start_pos + bc_length + m_right);
    }
    if (start_pos >= end_pos || start_pos < 0 || end_pos > (int)sequence.length()) {
        return "";
    }
    if(is_rc){
        std::string bc_seq = sequence.substr(start_pos, end_pos - start_pos);
        final_bc = seq_utils::revcomp(bc_seq);
    } else {
        final_bc = sequence.substr(start_pos, end_pos - start_pos);
    }
    return final_bc;
}

// Structure to hold a read for processing
struct read_chunk {
    std::string read_id;
    std::string sequence;
};

// Process a chunk of reads in parallel
std::vector<extracted_bc> process_read_chunk(const std::vector<read_chunk>& chunk,
                                               const std::string& adapter_seq,
                                               const std::string& adapter_seq_rc,
                                               int bc_length, int m_left, int m_right,
                                               double max_edit_distance_ratio) {
    std::vector<extracted_bc> chunk_results;

    const int max_edits_forward = compute_max_edit_distance(adapter_seq.length(), max_edit_distance_ratio);
    const int max_edits_rc = compute_max_edit_distance(adapter_seq_rc.length(), max_edit_distance_ratio);
    const bool exact_match_only = (max_edits_forward == 0 && max_edits_rc == 0);

    EdlibAlignConfig forward_config{};
    EdlibAlignConfig rc_config{};
    if (!exact_match_only) {
        forward_config = edlibNewAlignConfig(
            max_edits_forward,
            EDLIB_MODE_HW,
            EDLIB_TASK_LOC,
            kWildcardEqualities,
            8);
        rc_config = edlibNewAlignConfig(
            max_edits_rc,
            EDLIB_MODE_HW,
            EDLIB_TASK_LOC,
            kWildcardEqualities,
            8);
    }

    // Reserve space to avoid reallocations
    chunk_results.reserve(chunk.size() / 5);
    #pragma omp parallel
    {
        std::vector<extracted_bc> thread_results;
        thread_results.reserve(chunk.size() / (10 * omp_get_num_threads()));
        
        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < chunk.size(); ++i) {
            const auto& read = chunk[i];

            if (exact_match_only) {
                bool forward_found = false;
                auto pos = read.sequence.find(adapter_seq);
                if (pos != std::string::npos && pos <= static_cast<size_t>(std::numeric_limits<int>::max())) {
                    forward_found = true;
                    int adapter_seq_pos = static_cast<int>(pos);
                    std::string barcode = extract_barcode(
                        read.sequence,
                        adapter_seq_pos,
                        adapter_seq.length(),
                        bc_length,
                        m_left,
                        m_right,
                        false);

                    if (!barcode.empty()) {
                        thread_results.push_back({read.read_id, barcode, adapter_seq_pos, false, "", ""});
                    }
                }

                if (!forward_found) {
                    auto rc_pos = read.sequence.find(adapter_seq_rc);
                    if (rc_pos != std::string::npos && rc_pos <= static_cast<size_t>(std::numeric_limits<int>::max())) {
                        int adapter_seq_pos = static_cast<int>(rc_pos);
                        std::string barcode = extract_barcode(
                            read.sequence,
                        adapter_seq_pos,
                        adapter_seq_rc.length(),
                        bc_length,
                        m_left,
                        m_right,
                            true);
                        if (!barcode.empty()) {
                            thread_results.push_back({read.read_id, barcode, adapter_seq_pos, true, "", ""});
                        }
                    }
                }
                continue;
            }

            // Search for adapter_seq (forward)
            bool forward_found = false;
            EdlibAlignResult result = edlibAlign(
                adapter_seq.c_str(), adapter_seq.length(),
                read.sequence.c_str(), read.sequence.length(),
                forward_config);
            
            if (result.status == EDLIB_STATUS_OK && result.numLocations > 0) {
                forward_found = true;
                // Found adapter_seq - extract barcode
                int adapter_seq_pos = result.startLocations[0];
                std::string barcode = extract_barcode(read.sequence, adapter_seq_pos, adapter_seq.length(), 
                                                    bc_length, m_left, m_right, false);
                
                if (!barcode.empty()) {
                    thread_results.push_back({read.read_id, barcode, adapter_seq_pos, false, "", ""});
                }
            }
            edlibFreeAlignResult(result);

            if (forward_found) {
                continue;
            }
            
            // Search for adapter_seq (reverse complement)
            result = edlibAlign(adapter_seq_rc.c_str(), adapter_seq_rc.length(),
                              read.sequence.c_str(), read.sequence.length(),
                              rc_config);
            
            if (result.status == EDLIB_STATUS_OK && result.numLocations > 0) {
                // Found reverse complement adapter_seq - extract barcode
                int adapter_seq_pos = result.startLocations[0];
                std::string barcode = extract_barcode(read.sequence, adapter_seq_pos, adapter_seq_rc.length(), 
                                                    bc_length, m_left, m_right, true);
                if (!barcode.empty()) {
                    thread_results.push_back({read.read_id, barcode, adapter_seq_pos, true, "", ""});
                }
            }
            edlibFreeAlignResult(result);
        }
        
        // Merge thread results
        #pragma omp critical
        {
            chunk_results.insert(chunk_results.end(), thread_results.begin(), thread_results.end());
        }
    }
    return chunk_results;
}

// Parse FASTQ file (regular or gzipped) and extract barcodes
std::vector<extracted_bc> process_fastq(const std::string& input_path, 
                                          const std::string& adapter_seq,
                                          int bc_length, 
                                          int m_left, 
                                          int m_right,
                                          int max_reads, 
                                          double max_edit_distance_ratio = 0.1,
                                          int chunk_size = 10000, 
                                          int num_threads = 0) {
    std::vector<extracted_bc> results;

    if (!std::isfinite(max_edit_distance_ratio) || max_edit_distance_ratio < 0.0) {
        std::cerr << "Warning: invalid max_edit_distance_ratio (" << max_edit_distance_ratio
                  << "); clamping to 0 (exact matches only)\n";
        max_edit_distance_ratio = 0.0;
    } else if (max_edit_distance_ratio == 0.0) {
        std::cout << "Max edit distance ratio set to 0 - requiring exact adapter matches\n";
    }
    
    // Set number of OpenMP threads
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    
    std::cout << "Using " << omp_get_max_threads() << " threads for parallel processing" << std::endl;

    std::unique_ptr<file_streaming> files_ptr;
    std::unique_ptr<read_streaming> reader_ptr;

    try {
        files_ptr = std::make_unique<file_streaming>(input_path, std::max(1, num_threads));
        reader_ptr = std::make_unique<read_streaming>(*files_ptr);
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not open input " << input_path << ": " << e.what() << "\n";
        return results;
    }
    
    std::string adapter_seq_rc = seq_utils::revcomp(adapter_seq);
    int reads_processed = 0;
    int chunks_processed = 0;
    
    while (max_reads <= 0 || reads_processed < max_reads) {
        // Read a chunk of data
        std::vector<read_chunk> chunk;
        chunk.reserve(chunk_size);
        
        for (int i = 0; i < chunk_size; ++i) {
            if (max_reads > 0 && reads_processed >= max_reads) {
                break;
            }
            auto rec = reader_ptr->next_sequence();
            if (!rec) {
                break;
            }
            reads_processed++;
            chunk.push_back({rec->id, rec->seq});

        }
        
        if (chunk.empty()) break;
        
        chunks_processed++;
        if (chunks_processed % 100 == 0) {
            std::cout << "Processed " << reads_processed << " reads in " << chunks_processed << " chunks, found " << results.size() 
                      << " barcodes so far" << std::endl;
        }
        
        // Process chunk in parallel
        auto chunk_results = process_read_chunk(
            chunk, 
            adapter_seq, 
            adapter_seq_rc, 
            bc_length, 
            m_left, 
            m_right, 
            max_edit_distance_ratio);
        
        // Merge results
        results.insert(results.end(), chunk_results.begin(), chunk_results.end());
    };
    
    std::cout << "Processing complete: " << reads_processed 
              << " reads processed, " << results.size() << " barcodes extracted" << std::endl;
    
    return results;
}

// -----------------------------------------------------------------------
// Two-part barcode (BC1 + BC2) support
// -----------------------------------------------------------------------

// Load a two-part barcode whitelist CSV. The first row is a header; subsequent
// rows are bare barcode sequences.  Returns {sequence set, sorted vector
// of unique lengths found in the whitelist}.
struct split_whitelist {
    std::unordered_set<std::string> seqs;
    std::vector<int> lengths;
    std::vector<std::string> ordered;              // insertion-order list for index mapping
    std::unordered_map<std::string, uint16_t> idx; // sequence -> index
};

split_whitelist load_split_whitelist_csv(const std::string& path) {
    split_whitelist wl;

    std::vector<std::string> lines;
    try {
        lines = streaming_utils::import_text(path);
    } catch (const std::exception& e) {
        std::cerr << "Error: Could not open split-barcode whitelist: " << path
                  << " (" << e.what() << ")\n";
        return wl;
    }

    bool first_line = true;
    for (auto line : lines) {
        line = seq_utils::trim(line);
        if (line.empty()) continue;
        if (first_line) {
            first_line = false;
            if (seq_utils::trim(line) == "whitelist_bcs") continue; // header
        }
        auto comma = line.find(',');
        if (comma != std::string::npos) {
            line = seq_utils::trim(line.substr(0, comma));
        }
        if (line.empty()) continue;
        if (wl.seqs.insert(line).second) {
            wl.idx[line] = static_cast<uint16_t>(wl.ordered.size());
            wl.ordered.push_back(line);
        }
        wl.lengths.push_back((int)line.size());
    }

    std::sort(wl.lengths.begin(), wl.lengths.end());
    wl.lengths.erase(std::unique(wl.lengths.begin(), wl.lengths.end()), wl.lengths.end());
    return wl;
}

// Search for a valid BC1+BC2 pair in a read window that begins at
// adapter_end + umi_length.  The extra search offset range
// [search_extra_min, search_extra_max] handles the small positional
// jitter expected in two-part barcode protocols.
// Returns split barcode fields when a valid pair is found.
struct split_barcode_hit {
    std::string bc1;
    std::string bc2;
    std::string combined() const { return bc1 + bc2; }
};

std::optional<split_barcode_hit> extract_split_barcode(
    const std::string& seq,
    int adapter_end,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max)
{
    const int seq_len = (int)seq.size();
    const int window_base = adapter_end + umi_length;

    for (int extra = search_extra_min; extra <= search_extra_max; ++extra) {
        const int bc1_start = window_base + extra;
        if (bc1_start < 0 || bc1_start >= seq_len) continue;

        for (int len1 : bc1_lengths) {
            const int bc1_end = bc1_start + len1;
            if (bc1_end > seq_len) continue;

            if (bc1_set.count(seq.substr(bc1_start, len1))) {
                for (int len2 : bc2_lengths) {
                    if (bc1_end + len2 > seq_len) continue;
                    std::string bc2 = seq.substr(bc1_end, len2);
                    if (bc2_set.count(bc2)) {
                        split_barcode_hit hit;
                        hit.bc1 = seq.substr(bc1_start, len1);
                        hit.bc2 = std::move(bc2);
                        return hit;
                    }
                }
            }
        }
    }
    return std::nullopt;
}

std::vector<extracted_bc> process_read_chunk_split_barcode(
    const std::vector<read_chunk>& chunk,
    const std::string& adapter_seq,
    const std::string& adapter_seq_rc,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max,
    double max_edit_distance_ratio)
{
    std::vector<extracted_bc> chunk_results;
    chunk_results.reserve(chunk.size() / 5);

    const int max_edits_fwd = compute_max_edit_distance(adapter_seq.size(), max_edit_distance_ratio);
    const int max_edits_rc  = compute_max_edit_distance(adapter_seq_rc.size(), max_edit_distance_ratio);
    const bool exact_only   = (max_edits_fwd == 0 && max_edits_rc == 0);

    EdlibAlignConfig fwd_cfg{}, rc_cfg{};
    if (!exact_only) {
        fwd_cfg = edlibNewAlignConfig(max_edits_fwd, EDLIB_MODE_HW, EDLIB_TASK_LOC, kWildcardEqualities, 8);
        rc_cfg  = edlibNewAlignConfig(max_edits_rc,  EDLIB_MODE_HW, EDLIB_TASK_LOC, kWildcardEqualities, 8);
    }

    #pragma omp parallel
    {
        std::vector<extracted_bc> thread_results;
        thread_results.reserve(chunk.size() / (10 * omp_get_num_threads()));

        #pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < chunk.size(); ++i) {
            const auto& read = chunk[i];

            auto try_pair = [&](const std::string& seq, int adapter_end, bool is_rc) -> bool {
                auto hit = extract_split_barcode(
                    seq, adapter_end, umi_length,
                    bc1_set, bc1_lengths, bc2_set, bc2_lengths,
                    search_extra_min, search_extra_max);
                if (hit.has_value()) {
                    thread_results.push_back({
                        read.read_id,
                        hit->combined(),
                        adapter_end,
                        is_rc,
                        hit->bc1,
                        hit->bc2
                    });
                    return true;
                }
                return false;
            };

            if (exact_only) {
                auto pos = read.sequence.find(adapter_seq);
                if (pos != std::string::npos) {
                    try_pair(read.sequence, (int)(pos + adapter_seq.size()), false);
                } else {
                    pos = read.sequence.find(adapter_seq_rc);
                    if (pos != std::string::npos) {
                        std::string rc_seq = seq_utils::revcomp(read.sequence);
                        try_pair(rc_seq, (int)(rc_seq.size() - pos), true);
                    }
                }
                continue;
            }

            // Fuzzy forward
            bool found = false;
            EdlibAlignResult res = edlibAlign(
                adapter_seq.c_str(), adapter_seq.size(),
                read.sequence.c_str(), read.sequence.size(), fwd_cfg);
            if (res.status == EDLIB_STATUS_OK && res.numLocations > 0) {
                found = try_pair(read.sequence, res.endLocations[0] + 1, false);
            }
            edlibFreeAlignResult(res);

            if (!found) {
                res = edlibAlign(
                    adapter_seq_rc.c_str(), adapter_seq_rc.size(),
                    read.sequence.c_str(), read.sequence.size(), rc_cfg);
                if (res.status == EDLIB_STATUS_OK && res.numLocations > 0) {
                    std::string rc_seq = seq_utils::revcomp(read.sequence);
                    // In the RC'd read the adapter ends at (read_len - start_loc)
                    try_pair(rc_seq, (int)(rc_seq.size() - res.startLocations[0]), true);
                }
                edlibFreeAlignResult(res);
            }
        }

        #pragma omp critical
        {
            chunk_results.insert(chunk_results.end(), thread_results.begin(), thread_results.end());
        }
    }
    return chunk_results;
}

std::vector<extracted_bc> process_fastq_split_barcode(
    const std::string& input_path,
    const std::string& adapter_seq,
    int umi_length,
    const std::unordered_set<std::string>& bc1_set,
    const std::vector<int>& bc1_lengths,
    const std::unordered_set<std::string>& bc2_set,
    const std::vector<int>& bc2_lengths,
    int search_extra_min,
    int search_extra_max,
    int max_reads,
    double max_edit_distance_ratio = 0.1,
    int chunk_size = 10000,
    int num_threads = 0)
{
    std::vector<extracted_bc> results;

    if (num_threads > 0) omp_set_num_threads(num_threads);
    std::cout << "Using " << omp_get_max_threads() << " threads\n";

    std::unique_ptr<file_streaming> files_ptr;
    std::unique_ptr<read_streaming> reader_ptr;
    try {
        files_ptr  = std::make_unique<file_streaming>(input_path, std::max(1, num_threads));
        reader_ptr = std::make_unique<read_streaming>(*files_ptr);
    } catch (const std::exception& e) {
        std::cerr << "Error opening " << input_path << ": " << e.what() << "\n";
        return results;
    }

    std::string adapter_seq_rc = seq_utils::revcomp(adapter_seq);
    int reads_processed = 0, chunks_processed = 0;

    while (max_reads <= 0 || reads_processed < max_reads) {
        std::vector<read_chunk> chunk;
        chunk.reserve(chunk_size);
        for (int i = 0; i < chunk_size; ++i) {
            if (max_reads > 0 && reads_processed >= max_reads) break;
            auto rec = reader_ptr->next_sequence();
            if (!rec) break;
            reads_processed++;
            chunk.push_back({rec->id, rec->seq});
        }
        if (chunk.empty()) break;

        ++chunks_processed;
        if (chunks_processed % 100 == 0) {
            std::cout << "Processed " << reads_processed << " reads, "
                      << results.size() << " BC pairs found\n";
        }

        auto chunk_res = process_read_chunk_split_barcode(
            chunk, adapter_seq, adapter_seq_rc,
            umi_length, bc1_set, bc1_lengths, bc2_set, bc2_lengths,
            search_extra_min, search_extra_max, max_edit_distance_ratio);
        results.insert(results.end(), chunk_res.begin(), chunk_res.end());
    }

    std::cout << "Done: " << reads_processed << " reads, "
              << results.size() << " BC pairs extracted\n";
    return results;
}

// -----------------------------------------------------------------------

struct batch_entry {
    std::string input_fastq;
    std::string output_prefix;
};

static bool ensure_output_directory(const std::string& prefix) {
    std::filesystem::path p(prefix);
    std::filesystem::path dir = p.parent_path();
    if (dir.empty()) {
        return true; // current directory
    }

    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        std::cerr << "Error: Could not create output directory " << dir << ": " << ec.message() << "\n";
        return false;
    }
    return true;
}

void write_split_pair_counts(const std::vector<extracted_bc>& extracted_barcodes,
                             const std::string& output_csv) {
    std::unordered_map<std::string, uint64_t> counts;
    counts.reserve(extracted_barcodes.size());

    for (const auto& x : extracted_barcodes) {
        if (x.barcode_1.empty() || x.barcode_2.empty()) continue;
        std::string key = x.barcode_1;
        key.push_back('\t');
        key += x.barcode_2;
        ++counts[key];
    }

    std::vector<std::tuple<std::string, std::string, uint64_t>> rows;
    rows.reserve(counts.size());
    for (const auto& kv : counts) {
        auto tab = kv.first.find('\t');
        if (tab == std::string::npos) continue;
        rows.emplace_back(
            kv.first.substr(0, tab),
            kv.first.substr(tab + 1),
            kv.second
        );
    }

    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b) {
                  if (std::get<2>(a) != std::get<2>(b)) return std::get<2>(a) > std::get<2>(b);
                  if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
                  return std::get<1>(a) < std::get<1>(b);
              });

    std::ofstream out(output_csv);
    if (!out.is_open()) {
        std::cerr << "Warning: could not write split pair counts to " << output_csv << "\n";
        return;
    }
    out << "bc1,bc2,pair,count\n";
    for (const auto& row : rows) {
        const auto& bc1 = std::get<0>(row);
        const auto& bc2 = std::get<1>(row);
        const auto cnt = std::get<2>(row);
        out << bc1 << "," << bc2 << "," << bc1 << bc2 << "," << cnt << "\n";
    }
}

/// Write a dense NxM binary matrix CSV (0s and 1s) from the observed pairs.
/// Row = BC1 index, Col = BC2 index.  The file has no header.
void write_spat_mask_csv(
    const std::vector<extracted_bc>& extracted_barcodes,
    const split_whitelist& bc1_wl,
    const split_whitelist& bc2_wl,
    const std::string& output_path)
{
    const size_t nrows = bc1_wl.ordered.size();
    const size_t ncols = bc2_wl.ordered.size();
    if (nrows == 0 || ncols == 0) return;

    // Build the bit grid
    std::vector<bool> grid(nrows * ncols, false);
    size_t pairs_set = 0;
    for (const auto& x : extracted_barcodes) {
        if (x.barcode_1.empty() || x.barcode_2.empty()) continue;
        auto it1 = bc1_wl.idx.find(x.barcode_1);
        auto it2 = bc2_wl.idx.find(x.barcode_2);
        if (it1 == bc1_wl.idx.end() || it2 == bc2_wl.idx.end()) continue;
        size_t pos = static_cast<size_t>(it1->second) * ncols + it2->second;
        if (!grid[pos]) { grid[pos] = true; ++pairs_set; }
    }

    // Write row by row
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "Warning: could not write spat_mask to " << output_path << "\n";
        return;
    }
    for (size_t r = 0; r < nrows; ++r) {
        for (size_t c = 0; c < ncols; ++c) {
            if (c > 0) out << ',';
            out << (grid[r * ncols + c] ? '1' : '0');
        }
        out << '\n';
    }
    std::cout << "[spat_mask] Written " << nrows << "x" << ncols
              << " matrix (" << pairs_set << " valid pairs) to " << output_path << "\n";
}

std::vector<batch_entry> read_batch_csv(const std::string& filename) {
    std::vector<batch_entry> entries;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: Could not open batch CSV file " << filename << "\n";
        return entries;
    }

    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)) {
        line_no++;
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string fastq_path, prefix;
        std::getline(ss, fastq_path, ',');
        std::getline(ss, prefix);

        fastq_path = seq_utils::trim(fastq_path);
        prefix = seq_utils::trim(prefix);

        if (fastq_path.empty() || prefix.empty()) {
            std::cerr << "Warning: Skipping line " << line_no << " in " << filename
                      << " (need two non-empty columns)\n";
            continue;
        }

        entries.push_back({fastq_path, prefix});
    }

    if (entries.empty()) {
        std::cerr << "Warning: No valid entries found in batch CSV " << filename << "\n";
    }
    return entries;
}

// Barcode statistics structure
struct bc_wl_stats {
    std::string sequence;
    int count;
    double ncpm;
    double log1p_ncpm;
    double log1p_ncpm_ztpois;
    
    // Calculate NCPM per-barcode
    void calculate_bc_ncpm(int barcode_count, double total_reads) {
        if (total_reads > 0.0) {
            ncpm = (static_cast<double>(barcode_count) / total_reads) * 1e6;
        } else {
            ncpm = 0.0;
        }
    }
    
    // Calculate log1p_ncpm per-barcode
    void calculate_bc_log1p_ncpm() {
        log1p_ncpm = std::log1p(ncpm);
    }
    
    void calculate_bc_ztpois_pct(double k, double lambda) {
        if (k <= 0.0) {
            log1p_ncpm_ztpois = 0.0;
            return; // no zeroes allowed
        }
        if (lambda <= 0.0) {
            log1p_ncpm_ztpois = 0.0;
            return;
        }
        
        // Convert to integer for discrete distribution
        int k_int = static_cast<int>(std::floor(k));
        if (k_int <= 0) {
            log1p_ncpm_ztpois = 0.0;
            return;
        }
        
        // Pre-calculate zero probability for truncation
        double prob_zero = std::exp(-lambda);
        double truncation_denom = 1.0 - prob_zero;
        
        // Calculate regular Poisson CDF: P(X <= k)
        double regular_cdf = 0.0;
        for (int j = 0; j <= k_int; j++) {
            // Poisson PMF: P(X = j) = (lambda^j * exp(-lambda)) / j!
            double log_pmf = j * std::log(lambda) - lambda - std::lgamma(j + 1.0);
            regular_cdf += std::exp(log_pmf);
        }
        
        // Zero-truncated CDF: P(X <= k | X > 0) = [P(X <= k) - P(X = 0)] / [1 - P(X = 0)]
        double zt_cdf = (regular_cdf - prob_zero) / truncation_denom;
        
        // Clamp to [0, 1] and convert to percentage
        zt_cdf = std::max(0.0, std::min(1.0, zt_cdf));
        // Scale to percentage
        log1p_ncpm_ztpois = zt_cdf * 100.0;
    }

};

// Silverman's bandwidth calculation
double calculate_silverman_bandwidth(const std::vector<double>& data) {
    if (data.size() < 2) return 0.1;
    
    // Calculate standard deviation
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = 0.0;
    for (double val : data) {
        sq_sum += (val - mean) * (val - mean);
    }
    double std_dev = std::sqrt(sq_sum / (data.size() - 1));
    
    // Calculate IQR
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    double q1 = sorted_data[sorted_data.size() / 4];
    double q3 = sorted_data[3 * sorted_data.size() / 4];
    double iqr = q3 - q1;
    
    // Silverman's rule: 0.9 * min(sd, IQR/1.34) * n^(-1/5)
    double n = static_cast<double>(data.size());
    return 0.9 * std::min(std_dev, iqr / 1.34) * std::pow(n, -0.2);
}

// Simple kernel density estimation
std::vector<std::pair<double, double>> kde(const std::vector<double>& data, double bandwidth, int n_points = 512) {
    if (data.empty()) return {};
    
    double min_val = *std::min_element(data.begin(), data.end()) - 2 * bandwidth;
    double max_val = *std::max_element(data.begin(), data.end()) + 2 * bandwidth;
    double step = (max_val - min_val) / (n_points - 1);
    
    std::vector<std::pair<double, double>> density;
    
    for (int i = 0; i < n_points; ++i) {
        double x = min_val + i * step;
        double y = 0.0;
        
        // Gaussian kernel
        for (double data_point : data) {
            double z = (x - data_point) / bandwidth;
            y += std::exp(-0.5 * z * z);
        }
        
        y /= (data.size() * bandwidth * std::sqrt(2 * M_PI));
        density.push_back({x, y});
    }
    
    return density;
}

// Find peaks in density curve
std::vector<double> find_peaks(const std::vector<std::pair<double, double>>& density,
                               double min_height = 0.2, // default is 20%
                               bool debug = false) {
    std::vector<std::pair<double, double>> potential_peaks; // (x, height)

    // 1) Max and threshold
    double max_density = 0.0;
    for (const auto& point : density) max_density = std::max(max_density, point.second);
    double min_height_threshold = max_density * min_height; 
    if (debug) {
        std::cerr << "[find_peaks] max_density=" << std::setprecision(6) << max_density
                  << "  min_height_threshold=" << min_height_threshold << " (20%)\n";
    }

    // 2) Collect local maxima; print all local max w/ PASS/FAIL
    for (size_t i = 1; i + 1 < density.size(); ++i) {
        const double y0 = density[i-1].second;
        const double y1 = density[i].second;
        const double y2 = density[i+1].second;

        const bool is_local_max = (y1 > y0 && y1 > y2);
        if (debug && is_local_max) {
            std::cerr << "  local_max at x=" << std::fixed << std::setprecision(3) << density[i].first
                      << " y=" << std::setprecision(6) << y1
                      << " -> " << (y1 >= min_height_threshold ? "PASS" : "FAIL") << "\n";
        }
        if (is_local_max && y1 >= min_height_threshold) {
            potential_peaks.push_back({density[i].first, y1});
        }
    }

    // 3) Non-maximum suppression by separation
    std::vector<double> final_peaks;
    const double log_dist = 0.5; // doubled because floor doubled as well
    const double eps = 1e-9; // for x-comparison

    for (const auto& peak : potential_peaks) {
        bool keep_peak = true;
        for (const auto& other_peak : potential_peaks) {
            if (std::abs(peak.first - other_peak.first) < eps) continue; // same peak
            double distance = std::abs(peak.first - other_peak.first);
            if (distance <= log_dist && other_peak.second > peak.second) {
                keep_peak = false;
                break;
            }
        }
        if (keep_peak) final_peaks.push_back(peak.first);
    }

    if (debug) {
        std::cerr << "  kept peaks: ";
        for (double px : final_peaks) std::cerr << std::setprecision(3) << px << " ";
        std::cerr << "\n";
    }

    return final_peaks;
}

static inline std::vector<std::pair<double,double>>
kde_on_grid(
    const std::vector<double>& data, 
    double bandwidth,
    double min_val,
    double max_val, 
    int n_points = 512)
{
    std::vector<std::pair<double,double>> density;
    if (data.empty() || bandwidth <= 0.0 || n_points < 2) {
        return density;
    }
    density.reserve(n_points);
    const double step = (max_val - min_val) / (n_points - 1);
    const double inv_norm = 1.0 / (data.size() * bandwidth * std::sqrt(2.0 * M_PI));

    for (int i = 0; i < n_points; ++i) {
        const double x = min_val + i * step;
        double y = 0.0;
        for (double z0 : data) {
            const double z = (x - z0) / bandwidth;
            y += std::exp(-0.5 * z * z);
        }
        density.push_back({x, y * inv_norm});
    }
    return density;
}

static inline double
tau_intersection_on_grid(const std::vector<double>& xs, const std::vector<double>& LRw,
                         const std::vector<double>& BGw, double tau)
{
    const size_t n = xs.size();
    if (n == 0){
        return 0.0;
    } 
    const double c = (tau == 0.5) ? 1.0 : (1.0 - tau) / tau;

    auto diff_at = [&](size_t i){ 
        return LRw[i] - c * BGw[i]; 
    };

    // find first sign change
    for (size_t i = 0; i + 1 < n; ++i) {
        const double d0 = diff_at(i), d1 = diff_at(i+1);
        if ((d0 <= 0.0 && d1 >= 0.0) || (d0 >= 0.0 && d1 <= 0.0)) {
            // linear interp
            const double t = xs[i] + (xs[i+1] - xs[i]) * (0.0 - d0) / (d1 - d0 + 1e-300);
            return t;
        }
    }
    // no crossing: pick closest approach to zero
    size_t k = 0;
    double best = std::fabs(diff_at(0));
    for (size_t i = 1; i < n; ++i) {
        const double v = std::fabs(diff_at(i));
        if (v < best) { 
            best = v; k = i; 
        }
    }
    return xs[k];
}

static inline std::vector<double>
right_tail_area(const std::vector<double>& xs, const std::vector<double>& ys)
{
    const size_t n = xs.size();
    std::vector<double> tail(n, 0.0);
    if (n < 2) return tail;
    // segment areas
    std::vector<double> seg(n - 1, 0.0);
    for (size_t i = 0; i + 1 < n; ++i) {
        seg[i] = 0.5 * (ys[i] + ys[i+1]) * (xs[i+1] - xs[i]);
    }
    // cumulative from right
    double acc = 0.0;
    for (size_t i = n - 1; i-- > 0; ) {
        acc += seg[i];
        tail[i] = acc;
    }
    tail[n-1] = 0.0;
    return tail;
}

// ---------- AF saddle cut ----------
struct saddle_cut_result {
    bool   ok              = false;  // found >= 2 peaks and a valley
    double bw              = 0.0;
    int    n_peaks         = 0;
    double cut             = 0.0;    // raw valley x
    double final_cut       = 0.0;    // widened plateau cut (clamped by left_bound)
    bool   used_flat_widen = false;

    double valley_height   = 0.0;
    double left_peak       = 0.0;
    double right_peak      = 0.0;
    double plateau_left    = 0.0;    // widened plateau interval (for logging)
    double plateau_right   = 0.0;
};

static inline saddle_cut_result compute_saddle_cut_af(const std::vector<double>& af_x, double bw = -1.0, int n_points = 512,
                      double min_height_ratio = 0.20, bool debug = false,
                      double left_bound = -std::numeric_limits<double>::infinity()
                      ) {
    saddle_cut_result res;
    if (af_x.size() < 10) {
        return res;
    }
    // 1) Bandwidth (Silverman fallback)
    double use_bw = (bw > 0.0) ? bw : calculate_silverman_bandwidth(af_x);
    if (!(use_bw > 0.0)) use_bw = 0.1;
    res.bw = use_bw;

    // 2) KDE on AF range
    const double xmin = *std::min_element(af_x.begin(), af_x.end());
    const double xmax = *std::max_element(af_x.begin(), af_x.end());
    auto d = kde_on_grid(af_x, use_bw, xmin, xmax, n_points);
    if (d.size() < 3) return res;

    // 3) Peaks on AF KDE
    auto pkx = find_peaks(d, /*min_height=*/min_height_ratio, /*debug=*/debug);
    res.n_peaks = static_cast<int>(pkx.size());
    if (pkx.size() < 2) {
        if (debug) std::cerr << "[saddle] < 2 peaks → no saddle\n";
        return res;
    }

    auto idx_of_x = [&](double x){
        auto it = std::lower_bound(
            d.begin(), d.end(), std::make_pair(x, -INFINITY),
            [](const auto& a, const auto& b){ return a.first < b.first; });
        size_t i = (it == d.end()) ? (d.size()-1) : static_cast<size_t>(it - d.begin());
        if (i > 0 && std::fabs(d[i].first - x) > std::fabs(d[i-1].first - x)) --i;
        return i;
    };

    struct peak_h { 
        double x; 
        double h; 
        size_t i; 
    };
    std::vector<peak_h> ph;
    ph.reserve(pkx.size());
    for (double x : pkx) {
        size_t i = idx_of_x(x);
        ph.push_back({x, d[i].second, i});
    }
    std::sort(ph.begin(), ph.end(), [](const peak_h& a, const peak_h& b){ return a.h > b.h; });

    const peak_h p1 = ph[0];
    const peak_h p2 = ph[1];
    size_t i_left  = std::min(p1.i, p2.i);
    size_t i_right = std::max(p1.i, p2.i);

    // 4) Raw valley minimum between the two tallest peaks
    size_t i_valley = i_left;
    double y_min = d[i_left].second;
    for (size_t i = i_left + 1; i < i_right; ++i) {
        if (d[i].second < y_min) { 
            y_min = d[i].second; i_valley = i; 
        }
    }

    res.ok            = (i_right > i_left + 1);
    res.cut           = d[i_valley].first;
    res.valley_height = d[i_valley].second;
    res.left_peak     = d[i_left].first;
    res.right_peak    = d[i_right].first;

    if (!res.ok) return res;

    // 5) “Keep walking” across flat valley -> widen both directions,
    //    but: do not move left of left_bound; do not cross the right peak.
    const double valley_h = res.valley_height;

    // addressing peak switchoff--we expect p1 to be greater than p2 due to background
    // but in cases with spatial data, or whether a peak needs to be pruned more than it needs to be preserved,
    // the tolerance needs to be adjusted to shift left or right to whatever peak size is dominant.
    // Get heights of spatially left and right peaks
    double left_peak_height = d[i_left].second;
    double right_peak_height = d[i_right].second;

    double climb_left = left_peak_height - valley_h;
    double climb_right = right_peak_height - valley_h;

    // Stop threshold: the density level at which we stop walking
    // Whichever is lower: half peak height, or halfway up from valley
    double stop_at_left = std::min(0.5 * left_peak_height, valley_h + 0.5 * climb_left);
    double stop_at_right = std::min(0.5 * right_peak_height, valley_h + 0.5 * climb_right);

    // Left expansion - stop when density exceeds threshold
    const double left_peak_bound = std::max(left_bound, d[i_left].first);
    size_t L = i_valley;
    while (L > i_left + 1 && d[L-1].second <= stop_at_left && d[L-1].first > left_peak_bound) {
        --L;
    }

    // Right expansion - stop when density exceeds threshold
    size_t R = i_valley;
    while (R + 1 < d.size() && (R + 1) < i_right && d[R+1].second <= stop_at_right) {
        ++R;
    }


    res.plateau_left  = std::max(left_peak_bound, d[L].first);
    res.plateau_right = d[R].first;

    // Choose the leftmost point in the flat valley (“keep walking” left)
    res.final_cut = std::max(left_peak_bound, res.plateau_left);
    res.used_flat_widen = (res.final_cut != res.cut);

    if (debug) {
        std::cerr << "[saddle] peaks @ " << p1.x << ", " << p2.x
                  << "  raw valley = " << res.cut << " (h=" << valley_h << ")"
                  << "  plateau=[" << res.plateau_left << ", " << res.plateau_right << "]"
                  << "  final_cut = " << res.final_cut << "  (left_bound=" << left_peak_bound
                  << ", bw = " << res.bw << ")\n";
    }
    return res;
}

void count_perfect_matches_with_stats(const std::vector<extracted_bc>& extracted_barcodes,
                                      const std::unordered_set<int64_seq>& whitelist_set,
                                      const std::string& output_csv, const std::string text_out,
                                      bool verbose = false)
{
    std::cout << "Building hash map of extracted sequences...\n";

    // 1) Collapse read-level to barcode counts
    std::unordered_map<std::string, int> extracted_counts;
    extracted_counts.reserve(extracted_barcodes.size());

    for (const auto& x : extracted_barcodes){
        extracted_counts[x.sequence]++;
    }

    std::cout << "Found " << extracted_counts.size() << " unique sequences from "
              << extracted_barcodes.size() << " total extractions\n";

    // 2) Whitelist lookup (optional - if empty, analyze all sequences)
    bool use_whitelist = !whitelist_set.empty();
    if (use_whitelist) {
        std::cout << "Using whitelist with " << whitelist_set.size() << " barcodes\n";
    } else {
        std::cout << "No whitelist provided - analyzing all sequences\n";
    }

    // 3) Build stats (WL barcodes if provided, otherwise all)
    std::vector<bc_wl_stats> barcode_stats;
    barcode_stats.reserve(extracted_counts.size());

    double total_reads = static_cast<double>(extracted_barcodes.size());
    int total_perfect_matches = 0;
    int unique_matches = 0;

    for (const auto& kv : extracted_counts) {
        const std::string& seq = kv.first;
        int cnt = kv.second;
        
        // Skip if whitelist provided and sequence not in whitelist
        if (use_whitelist) {
            int64_seq encoded(seq);
            if (whitelist_set.find(encoded) == whitelist_set.end()) {
                continue;
            }
        }
        bc_wl_stats s;
        s.sequence = seq;
        s.count    = cnt;
        s.calculate_bc_ncpm(cnt, total_reads);
        s.calculate_bc_log1p_ncpm();
        barcode_stats.push_back(s);

        total_perfect_matches += cnt;
        unique_matches++;
    }

    // 4) zt_poisson on log1p_ncpm
    double lambda = 0.0;
    if (!barcode_stats.empty()) {
        double sum = 0.0;
        for (const auto& s : barcode_stats) sum += s.log1p_ncpm;
        lambda = sum / barcode_stats.size();
    }
    std::cout << "Calculated lambda for zero-truncated Poisson (log1p_ncpm): " << lambda << "\n";

    // ZT-Poisson percentile for each barcode
    for (auto& s : barcode_stats) {
        s.calculate_bc_ztpois_pct(s.log1p_ncpm, lambda);
    }

    // 5) Above-floor population (log1p scale)
    const double FLOOR = 2.0; // this number used to be 1,
    // and then i modified the code to handle reverse complements and now 2 works great. why? same reason as 1.
    const double ZTPOIS_RULE = 95.0; // fallback percentile if saddle/50% fail

    std::vector<double> af_x; // log1p_ncpm for AF
    af_x.reserve(barcode_stats.size());
    for (const auto& s : barcode_stats){
        if (s.log1p_ncpm >= FLOOR) {
            af_x.push_back(s.log1p_ncpm);
        }
    }

    // 6) PRIMARY THRESHOLD: AF-KDE saddle point on AF distribution
    double threshold = FLOOR;                // numeric gate (log1p_ncpm)
    std::string rule = "af_kde_saddle";
    bool have_threshold = false;

    double t_saddle = std::numeric_limits<double>::quiet_NaN();
    bool   have_saddle = false;
    saddle_cut_result sd;

    if (af_x.size() >= 10) {
        sd = compute_saddle_cut_af(af_x, /*bw=*/-1.0, /*n_points=*/512, /*min_height_ratio=*/0.20, /*debug=*/verbose);
        if (sd.ok) {
            t_saddle = sd.final_cut;
            threshold = std::max(sd.final_cut, FLOOR);
            have_threshold = true;
            have_saddle = true;
            std::cout << "\n[AF-KDE saddle] peaks @ " << sd.left_peak << " & " << sd.right_peak << "  → saddle cut=" << sd.final_cut << " (bw = " << sd.bw << ")";
        } else {
            std::cout << "\n[AF-KDE saddle] No stable valley; will try fallback.\n";
        }
    } else {
        std::cout << "\n[AF-KDE saddle] Not enough AF points (<10); will try fallback.\n";
    }

    // 7) FALLBACK THRESHOLD: ZT-Poisson 95th percentile → smallest x meeting it
    if (!have_threshold) {
        rule = "ztpois_95pct";
        double best = std::numeric_limits<double>::infinity();
        for (const auto& s : barcode_stats) {
            if (s.log1p_ncpm >= FLOOR && s.log1p_ncpm_ztpois >= ZTPOIS_RULE) {
                if (s.log1p_ncpm < best) best = s.log1p_ncpm;
            }
        }
        if (std::isfinite(best)) {
            threshold = best;
            have_threshold = true;
            std::cout << "[ZT-Poisson fallback] threshold = " << threshold << " (first x with ≥ "
                      << ZTPOIS_RULE << "%)\n";
        } else {
            threshold = FLOOR;
            std::cout << "[ZT-Poisson fallback] No items ≥ " << ZTPOIS_RULE
                      << "% — using FLOOR = " << FLOOR << "\n";
        }
    }

    // 8) DENSITY CALCS : AF mixture using a provisional LR/BG split (ZT-Poisson 95)
    //    Also compute the 50% posterior intersection (t_purity50).
    std::vector<double> lr_x, bg_x;
    lr_x.reserve(af_x.size());
    bg_x.reserve(af_x.size());
    for (const auto& s : barcode_stats) {
        if (s.log1p_ncpm <= FLOOR) {
            continue;
        }
        if (s.log1p_ncpm_ztpois >= ZTPOIS_RULE) {
            lr_x.push_back(s.log1p_ncpm);
        } else {
            bg_x.push_back(s.log1p_ncpm);
        }
    }

    double t_purity50 = std::numeric_limits<double>::quiet_NaN();
    bool   have_tpurity50 = false;

    if(af_x.size() >= 10 && lr_x.size() >= 10 && bg_x.size() >= 10) {
        double bw = calculate_silverman_bandwidth(lr_x);

        if (!(bw > 0.0)){
            bw = 0.1;
        }

        const int KDE_N = 512;
        const double xmin = *std::min_element(af_x.begin(), af_x.end());
        const double xmax = *std::max_element(af_x.begin(), af_x.end());

        auto dR = kde_on_grid(lr_x, bw, xmin, xmax, KDE_N);  // fR|AF
        auto dB = kde_on_grid(bg_x, bw, xmin, xmax, KDE_N);  // fB|AF

        std::vector<double> xs, fR, fB;
        xs.reserve(dR.size()); fR.reserve(dR.size()); fB.reserve(dB.size());
        for (size_t i = 0; i < dR.size(); ++i) {
            xs.push_back(dR[i].first);
            fR.push_back(dR[i].second);
            fB.push_back(dB[i].second);
        }

        const double piR = static_cast<double>(lr_x.size()) /
                           static_cast<double>(lr_x.size() + bg_x.size());

        // mixture-weighted components and AF
        std::vector<double> LRw(fR.size()), BGw(fB.size()), AFw(fR.size());
        for (size_t i = 0; i < fR.size(); ++i) {
            LRw[i] = piR * fR[i];
            BGw[i] = (1.0 - piR) * fB[i];
            AFw[i] = LRw[i] + BGw[i];
        }

        std::cout << "\n[AF mixture diagnostics] bw = " << bw
                  << "  |AF|=" << af_x.size()
                  << "  |LR_in_AF|=" << lr_x.size()
                  << "  |BG_in_AF|=" << bg_x.size()
                  << "  piR=" << std::setprecision(6) << piR << "\n";

        // 50% posterior intersection
        t_purity50 = tau_intersection_on_grid(xs, LRw, BGw, /*tau=*/0.5);
        have_tpurity50 = std::isfinite(t_purity50);
        if (have_tpurity50) {
            std::cout << "[t_purity50] 50% boundary at x = " << t_purity50 << "\n";
        } else {
            std::cout << "[t_purity50] could not compute.\n";
        }
        if (t_purity50 > sd.right_peak) {
            std::cout << "[t_purity50] " << t_purity50 << " exceeds right peak " 
                  << sd.right_peak << " — invalid, using saddle only\n";
            have_tpurity50 = false;
        }

        // If we have t_purity50 (and maybe a saddle), choose the greater of the two.
        if (have_tpurity50 && have_saddle) {
            double prev = threshold;
            threshold = std::max(t_saddle, t_purity50);
            rule = "max(saddle, t_purity50)";
            std::cout << "[gate] prev=" << prev << "  → chosen=max(" << t_saddle
                      << ", " << t_purity50 << ") = " << threshold << "\n";
        } else if (have_tpurity50 && !have_saddle) {
            double prev = threshold;
            threshold = std::max({t_purity50, FLOOR, prev});
            std::cout << "[gate] prev=" << prev << "  -> chosen=" << threshold << "\n";
        }
    } else {
        std::cout << "\n[AF mixture diagnostics] " << "  |AF|=" << af_x.size()
                  << "  |LR_in_AF|=" << lr_x.size() << "  |BG_in_AF|=" << bg_x.size() <<   "\n";
        std::cout << "\n[AF mixture diagnostics] Skipped (need ≥10 in AF/LR/BG).\n";
    }

    // Prepare the "between" band for annotation (only when both lines exist)
    double band_lo = std::numeric_limits<double>::infinity();
    double band_hi = -std::numeric_limits<double>::infinity();
    bool   have_band = false;
    if (have_saddle && have_tpurity50) {
        band_lo = std::min(t_saddle, t_purity50);
        band_hi = std::max(t_saddle, t_purity50);
        have_band = true;
    }

    // 9) COLLAPSE DECISION: combine left-tail fraction + spread check
    //    (a) left-tail %: if AF below threshold is small (≤10%), return ALL AF
    //    (b) spread: if both sides are narrow (≤1.0) OR full range ≤ 2.0, return ALL AF
    // 9b) LEFT_BUDGET fallback: if collapse did not fire but left_frac is only
    //     marginally over the collapse threshold (10% < left_frac ≤ 20%),
    //     rescue the top LEFT_TAIL_MAX_FRAC * af_total barcodes from AF_left
    //     (sorted by count desc). This recovers real cells that fall just
    //     under the ZT-Poisson threshold without admitting the deep noise tail.
    //     Above LEFT_TAIL_BUDGET_MAX_FRAC (20%) the threshold is treated as
    //     unreliable / the data is too noisy and we fall back to the strict
    //     "drop AF_left" behavior.
    const double LEFT_TAIL_MAX_FRAC        = 0.10;  // 10%
    const double LEFT_TAIL_BUDGET_MAX_FRAC = 0.20;  // 20% — upper bound for budget fallback
    const double SPREAD_EPS = 1.0;           // in log1p_ncpm units
    const int    MIN_SIDE_N = 5;             // avoid tiny/empty groups

    size_t af_total = 0, af_left = 0, af_right = 0;
    double left_min  =  std::numeric_limits<double>::infinity();
    double right_min =  std::numeric_limits<double>::infinity();
    double left_max  = -std::numeric_limits<double>::infinity();
    double right_max = -std::numeric_limits<double>::infinity();

    for (const auto& s : barcode_stats) {
        if (s.log1p_ncpm < FLOOR) continue;
        af_total++;
        if (s.log1p_ncpm < threshold) {
            af_left++;
            left_min  = std::min(left_min,  s.log1p_ncpm);
            left_max  = std::max(left_max,  s.log1p_ncpm);
        } else {
            af_right++;
            right_min = std::min(right_min, s.log1p_ncpm);
            right_max = std::max(right_max, s.log1p_ncpm);
        }
    }

    const double left_frac   = (af_total > 0) ? (static_cast<double>(af_left) / af_total) : 1.0;
    const bool   have_left   = (af_left  >= MIN_SIDE_N);
    const bool   have_right  = (af_right >= MIN_SIDE_N);
    const double left_range  = (have_left  ? (left_max  - left_min)  : 0.0);
    const double right_range = (have_right ? (right_max - right_min) : 0.0);
    const double full_range  = (std::isfinite(left_min) && std::isfinite(right_max))
                                 ? (right_max - left_min)
                                 : 0.0;

    const bool narrow_both   = have_left && have_right &&
                               (left_range <= SPREAD_EPS) && (right_range <= SPREAD_EPS);
    const bool narrow_full   = (full_range > 0.0) && (full_range <= 2.0 * SPREAD_EPS);
    const bool collapse_all_af = (af_total >= 10) &&
                                 ((left_frac <= LEFT_TAIL_MAX_FRAC) || narrow_both || narrow_full);

    // Budget rescue: only fires when collapse didn't, and left_frac is in the
    // borderline (LEFT_TAIL_MAX_FRAC, LEFT_TAIL_BUDGET_MAX_FRAC] window.
    const bool use_left_budget = (af_total >= 10) && !collapse_all_af
                              && (left_frac >  LEFT_TAIL_MAX_FRAC)
                              && (left_frac <= LEFT_TAIL_BUDGET_MAX_FRAC);
    std::unordered_set<std::string> budget_rescue;
    size_t budget_n = 0;
    if (use_left_budget) {
        std::vector<size_t> af_left_idx;
        af_left_idx.reserve(af_left);
        for (size_t i = 0; i < barcode_stats.size(); ++i) {
            const auto& s = barcode_stats[i];
            if (s.log1p_ncpm >= FLOOR && s.log1p_ncpm < threshold) {
                af_left_idx.push_back(i);
            }
        }
        budget_n = static_cast<size_t>(LEFT_TAIL_MAX_FRAC * static_cast<double>(af_total));
        if (budget_n > af_left_idx.size()) budget_n = af_left_idx.size();
        if (budget_n > 0) {
            std::partial_sort(
                af_left_idx.begin(),
                af_left_idx.begin() + budget_n,
                af_left_idx.end(),
                [&](size_t a, size_t b) {
                    return barcode_stats[a].count > barcode_stats[b].count;
                });
            for (size_t i = 0; i < budget_n; ++i) {
                budget_rescue.insert(barcode_stats[af_left_idx[i]].sequence);
            }
        }
    }

    std::cout << "\n[AF split] floor=" << FLOOR
              << "  threshold=" << threshold
              << "  rule=" << rule
              << "\n         AF_total=" << af_total
              << "  AF_left=" << af_left << " (" << left_frac * 100.0 << "%)"
              << "  AF_right=" << af_right
              << "\n[Spread]  left_range="  << left_range
              << "  right_range=" << right_range
              << "  full_range="  << full_range
              << "\n[Collapse] "
              << (collapse_all_af
                    ? "YES → return ALL AF"
                    : (use_left_budget
                          ? ("NO; LEFT_BUDGET fallback → +"
                             + std::to_string(budget_n)
                             + " AF_left rescues (top "
                             + std::to_string(static_cast<int>(LEFT_TAIL_MAX_FRAC * 100.0))
                             + "% of AF_total by count)")
                          : "NO"))
              << "\n";

    // 10) Write CSV — explicit outputs + annotation band
    //     final_barcode = TRUE iff:
    //       - collapse_all_af: (log1p_ncpm >= FLOOR)
    //       - else:            (log1p_ncpm >= threshold)
    std::ofstream csv_out(output_csv);
    std::vector<std::string> high_conf_barcodes;

    csv_out << "barcode,count,ncpm,log1p_ncpm,ztpois_percentile,above_floor,over_threshold,"
               "final_barcode,final_bc_annotation\n";

    size_t n_final = 0, n_above_floor = 0;
    for (const auto& s : barcode_stats) {
        const bool above_floor   = (s.log1p_ncpm >= FLOOR);
        if (above_floor) n_above_floor++;

        const bool over_threshold = (s.log1p_ncpm >= threshold);
        const bool budget_pass    = use_left_budget &&
                                    above_floor && !over_threshold &&
                                    (budget_rescue.find(s.sequence) != budget_rescue.end());
        const bool final_barcode  = collapse_all_af
                                      ? above_floor
                                      : (over_threshold || budget_pass);
        if (final_barcode) n_final++;

        // Annotation:
        // - collapse mode: AF that pass are "high_confidence"
        // - otherwise:
        //     * ≥ chosen threshold: "high_confidence"
        //     * budget rescue (10-20% borderline window): "high_sensitivity"
        //     * between t_saddle and t_purity50: "high_sensitivity"
        std::string ann;
        if (collapse_all_af) {
            if (final_barcode){
                ann = "high_confidence";
                high_conf_barcodes.push_back(s.sequence);
            }
        } else {
            if (over_threshold && above_floor) {
                ann = "high_confidence";
                high_conf_barcodes.push_back(s.sequence);
            } else if (budget_pass) {
                ann = "high_sensitivity";
                high_conf_barcodes.push_back(s.sequence);
            } else if (have_band && above_floor && s.log1p_ncpm >= band_lo && s.log1p_ncpm < band_hi) {
                ann = "high_sensitivity";
            } else {
                ann = "low_confidence";
            }
        }

        csv_out << s.sequence << ","
                << s.count    << ","
                << std::fixed << std::setprecision(6) << s.ncpm << ","
                << s.log1p_ncpm << ","
                << s.log1p_ncpm_ztpois << ","
                << (above_floor   ? "TRUE" : "FALSE") << ","
                << (over_threshold? "TRUE" : "FALSE") << ","
                << (final_barcode ? "TRUE" : "FALSE") << ","
                << ann
                << "\n";
    }

    // 11) Console summary
    std::cout << "\nStatistical analysis complete!\n"
              << "Total extracted sequences: " << extracted_barcodes.size() << "\n"
              << "Unique extracted sequences: " << extracted_counts.size() << "\n"
              << "Total perfect matches: " << total_perfect_matches << "\n"
              << "Unique barcodes with perfect matches: " << unique_matches << "\n"
              << "Match rate: " << std::fixed << std::setprecision(2)
              << (100.0 * total_perfect_matches / std::max<size_t>(1, extracted_barcodes.size())) << "%\n"
              << "Above-floor barcodes: " << n_above_floor << "\n"
              << "Final barcodes (TRUE): " << n_final << "\n"
              << "Results written to: " << output_csv << "\n";

    if (!text_out.empty()) {
        std::ofstream hc_out(text_out);
        for (const auto& bc : high_conf_barcodes) {
            hc_out << bc << "\n";
        }
        std::cout << "High-confidence barcodes written to: " << text_out << "\n";
    }
}


void usage_scan_wl(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  -i, --input FILE            Input FASTQ file (use this for single runs)\n"
              << "  -b, --batch-csv FILE        CSV with FASTQ path and output prefix (two columns, optional)\n"
              << "  -o, --output-prefix PREFIX  Output prefix for generated files (.txt and .csv)\n"
              << "  -p, --adapter_seq SEQ       Adapter/primer sequence to search for\n"
              << "  -n, --barcode-length INT    Number of bases to extract (barcode length)\n"
              << "  -l, --left-margin INT       Bases to include on left side [default: 0]\n"
              << "  -r, --right-margin INT      Bases to include on right side [default: 0]\n"
              << "  -m, --max-reads INT         Maximum number of reads to process [default: all]\n"
              << "  -e, --max-error FLOAT       Maximum edit distance ratio [default: 0.1]\n"
              << "  -w, --whitelist FILE        Barcode whitelist kit key or path (optional)\n"
              << "  -t, --threads INT           Number of threads for parallel processing [default: auto]\n"
              << "  -k, --chunk-size INT        Chunk size for parallel processing [default: 10000]\n"
              << "  -v, --verbose               Enable verbose/debug output\n"
              << "  -h, --help                  Show this help message\n"
              << "\nTwo-part barcode options (BC1+BC2 mode):\n"
              << "  -1, --bc1-whitelist FILE    BC1 whitelist kit key or path (activates two-part mode)\n"
              << "  -2, --bc2-whitelist FILE    BC2 whitelist kit key or path\n"
              << "  -u, --umi-length INT        UMI length between primer and BC1 [default: 9]\n"
              << "      --offset-min INT        Min extra bases after UMI to start BC1 search [default: 0]\n"
              << "      --offset-max INT        Max extra bases after UMI to start BC1 search [default: 3]\n";
}

inline int cmd_scan_wl(int argc, char* argv[]) {
    optind = 1; // reset getopt state for subcommand dispatch
    std::string input_file, output_prefix, adapter_seq, whitelist_file, batch_csv_file;
    std::string bc1_whitelist_file, bc2_whitelist_file;
    int bc_length = 16, m_left = 0, m_right = 0, max_reads = 0, num_threads = 0, chunk_size = 10000;
    int umi_length = 9, offset_min = 0, offset_max = 3;
    double max_error = 0.1;
    bool verbose = false;

    static struct option long_options[] = {
        {"input",          required_argument, 0, 'i'},
        {"batch-csv",      required_argument, 0, 'b'},
        {"output-prefix",  required_argument, 0, 'o'},
        {"adapter_seq",    required_argument, 0, 'p'},
        {"barcode-length", required_argument, 0, 'n'},
        {"left-margin",    required_argument, 0, 'l'},
        {"right-margin",   required_argument, 0, 'r'},
        {"max-reads",      required_argument, 0, 'm'},
        {"max-error",      required_argument, 0, 'e'},
        {"whitelist",      optional_argument, 0, 'w'},
        {"threads",        required_argument, 0, 't'},
        {"chunk-size",     required_argument, 0, 'k'},
        {"bc1-whitelist",  required_argument, 0, '1'},
        {"bc2-whitelist",  required_argument, 0, '2'},
        {"umi-length",     required_argument, 0, 'u'},
        {"offset-min",     required_argument, 0,  3 },
        {"offset-max",     required_argument, 0,  4 },
        // Backward-compatible aliases
        {"hd-offset-min",  required_argument, 0,  3 },
        {"hd-offset-max",  required_argument, 0,  4 },
        {"verbose",        no_argument,       0, 'v'},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:b:o:p:n:l:r:m:e:w:t:k:1:2:u:hv", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'b': batch_csv_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 'p': adapter_seq = optarg; break;
            case 'n': bc_length = std::atoi(optarg); break;
            case 'l': m_left = std::atoi(optarg); break;
            case 'r': m_right = std::atoi(optarg); break;
            case 'm': max_reads = std::atoi(optarg); break;
            case 'e': max_error = std::atof(optarg); break;
            case 'w': whitelist_file = optarg; break;
            case 't': num_threads = std::atoi(optarg); break;
            case 'k': chunk_size = std::atoi(optarg); break;
            case '1': bc1_whitelist_file = optarg; break;
            case '2': bc2_whitelist_file = optarg; break;
            case 'u': umi_length = std::atoi(optarg); break;
            case  3 : offset_min = std::atoi(optarg); break;
            case  4 : offset_max = std::atoi(optarg); break;
            case 'v': verbose = true; break;
            case 'h': usage_scan_wl(argv[0]); return 0;
            default: usage_scan_wl(argv[0]); return 1;
        }
    }
    
    if (adapter_seq.empty()) {
        std::cerr << "Error: Missing required adapter sequence\n";
        usage_scan_wl(argv[0]);
        return 1;
    }

    const bool using_single = !input_file.empty();
    const bool using_batch  = !batch_csv_file.empty();

    if (!using_single && !using_batch) {
        std::cerr << "Error: Provide either --input (single run) or --batch-csv (batch mode)\n";
        usage_scan_wl(argv[0]);
        return 1;
    }

    if (using_batch && using_single) {
        std::cout << "Batch CSV provided - ignoring single --input value\n";
    }

    const bool split_barcode_mode = !bc1_whitelist_file.empty() && !bc2_whitelist_file.empty();

    auto run_single = [&](const std::string& fastq_path, const std::string& output_prefix) {
        std::string resolved_prefix = output_prefix.empty() ? "barcodes" : output_prefix;
        if (!ensure_output_directory(resolved_prefix)) {
            return;
        }

        std::string csv_output  = resolved_prefix + ".csv";
        std::string text_output = resolved_prefix + ".txt";

        if (split_barcode_mode) {
            std::string bc1_path, bc2_path;
            try {
                bc1_path = whitelist_utils::kit_to_path(bc1_whitelist_file);
                bc2_path = whitelist_utils::kit_to_path(bc2_whitelist_file);
            } catch (const std::exception& e) {
                std::cerr << "Error: could not resolve split whitelist path(s): " << e.what() << "\n";
                return;
            }

            std::cout << "Two-part barcode mode: " << fastq_path << "\n"
                      << "  Adapter: " << adapter_seq << "  UMI length: " << umi_length
                      << "  BC1 search offset: [" << offset_min << ", " << offset_max << "]\n";
            if (verbose) {
                std::cout << "  BC1 whitelist: " << bc1_path << "\n"
                          << "  BC2 whitelist: " << bc2_path << "\n";
            }

            auto bc1_wl = load_split_whitelist_csv(bc1_path);
            auto bc2_wl = load_split_whitelist_csv(bc2_path);
            if (bc1_wl.seqs.empty() || bc2_wl.seqs.empty()) {
                std::cerr << "Error: failed to load two-part barcode whitelists\n";
                return;
            }
            std::cout << "Loaded " << bc1_wl.seqs.size() << " BC1 / " << bc2_wl.seqs.size() << " BC2 sequences\n";

            auto barcodes = process_fastq_split_barcode(
                fastq_path, adapter_seq,
                umi_length, bc1_wl.seqs, bc1_wl.lengths, bc2_wl.seqs, bc2_wl.lengths,
                offset_min, offset_max,
                max_reads, max_error, chunk_size, num_threads);

            std::string pair_counts_output = resolved_prefix + "_valid_pairs.csv";
            write_split_pair_counts(barcodes, pair_counts_output);
            std::cout << "Raw valid BC1+BC2 pair counts written to: " << pair_counts_output << "\n";

            std::string spat_mask_output = resolved_prefix + "_spat_mask.csv";
            write_spat_mask_csv(barcodes, bc1_wl, bc2_wl, spat_mask_output);

            std::unordered_set<int64_seq> empty_wl;
            count_perfect_matches_with_stats(barcodes, empty_wl, csv_output, text_output, verbose);
        } else {
            std::cout << "Processing " << fastq_path << "...\n";
            std::cout << "Adapter: " << adapter_seq << "\n";
            std::cout << "Extracting " << bc_length << " bases with margins: left=" << m_left << ", right=" << m_right << "\n";
            std::cout << "Chunk size: " << chunk_size << " reads per chunk\n";

            auto barcodes = process_fastq(
                fastq_path, adapter_seq, bc_length, m_left, m_right,
                max_reads, max_error, chunk_size, num_threads);

            std::cout << "Found " << barcodes.size() << " barcodes\n";

            if (!whitelist_file.empty()) {
                std::string resolved_wl;
                try {
                    resolved_wl = whitelist_utils::kit_to_path(whitelist_file);
                } catch (const std::exception& e) {
                    std::cerr << "Error: could not resolve whitelist path: " << e.what() << "\n";
                    return;
                }
                std::cout << "\nLoading whitelist and performing perfect match analysis...\n";
                if (verbose) {
                    std::cout << "Whitelist source: " << whitelist_file << "\n"
                              << "Resolved path: " << resolved_wl << "\n";
                }
                auto wl_set = ::whitelist::load_barcodes(resolved_wl, bc_length, verbose);
                std::cout << "Loaded " << wl_set.size() << " whitelist barcodes\n";
                count_perfect_matches_with_stats(barcodes, wl_set, csv_output, text_output, verbose);
            } else {
                std::unordered_set<int64_seq> empty_whitelist;
                std::cout << " No whitelist provided; generating whitelist de novo.\n";
                count_perfect_matches_with_stats(barcodes, empty_whitelist, csv_output, text_output, verbose);
            }
        }
    };

    if (using_batch) {
        auto batch_entries = read_batch_csv(batch_csv_file);
        if (batch_entries.empty()) {
            return 1;
        }
        std::cout << "Running batch of " << batch_entries.size() << " FASTQ files from " << batch_csv_file << "\n";
        for (size_t idx = 0; idx < batch_entries.size(); ++idx) {
            const auto& entry = batch_entries[idx];
            std::cout << "\n[Batch " << (idx + 1) << "/" << batch_entries.size() << "] Prefix: "
                      << entry.output_prefix << "\n";
            run_single(entry.input_fastq, entry.output_prefix);
        }
    } else {
        run_single(input_file, output_prefix);
    }
    return 0;
}
