#include "include/rad/rad_headers.h"


class BarcodeErrorModel {
private:
    int k_size;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> kmer_alternatives;
    std::mutex kmer_mutex;  // For thread-safe k-mer updates
    
    // Helper function to check if sequence contains only ACGT
    bool is_valid_sequence(const std::string& seq) {
        for (char c : seq) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                return false;
            }
        }
        return true;
    }

private:
    // Random subsampling helper
    template<typename T>
    std::vector<T> random_subsample(const std::vector<T>& input, int max_size, unsigned int seed = 42) {
        if (max_size <= 0 || input.size() <= max_size) {
            return input;  // Return all if max_size is larger than input
        }
        
        std::vector<T> result = input;  // Copy the input
        
        // Use random_device for better randomness, but fall back to seed if needed
        std::random_device rd;
        std::mt19937 gen(rd.entropy() ? rd() : seed);
        
        // Shuffle and resize
        std::shuffle(result.begin(), result.end(), gen);
        result.resize(max_size);
        
        return result;
    }
    // Parse CIGAR string and create aligned sequences
    void parse_cigar_to_alignment(const std::string& cigar, const std::string& ref, 
                                 const std::string& read, std::string& aligned_ref, 
                                 std::string& aligned_read) {
        int ref_pos = 0, read_pos = 0;
        aligned_ref.clear();
        aligned_read.clear();
        
        for (size_t i = 0; i < cigar.length(); ) {
            int num = 0;
            while (i < cigar.length() && std::isdigit(cigar[i])) {
                num = num * 10 + (cigar[i] - '0');
                i++;
            }
            
            if (i >= cigar.length()) break;
            char op = cigar[i++];
            
            switch (op) {
                case '=':  // Match
                case 'X':  // Mismatch
                case 'M':  // Match/Mismatch (ambiguous)
                    if (ref_pos + num <= ref.length() && read_pos + num <= read.length()) {
                        aligned_ref += ref.substr(ref_pos, num);
                        aligned_read += read.substr(read_pos, num);
                        ref_pos += num;
                        read_pos += num;
                    }
                    break;
                    
                case 'I':  // Insertion in read
                    if (read_pos + num <= read.length()) {
                        aligned_ref += std::string(num, '-');
                        aligned_read += read.substr(read_pos, num);
                        read_pos += num;
                    }
                    break;
                    
                case 'D':  // Deletion in read
                    if (ref_pos + num <= ref.length()) {
                        aligned_ref += ref.substr(ref_pos, num);
                        aligned_read += std::string(num, '-');
                        ref_pos += num;
                    }
                    break;
            }
        }
    }
    
    // Align two sequences using edlib
    std::pair<std::string, std::string> align_with_edlib(const std::string& ref, const std::string& read) {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);
        EdlibAlignResult result = edlibAlign(read.c_str(), read.length(), 
                                           ref.c_str(), ref.length(), config);
        
        if (result.status != EDLIB_STATUS_OK || result.alignmentLength == 0) {
            edlibFreeAlignResult(result);
            return {"", ""};
        }
        
        // Convert alignment to CIGAR
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        
        if (!cigar) {
            edlibFreeAlignResult(result);
            return {"", ""};
        }
        
        std::string cigar_str(cigar);
        std::string aligned_ref, aligned_read;
        parse_cigar_to_alignment(cigar_str, ref, read, aligned_ref, aligned_read);
        
        free(cigar);
        edlibFreeAlignResult(result);
        
        return {aligned_ref, aligned_read};
    }
    
    // Thread-safe k-mer error extraction
    void extract_kmer_errors_threadsafe(const std::string& aligned_ref, const std::string& aligned_read) {
        // Build local k-mer map first (no locks needed)
        std::unordered_map<std::string, std::unordered_map<std::string, int>> local_kmer_alternatives;
        
        int start = 0;
        while (start < static_cast<int>(aligned_ref.length())) {
            // Find k_size non-gap characters in reference
            std::string ref_kmer = "";
            std::string read_kmer = "";
            int end = start;
            int ref_bases = 0;
            
            // Expand window until we have k_size reference bases
            while (end < static_cast<int>(aligned_ref.length()) && ref_bases < k_size) {
                if (aligned_ref[end] != '-') {
                    ref_kmer += aligned_ref[end];
                    ref_bases++;
                }
                if (aligned_read[end] != '-') {
                    read_kmer += aligned_read[end];
                }
                end++;
            }
            
            // Only count if we got a complete k-mer
            if (ref_bases == k_size && !ref_kmer.empty() && !read_kmer.empty() && read_kmer.length() > 1) {
                // Quality checks: first/last bases match, only ACGT
                if (ref_kmer.front() == read_kmer.front() && 
                    ref_kmer.back() == read_kmer.back() &&
                    is_valid_sequence(ref_kmer) && is_valid_sequence(read_kmer)) {
                    
                    local_kmer_alternatives[ref_kmer][read_kmer]++;
                }
            }
            
            // Move to next position (skip gaps)
            start++;
            while (start < static_cast<int>(aligned_ref.length()) && aligned_ref[start] == '-') {
                start++;
            }
        }
        
        // Now update global map with lock
        if (!local_kmer_alternatives.empty()) {
            std::lock_guard<std::mutex> lock(kmer_mutex);
            for (const auto& [ref_kmer, alternatives] : local_kmer_alternatives) {
                for (const auto& [read_kmer, count] : alternatives) {
                    kmer_alternatives[ref_kmer][read_kmer] += count;
                }
            }
        }
    }
    
    // Load FASTA sequences
    std::vector<std::pair<std::string, std::string>> load_fasta(const std::string& filename) {
        std::vector<std::pair<std::string, std::string>> sequences;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open FASTA file: " << filename << std::endl;
            return sequences;
        }
        
        std::string line, current_id, current_seq;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous sequence if exists
                if (!current_id.empty() && !current_seq.empty()) {
                    sequences.push_back({current_id, current_seq});
                }
                
                // Start new sequence
                current_id = line.substr(1);
                // Remove everything after first space
                size_t space_pos = current_id.find(' ');
                if (space_pos != std::string::npos) {
                    current_id = current_id.substr(0, space_pos);
                }
                current_seq.clear();
            } else {
                // Convert to uppercase and add to sequence
                for (char c : line) {
                    if (!std::isspace(c)) {
                        current_seq += std::toupper(c);
                    }
                }
            }
        }
        
        // Don't forget the last sequence
        if (!current_id.empty() && !current_seq.empty()) {
            sequences.push_back({current_id, current_seq});
        }
        
        file.close();
        return sequences;
    }
    
    // Load FASTQ sequences
    std::vector<std::pair<std::string, std::string>> load_fastq(const std::string& filename) {
        std::vector<std::pair<std::string, std::string>> sequences;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open FASTQ file: " << filename << std::endl;
            return sequences;
        }
        
        std::string line, read_id, read_seq, plus_line, quality;
        int line_count = 0;
        
        while (std::getline(file, line)) {
            line_count++;
            
            switch (line_count % 4) {
                case 1:  // Header line (@)
                    if (!line.empty() && line[0] == '@') {
                        read_id = line.substr(1);
                        // Remove everything after first space
                        size_t space_pos = read_id.find(' ');
                        if (space_pos != std::string::npos) {
                            read_id = read_id.substr(0, space_pos);
                        }
                    }
                    break;
                case 2:  // Sequence line
                    read_seq = line;
                    // Convert to uppercase
                    std::transform(read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper);
                    break;
                case 3:  // Plus line (+)
                    plus_line = line;
                    break;
                case 0:  // Quality line
                    quality = line;
                    // Store the sequence (we ignore quality for now)
                    if (!read_id.empty() && !read_seq.empty()) {
                        sequences.push_back({read_id, read_seq});
                    }
                    break;
            }
        }
        
        file.close();
        return sequences;
    }
    
    // Calculate alignment score (simple: count matches)
    int calculate_alignment_score(const std::string& aligned_ref, const std::string& aligned_read) {
        int score = 0;
        int min_len = std::min(aligned_ref.length(), aligned_read.length());
        
        for (int i = 0; i < min_len; i++) {
            if (aligned_ref[i] != '-' && aligned_read[i] != '-' && aligned_ref[i] == aligned_read[i]) {
                score++;
            }
        }
        return score;
    }
    
    // Progress bar helper
    void print_progress_bar(int current, int total, const std::string& current_ref = "", int bar_width = 50) {
        double percentage = static_cast<double>(current) / total;
        int filled = static_cast<int>(bar_width * percentage);
        
        std::cout << "\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < filled) std::cout << "█";
            else std::cout << " ";
        }
        std::cout << "] " << std::fixed << std::setprecision(1) 
                  << percentage * 100.0 << "% (" << current << "/" << total << ")";
        
        if (!current_ref.empty()) {
            std::cout << " | Current: " << current_ref;
        }
        
        std::cout.flush();
    }
    
public:
    BarcodeErrorModel(int k) : k_size(k) {}
    
    // Process dataset by iterating through reference sequences
void process_dataset(const std::string& read_file, const std::string& ref_fasta, int num_threads,
                    int max_reads = -1, bool is_fastq = true,
                    int min_alignment_score = 5) {
    
    // Set number of threads
    omp_set_num_threads(num_threads);
    
    // Load reference sequences
    auto ref_sequences = load_fasta(ref_fasta);
    
    if (ref_sequences.empty()) {
        std::cerr << "Error: No reference sequences found" << std::endl;
        return;
    }
    
    // Load ALL read sequences (we'll subsample per reference)
    std::vector<std::pair<std::string, std::string>> all_read_sequences;
    if (is_fastq) {
        all_read_sequences = load_fastq(read_file);
    } else {
        all_read_sequences = load_fasta(read_file);
    }
    
    if (all_read_sequences.empty()) {
        std::cerr << "Error: No read sequences found" << std::endl;
        return;
    }
    
    std::cout << "Loaded " << ref_sequences.size() << " reference sequences" << std::endl;
    std::cout << "Loaded " << all_read_sequences.size() << " total read sequences" << std::endl;
    std::cout << "Using " << num_threads << " threads" << std::endl;
    std::cout << "Minimum alignment score: " << min_alignment_score << std::endl;
    
    if (max_reads > 0) {
        std::cout << "Will subsample " << max_reads << " reads per reference sequence" << std::endl;
    } else {
        std::cout << "Using all reads for each reference sequence" << std::endl;
    }
    
    int total_refs = ref_sequences.size();
    int total_successful_alignments = 0;
    int total_processed_alignments = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "\nProcessing reference sequences..." << std::endl;
    print_progress_bar(0, total_refs);
    
    // Iterate through each reference sequence
    for (int ref_idx = 0; ref_idx < total_refs; ++ref_idx) {
        // Extract values BEFORE parallel region to avoid structured binding capture issue
        std::string current_ref_id = ref_sequences[ref_idx].first;
        std::string current_ref_seq = ref_sequences[ref_idx].second;
        
        // Subsample reads for THIS reference
        std::vector<std::pair<std::string, std::string>> reads_for_this_ref;
        if (max_reads > 0 && all_read_sequences.size() > max_reads) {
            reads_for_this_ref = random_subsample(all_read_sequences, max_reads);
        } else {
            reads_for_this_ref = all_read_sequences;  // Use all reads
        }
        
        int successful_alignments_for_ref = 0;
        int processed_alignments_for_ref = 0;
        
        // Process the subsampled reads against this reference in parallel
        #pragma omp parallel
        {
            int local_successful = 0;
            int local_processed = 0;
            
            #pragma omp for schedule(dynamic)
            for (int read_idx = 0; read_idx < static_cast<int>(reads_for_this_ref.size()); ++read_idx) {
                // Extract read values (avoid structured binding in parallel region)
                std::string current_read_id = reads_for_this_ref[read_idx].first;
                std::string current_read_seq = reads_for_this_ref[read_idx].second;
                
                local_processed++;
                
                // Align this read to current reference
                auto alignment = align_with_edlib(current_ref_seq, current_read_seq);
                
                if (!alignment.first.empty()) {
                    int score = calculate_alignment_score(alignment.first, alignment.second);
                    
                    if (score >= min_alignment_score) {
                        extract_kmer_errors_threadsafe(alignment.first, alignment.second);
                        local_successful++;
                    }
                }
            }
            
            // Accumulate thread-local counts
            #pragma omp atomic
            successful_alignments_for_ref += local_successful;
            
            #pragma omp atomic
            processed_alignments_for_ref += local_processed;
        }
        
        // Update global counters
        total_successful_alignments += successful_alignments_for_ref;
        total_processed_alignments += processed_alignments_for_ref;
        
        // Update progress bar continuously
        print_progress_bar(ref_idx + 1, total_refs, current_ref_id);
    }
    
    std::cout << std::endl;  // New line after progress bar
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\nProcessing complete!" << std::endl;
    std::cout << "Total references processed: " << total_refs << std::endl;
    std::cout << "Total alignments attempted: " << total_processed_alignments << std::endl;
    std::cout << "Successful alignments: " << total_successful_alignments << std::endl;
    std::cout << "Success rate: " << std::fixed << std::setprecision(1) 
              << (100.0 * total_successful_alignments / total_processed_alignments) << "%" << std::endl;
    std::cout << "Total time: " << total_elapsed.count() << " seconds" << std::endl;
    std::cout << "Average rate: " << std::fixed << std::setprecision(2) 
              << static_cast<double>(total_refs) / total_elapsed.count() << " refs/sec" << std::endl;
    std::cout << "K-mers found: " << kmer_alternatives.size() << std::endl;
    
    // Print some statistics
    int total_observations = 0;
    for (const auto& [kmer, alternatives] : kmer_alternatives) {
        for (const auto& [alt, count] : alternatives) {
            total_observations += count;
        }
    }
    std::cout << "Total k-mer observations: " << total_observations << std::endl;
    
    // Calculate expected total alignments for verification
    long long expected_total = static_cast<long long>(total_refs) * 
                              (max_reads > 0 ? std::min(max_reads, static_cast<int>(all_read_sequences.size())) 
                                            : static_cast<int>(all_read_sequences.size()));
    std::cout << "Expected total alignments: " << expected_total << std::endl;
}

    // Output error model in BadRead format
    void output_error_model(const std::string& output_file) {
        std::ofstream out(output_file);
        
        if (!out.is_open()) {
            std::cerr << "Error: Cannot open output file: " << output_file << std::endl;
            return;
        }
        
        int kmers_written = 0;
        
        for (const auto& [kmer, alternatives] : kmer_alternatives) {
            if (alternatives.empty()) continue;
            
            int total = 0;
            for (const auto& [alt, count] : alternatives) {
                total += count;
            }
            
            if (total == 0) continue;
            
            // Output format: KMER,correct_prob;ALT1,prob1;ALT2,prob2;
            out << kmer << ",";
            
            // Correct probability (if kmer maps to itself)
            auto it = alternatives.find(kmer);
            if (it != alternatives.end()) {
                out << std::fixed << std::setprecision(6) << (double)it->second / total;
            } else {
                out << "0.000000";
            }
            
            // Sort alternatives by frequency (excluding self-mapping)
            std::vector<std::pair<std::string, double>> alt_probs;
            for (const auto& [alt, count] : alternatives) {
                if (alt != kmer) {
                    alt_probs.push_back({alt, (double)count / total});
                }
            }
            
            std::sort(alt_probs.begin(), alt_probs.end(), 
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            for (const auto& [alt, prob] : alt_probs) {
                out << ";" << alt << "," << std::fixed << std::setprecision(6) << prob;
            }
            out << std::endl;
            kmers_written++;
        }
        
        out.close();
        std::cout << "Error model written to: " << output_file << std::endl;
        std::cout << "K-mers written: " << kmers_written << std::endl;
    }
};

// Main function for dataset processing
void process_barcode_dataset(const std::string& read_file, const std::string& ref_fasta, 
                           int k_size, const std::string& output_file, 
                           int max_reads = -1, bool is_fastq = true, int num_threads = 4,
                           int min_alignment_score = 5) {
    
    std::cout << "=== Barcode Error Model - Reference-Based Processing ===" << std::endl;
    std::cout << "Read file: " << read_file << std::endl;
    std::cout << "Reference FASTA: " << ref_fasta << std::endl;
    std::cout << "K-mer size: " << k_size << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "Max reads: " << (max_reads > 0 ? std::to_string(max_reads) : "unlimited") << std::endl;
    std::cout << "File type: " << (is_fastq ? "FASTQ" : "FASTA") << std::endl;
    std::cout << "Threads: " << (num_threads > 0 ? std::to_string(num_threads) : "auto") << std::endl;
    std::cout << "Min alignment score: " << min_alignment_score << std::endl;
    std::cout << "=========================================================" << std::endl;
    
    BarcodeErrorModel model(k_size);
    model.process_dataset(read_file, ref_fasta, num_threads, max_reads, is_fastq, min_alignment_score);
    model.output_error_model(output_file);
}

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 9) {
        std::cout << "Usage: " << argv[0] << " <read_file> <ref_fasta> <k_size> <output_file> [max_reads] [file_type] [threads] [min_score]" << std::endl;
        std::cout << "  read_file: FASTQ/FASTA file with reads to analyze" << std::endl;
        std::cout << "  ref_fasta: FASTA file with reference sequences" << std::endl;
        std::cout << "  k_size: K-mer size (e.g., 5)" << std::endl;
        std::cout << "  output_file: Output error model file" << std::endl;
        std::cout << "  max_reads: Maximum reads to process (optional, -1 for all)" << std::endl;
        std::cout << "  file_type: 'fastq' or 'fasta' (optional, default: fastq)" << std::endl;
        std::cout << "  threads: Number of threads (optional, 0 for auto)" << std::endl;
        std::cout << "  min_score: Minimum alignment score (optional, default: 5)" << std::endl;
        std::cout << std::endl;
        std::cout << "Example: " << argv[0] << " reads.fq barcodes.fa 5 error_model.txt 100000 fastq 8 10" << std::endl;
        return 1;
    }
    
    std::string read_file = argv[1];
    std::string ref_fasta = argv[2];
    int k_size = std::stoi(argv[3]);
    std::string output_file = argv[4];
    
    int max_reads = -1;
    if (argc >= 6) {
        max_reads = std::stoi(argv[5]);
    }
    
    bool is_fastq = true;
    if (argc >= 7) {
        std::string file_type = argv[6];
        std::transform(file_type.begin(), file_type.end(), file_type.begin(), ::tolower);
        is_fastq = (file_type == "fastq");
    }
    
    int num_threads = 0;
    if (argc >= 8) {
        num_threads = std::stoi(argv[7]);
    }
    
    int min_alignment_score = 5;
    if (argc >= 9) {
        min_alignment_score = std::stoi(argv[8]);
    }
    
    process_barcode_dataset(read_file, ref_fasta, k_size, output_file, max_reads, is_fastq, num_threads, min_alignment_score);
    
    return 0;
}