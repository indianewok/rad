#include "include/rad/rad_headers.h"

// ============================================================================
// File I/O Utilities
// ============================================================================

struct GzWriter {
    gzFile fp = nullptr;
    
    explicit GzWriter(const std::string& path, int level = Z_DEFAULT_COMPRESSION) {
        fp = gzopen(path.c_str(), "wb");
        if (!fp) throw std::runtime_error("Failed to open gzip file: " + path);
        gzbuffer(fp, 1 << 20);
        gzsetparams(fp, level, Z_DEFAULT_STRATEGY);
    }
    
    void write(const char* data, size_t n) {
        if (!fp) throw std::runtime_error("gzFile is not open");
        int wrote = gzwrite(fp, data, static_cast<unsigned>(n));
        if (wrote == 0 || (size_t)wrote != n) 
            throw std::runtime_error("gzwrite failed");
    }
    
    void write(const std::string& s) { write(s.data(), s.size()); }
    
    ~GzWriter() { if (fp) gzclose(fp); }
};

static inline bool ends_with(const std::string& s, const char* suffix) {
    const size_t n = std::strlen(suffix);
    return s.size() >= n && std::equal(suffix, suffix + n, s.end() - n);
}

static std::vector<std::string> read_lines(const std::string& path) {
    std::vector<std::string> lines;
    
    if (ends_with(path, ".gz")) {
        gzFile f = gzopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("Failed to open: " + path);
        
        const int BUFSZ = 1 << 15;
        std::unique_ptr<char[]> buf(new char[BUFSZ]);
        std::string cur;
        cur.reserve(BUFSZ);
        
        int r;
        while ((r = gzread(f, buf.get(), BUFSZ)) > 0) {
            for (int i = 0; i < r; ++i) {
                char c = buf[i];
                if (c == '\n') {
                    lines.push_back(cur);
                    cur.clear();
                } else if (c != '\r') {
                    cur.push_back(c);
                }
            }
        }
        if (!cur.empty()) lines.push_back(cur);
        gzclose(f);
    } else {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Failed to open: " + path);
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            lines.push_back(std::move(line));
        }
    }
    return lines;
}

// ============================================================================
// Barcode Loading
// ============================================================================

static std::string trim(std::string s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    size_t e = s.find_last_not_of(" \t\r\n");
    return (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
}

// Extract first column from CSV/TSV, skip header if it contains "barcode"
static std::vector<std::string> parse_first_column(const std::vector<std::string>& lines) {
    if (lines.empty()) return {};
    
    size_t start_idx = 0;
    if (lines[0].find("barcode") != std::string::npos || 
        lines[0].find("Barcode") != std::string::npos) {
        start_idx = 1;
    }
    
    std::vector<std::string> out;
    out.reserve(lines.size() - start_idx);
    
    for (size_t i = start_idx; i < lines.size(); ++i) {
        const std::string& line = lines[i];
        if (line.empty()) continue;
        
        size_t cut = std::min(line.find(','), line.find('\t'));
        if (cut == std::string::npos) cut = line.size();
        
        std::string field = trim(line.substr(0, cut));
        if (!field.empty()) out.push_back(std::move(field));
    }
    return out;
}

// Parse "path" or "path:N" format
static std::pair<std::string, size_t> parse_path_with_limit(const std::string& spec) {
    std::string s = trim(spec);
    if (s.empty()) return {"", 0};
    
    size_t colon = s.rfind(':');
    if (colon == std::string::npos) return {s, 0};
    
    std::string path = trim(s.substr(0, colon));
    std::string tail = trim(s.substr(colon + 1));
    
    if (!tail.empty() && std::all_of(tail.begin(), tail.end(), ::isdigit)) {
        return {path, std::stoull(tail)};
    }
    return {s, 0};
}

static std::vector<std::string> load_barcodes(const std::string& spec, bool verbose, 
                                               std::mt19937* rng, bool random_sample) {
    auto [path, limit] = parse_path_with_limit(spec);
    if (path.empty()) return {};
    
    auto barcodes = parse_first_column(read_lines(path));
    
    if (limit > 0 && barcodes.size() > limit) {
        if (random_sample && rng) {
            std::shuffle(barcodes.begin(), barcodes.end(), *rng);
        }
        barcodes.resize(limit);
    }
    
    if (verbose) {
        std::cout << "[synth_gen] Loaded " << barcodes.size() 
                  << " barcodes from " << path;
        if (limit > 0) std::cout << " (limited to " << limit << ")";
        std::cout << "\n";
    }
    
    return barcodes;
}

// ============================================================================
// Synthetic Read Generator
// ============================================================================

class SyntheticGenerator {
private:
    std::mt19937 rng_;
    std::uniform_int_distribution<int> base_dist_;
    std::uniform_int_distribution<int> read_len_dist_;
    const std::vector<char> bases_ = {'A', 'C', 'G', 'T'};
    
    // Strand orientation settings
    bool rc_all_ = false;        // Force all reads to RC
    double rc_frac_ = 0.5;       // Default 50% RC rate
    
public:
    SyntheticGenerator(int seed = 42, int read_min = 300, int read_max = 600)
        : rng_(seed)
        , base_dist_(0, 3)
        , read_len_dist_(read_min, read_max) 
    {
        if (read_min <= 0 || read_max < read_min) {
            throw std::runtime_error("Invalid read length range");
        }
    }
    
    void set_read_length_range(int min_len, int max_len) {
        if (min_len <= 0 || max_len < min_len) {
            throw std::runtime_error("Invalid read length range");
        }
        read_len_dist_ = std::uniform_int_distribution<int>(min_len, max_len);
    }
    
    void set_rc_all(bool v) { rc_all_ = v; }
    void set_rc_frac(double f) { rc_frac_ = std::clamp(f, 0.0, 1.0); }
    
    // Generate in "reads mode" - sample barcodes independently per read
    void generate_reads_mode(const ReadLayout& layout,
                            const std::unordered_map<std::string, std::vector<std::string>>& barcode_lists,
                            const std::string& output_prefix,
                            bool output_fastq,
                            uint64_t num_reads,
                            bool verbose) {
        if (num_reads == 0) {
            throw std::runtime_error("num_reads must be > 0");
        }
        
        const std::string ext = output_fastq ? ".fq.gz" : ".fa.gz";
        GzWriter out_reads(output_prefix + ext);
        std::ofstream out_meta(output_prefix + "_metadata.csv");
        if (!out_meta) throw std::runtime_error("Failed to create metadata CSV");
        
        out_meta << "read_id,original_barcode,original_umi\n";
        
        std::bernoulli_distribution rc_flip(rc_frac_);
        std::string record, qualities;
        record.reserve(1024);
        qualities.reserve(512);
        
        auto t0 = std::chrono::steady_clock::now();
        
        for (uint64_t i = 0; i < num_reads; ++i) {
            // Assemble read in forward orientation
            auto [seq, bc_str, umi] = assemble_read(layout, barcode_lists);
            
            // Apply strand orientation
            bool flip = rc_all_ || (rc_frac_ > 0.0 && rc_flip(rng_));
            if (flip) seq = seq_utils::revcomp(seq);
            
            // Write record
            const std::string read_id = "synthetic_" + std::to_string(i + 1);
            write_record(out_reads, read_id, seq, output_fastq, record, qualities);
            out_meta << read_id << "," << bc_str << "," << umi << "\n";
            
            // Progress logging
            if (((i + 1) % 10000) == 0) {
                log_progress(i + 1, t0);
            }
        }
        
        if (verbose) {
            std::cout << "[synth_gen] Generated " << num_reads << " reads\n";
        }
    }
    
    // Generate in "cells mode" - sample barcodes once per cell, multiple reads per cell
    void generate_cells_mode(const ReadLayout& layout,
                            const std::unordered_map<std::string, std::vector<std::string>>& barcode_lists,
                            const std::string& output_prefix,
                            bool output_fastq,
                            uint64_t num_cells,
                            uint64_t reads_per_cell_fixed,
                            std::optional<std::pair<uint64_t,uint64_t>> reads_per_cell_range,
                            bool verbose) {
        if (num_cells == 0) {
            throw std::runtime_error("num_cells must be > 0");
        }
        
        const std::string ext = output_fastq ? ".fq.gz" : ".fa.gz";
        GzWriter out_reads(output_prefix + ext);
        std::ofstream out_meta(output_prefix + "_metadata.csv");
        if (!out_meta) throw std::runtime_error("Failed to create metadata CSV");
        
        out_meta << "read_id,cell_id,original_barcode,original_umi\n";
        
        std::bernoulli_distribution rc_flip(rc_frac_);
        
        // Setup reads-per-cell distribution
        std::uniform_int_distribution<uint64_t> rpc_dist(1, 1);
        bool use_range = (reads_per_cell_fixed == 0 && reads_per_cell_range.has_value());
        
        if (use_range) {
            auto [a, b] = *reads_per_cell_range;
            if (a == 0 || b < a) {
                throw std::runtime_error("reads_per_cell_range must be A:B with 1 <= A <= B");
            }
            rpc_dist = std::uniform_int_distribution<uint64_t>(a, b);
        }
        
        std::string record, qualities;
        record.reserve(1024);
        qualities.reserve(512);
        
        auto t0 = std::chrono::steady_clock::now();
        uint64_t total_reads = 0;
        
        for (uint64_t c = 0; c < num_cells; ++c) {
            // Sample cell barcode once (for all forward-oriented barcode elements)
            std::unordered_map<std::string, std::string> cell_barcodes = sample_cell_barcodes(layout, barcode_lists);
            
            // Determine reads for this cell
            uint64_t rpc = reads_per_cell_fixed > 0 ? reads_per_cell_fixed : rpc_dist(rng_);
            
            for (uint64_t r = 0; r < rpc; ++r) {
                // Assemble read with cell barcodes
                auto [seq, bc_str, umi] = assemble_read_with_cell_barcodes(layout, cell_barcodes);
                
                // Apply strand orientation
                bool flip = rc_all_ || (rc_frac_ > 0.0 && rc_flip(rng_));
                if (flip) seq = seq_utils::revcomp(seq);
                
                // Write record
                const std::string read_id = "synthetic_c" + std::to_string(c + 1) + 
                                           "_r" + std::to_string(r + 1);
                write_record(out_reads, read_id, seq, output_fastq, record, qualities);
                out_meta << read_id << ",c" << (c + 1) << "," << bc_str << "," << umi << "\n";
                
                ++total_reads;
                if ((total_reads % 10000) == 0) {
                    log_progress(total_reads, t0);
                }
            }
        }
        
        if (verbose) {
            std::cout << "[synth_gen] Generated " << num_cells << " cells, " 
                      << total_reads << " total reads\n";
        }
    }
    
private:
    std::string generate_random_sequence(int length) {
        std::string s;
        s.reserve(length);
        for (int i = 0; i < length; ++i) {
            s.push_back(bases_[base_dist_(rng_)]);
        }
        return s;
    }
    
    // Sample cell barcodes once per cell (only for forward barcode elements)
    std::unordered_map<std::string, std::string> 
    sample_cell_barcodes(const ReadLayout& layout,
                        const std::unordered_map<std::string, std::vector<std::string>>& barcode_lists) {
        std::unordered_map<std::string, std::string> cell_barcodes;
        
        for (const auto& elem : layout.by_order()) {
            if (elem.global_class != "barcode") continue;
            
            // Skip reverse elements - only sample forward barcodes
            if (elem.direction != "forward") continue;
            
            // Sample barcode for this element
            auto it = barcode_lists.find(elem.class_id);
            if (it != barcode_lists.end() && !it->second.empty()) {
                std::uniform_int_distribution<size_t> dist(0, it->second.size() - 1);
                cell_barcodes[elem.class_id] = it->second[dist(rng_)];
            } else {
                // Generate random barcode
                int len = elem.expected_length.value_or(16);
                cell_barcodes[elem.class_id] = generate_random_sequence(len);
            }
        }
        
        return cell_barcodes;
    }
    
    std::tuple<std::string, std::string, std::string>
    assemble_read(const ReadLayout& layout,
                  const std::unordered_map<std::string, std::vector<std::string>>& barcode_lists) {
        
        std::string full_seq, bc_concat, umi;
        bool first_bc = true;
        const int read_length = read_len_dist_(rng_);
        
        for (const auto& elem : layout.by_order()) {
            if (elem.global_class == "start" || elem.global_class == "stop") {
                continue;
            }
            
            // Only process forward elements
            if (elem.direction != "forward") continue;
            
            // Static elements
            if (elem.type == "static") {
                if (elem.global_class == "poly_tail") {
                    int L = elem.expected_length.value_or(16);
                    char base = (elem.seq.find('A') != std::string::npos || 
                               elem.class_id.find("poly_a") != std::string::npos) ? 'A' : 'T';
                    full_seq.append(L, base);
                } else {
                    full_seq.append(elem.seq);
                }
                continue;
            }
            
            // Barcode elements
            if (elem.global_class == "barcode") {
                std::string barcode;
                
                // Get barcode from list or generate random
                auto it = barcode_lists.find(elem.class_id);
                if (it != barcode_lists.end() && !it->second.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, it->second.size() - 1);
                    barcode = it->second[dist(rng_)];
                } else {
                    int len = elem.expected_length.value_or(16);
                    barcode = generate_random_sequence(len);
                }
                
                // Add barcode to BOTH sequence AND metadata
                full_seq.append(barcode);
                
                if (!barcode.empty()) {
                    if (!first_bc) bc_concat.push_back('+');
                    bc_concat.append(barcode);
                    first_bc = false;
                }
                continue;
            }
            
            // UMI elements
            if (elem.global_class == "umi") {
                int L = elem.expected_length.value_or(10);
                std::string u = generate_random_sequence(L);
                full_seq.append(u);
                if (umi.empty()) umi = u;
                continue;
            }
            
            // Read elements
            if (elem.global_class == "read") {
                full_seq.append(generate_random_sequence(read_length));
                continue;
            }
            
            // Default: random sequence
            full_seq.append(generate_random_sequence(elem.expected_length.value_or(20)));
        }
        
        return {full_seq, bc_concat, umi};
    }
    
    std::tuple<std::string, std::string, std::string>
    assemble_read_with_cell_barcodes(const ReadLayout& layout,
                                     const std::unordered_map<std::string, std::string>& cell_barcodes) {
        
        std::string full_seq, bc_concat, umi;
        bool first_bc = true;
        const int read_length = read_len_dist_(rng_);
        
        for (const auto& elem : layout.by_order()) {
            if (elem.global_class == "start" || elem.global_class == "stop") {
                continue;
            }
            
            // Only process forward elements
            if (elem.direction != "forward") continue;
            
            // Static elements
            if (elem.type == "static") {
                if (elem.global_class == "poly_tail") {
                    int L = elem.expected_length.value_or(16);
                    char base = (elem.seq.find('A') != std::string::npos || 
                               elem.class_id.find("poly_a") != std::string::npos) ? 'A' : 'T';
                    full_seq.append(L, base);
                } else {
                    full_seq.append(elem.seq);
                }
                continue;
            }
            
            // Barcode elements
            if (elem.global_class == "barcode") {
                // Get cell barcode
                auto it = cell_barcodes.find(elem.class_id);
                if (it == cell_barcodes.end()) {
                    throw std::runtime_error("Cell barcode not found for: " + elem.class_id);
                }
                
                std::string barcode = it->second;
                
                // Add barcode to BOTH sequence AND metadata
                full_seq.append(barcode);
                
                if (!barcode.empty()) {
                    if (!first_bc) bc_concat.push_back('+');
                    bc_concat.append(barcode);
                    first_bc = false;
                }
                continue;
            }
            
            // UMI elements (unique per read)
            if (elem.global_class == "umi") {
                int L = elem.expected_length.value_or(10);
                std::string u = generate_random_sequence(L);
                full_seq.append(u);
                if (umi.empty()) umi = u;
                continue;
            }
            
            // Read elements
            if (elem.global_class == "read") {
                full_seq.append(generate_random_sequence(read_length));
                continue;
            }
            
            // Default: random sequence
            full_seq.append(generate_random_sequence(elem.expected_length.value_or(20)));
        }
        
        return {full_seq, bc_concat, umi};
    }
    
    void write_record(GzWriter& writer, 
                     const std::string& read_id,
                     const std::string& seq,
                     bool is_fastq,
                     std::string& record_buf,
                     std::string& qual_buf) {
        record_buf.clear();
        
        if (is_fastq) {
            qual_buf.assign(seq.size(), '5');  // Dummy quality score
            record_buf.append("@").append(read_id).append("\n");
            record_buf.append(seq).append("\n+\n");
            record_buf.append(qual_buf).append("\n");
        } else {
            record_buf.append(">").append(read_id).append("\n");
            record_buf.append(seq).append("\n");
        }
        
        writer.write(record_buf);
    }
    
    void log_progress(uint64_t count, std::chrono::steady_clock::time_point t0) {
        auto t1 = std::chrono::steady_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();
        double rate = count / std::max(1e-9, sec);
        
        std::cerr << "[synth_gen] " << count << " reads  |  "
                  << std::fixed << std::setprecision(2) << sec << " s  |  "
                  << std::fixed << std::setprecision(1) << rate << " reads/s\n";
    }
};

// ============================================================================
// Main Driver
// ============================================================================

void run_synthetic_generator(const std::string& layout_csv_path,
                             const std::string& output_prefix,
                             bool output_fastq,
                             bool verbose,
                             int seed,
                             int read_min,
                             int read_max,
                             bool rc_all,
                             double rc_frac,
                             // Reads mode
                             uint64_t num_reads,
                             // Cells mode
                             uint64_t num_cells,
                             uint64_t reads_per_cell_fixed,
                             std::optional<std::pair<uint64_t,uint64_t>> reads_per_cell_range) {
    
    if (verbose) {
        std::cout << "[synth_gen] Layout: " << layout_csv_path << "\n";
        std::cout << "[synth_gen] Output: " << output_prefix 
                  << (output_fastq ? ".fq.gz" : ".fa.gz") << "\n";
        std::cout << "[synth_gen] Read length range: [" << read_min << ", " << read_max << "]\n";
        
        if (rc_all) {
            std::cout << "[synth_gen] Orientation: ALL reverse complement\n";
        } else {
            std::cout << "[synth_gen] Orientation: " << (rc_frac * 100.0) 
                      << "% reverse complement\n";
        }
    }
    
    // Load layout
    ReadLayout layout;
    layout.prep_new_layout(layout_csv_path, verbose);
    
    // Load barcode whitelists for forward barcode elements only
    std::unordered_map<std::string, std::vector<std::string>> barcode_lists;
    std::mt19937 rng_for_loading(seed);
    
    for (const auto& elem : layout.by_order()) {
        if (elem.global_class != "barcode") continue;
        
        // Only load for forward elements
        if (elem.direction != "forward") continue;
        
        if (!elem.whitelist_path.empty()) {
            barcode_lists[elem.class_id] = load_barcodes(elem.whitelist_path, verbose, 
                                                         &rng_for_loading, false);
        }
    }
    
    // Create generator
    SyntheticGenerator gen(seed, read_min, read_max);
    gen.set_rc_all(rc_all);
    gen.set_rc_frac(rc_frac);
    
    // Choose mode
    if (num_cells > 0) {
        gen.generate_cells_mode(layout, barcode_lists, output_prefix, output_fastq,
                               num_cells, reads_per_cell_fixed, reads_per_cell_range, 
                               verbose);
    } else {
        if (num_reads == 0) {
            throw std::runtime_error("Either --num-reads or --cells must be specified");
        }
        gen.generate_reads_mode(layout, barcode_lists, output_prefix, output_fastq,
                               num_reads, verbose);
    }
}

// ============================================================================
// CLI
// ============================================================================

static void print_usage(const char* prog) {
    std::cerr
      << "Usage: " << prog << " [options]\n\n"
      << "Generate synthetic reads from a read layout specification.\n\n"
      << "Required:\n"
      << "  -l, --layout FILE           Read layout CSV file\n"
      << "  -o, --output PREFIX         Output file prefix\n\n"
      << "Generation Mode (choose one):\n"
      << "  -n, --num-reads N           Generate N reads (reads mode)\n"
      << "      --cells C               Generate C cells (cells mode)\n\n"
      << "Cells Mode Options (requires --cells):\n"
      << "      --reads-per-cell N      Fixed reads per cell\n"
      << "      --reads-per-cell-range A:B  Random reads per cell in [A,B]\n\n"
      << "Common Options:\n"
      << "  -s, --seed SEED             Random seed (default: 42)\n"
      << "  -f, --format FMT            Output format: fastq|fasta (default: fastq)\n"
      << "      --read-min L            Min read length (default: 300)\n"
      << "      --read-max U            Max read length (default: 600)\n"
      << "      --rc-reads              Force ALL reads to reverse complement\n"
      << "      --rc-frac P             Fraction to RC, 0-1 (default: 0.5)\n"
      << "  -v, --verbose               Verbose output\n"
      << "  -h, --help                  Show this help\n\n"
      << "Examples:\n"
      << "  # Generate 10k reads\n"
      << "  " << prog << " -l layout.csv -o out -n 10000\n\n"
      << "  # Generate 1k cells with 50-100 reads each\n"
      << "  " << prog << " -l layout.csv -o out --cells 1000 --reads-per-cell-range 50:100\n";
}

int main_synth_gen(int argc, char* argv[]) {
    // Default parameters
    std::string layout_path, output_prefix, format = "fastq";
    int seed = 42;
    int read_min = 300, read_max = 600;
    bool verbose = false;
    bool rc_all = false;
    double rc_frac = 0.5;
    
    // Mode selection
    uint64_t num_reads = 0;
    uint64_t num_cells = 0;
    uint64_t reads_per_cell_fixed = 0;
    std::optional<std::pair<uint64_t,uint64_t>> reads_per_cell_range;
    
    // Parse command line
    const char* optstring = "l:o:s:f:n:vh";
    struct option longopts[] = {
        {"layout",               required_argument, nullptr, 'l'},
        {"output",               required_argument, nullptr, 'o'},
        {"seed",                 required_argument, nullptr, 's'},
        {"format",               required_argument, nullptr, 'f'},
        {"num-reads",            required_argument, nullptr, 'n'},
        {"read-min",             required_argument, nullptr, 1000},
        {"read-max",             required_argument, nullptr, 1001},
        {"rc-reads",             no_argument,       nullptr, 1002},
        {"rc-frac",              required_argument, nullptr, 1003},
        {"cells",                required_argument, nullptr, 1004},
        {"reads-per-cell",       required_argument, nullptr, 1005},
        {"reads-per-cell-range", required_argument, nullptr, 1006},
        {"verbose",              no_argument,       nullptr, 'v'},
        {"help",                 no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'l': layout_path = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 's': seed = std::stoi(optarg); break;
            case 'f': format = optarg; break;
            case 'n': num_reads = std::stoull(optarg); break;
            case 1000: read_min = std::stoi(optarg); break;
            case 1001: read_max = std::stoi(optarg); break;
            case 1002: rc_all = true; break;
            case 1003: rc_frac = std::stod(optarg); break;
            case 1004: num_cells = std::stoull(optarg); break;
            case 1005: reads_per_cell_fixed = std::stoull(optarg); break;
            case 1006: {
                std::string s = optarg;
                auto pos = s.find(':');
                if (pos == std::string::npos) {
                    std::cerr << "Error: --reads-per-cell-range must be A:B\n";
                    return 1;
                }
                uint64_t a = std::stoull(s.substr(0, pos));
                uint64_t b = std::stoull(s.substr(pos + 1));
                reads_per_cell_range = {a, b};
                break;
            }
            case 'v': verbose = true; break;
            case 'h': print_usage(argv[0]); return 0;
            default: std::cerr << "Use -h for help\n"; return 1;
        }
    }
    
    // Validate required arguments
    if (layout_path.empty() || output_prefix.empty()) {
        std::cerr << "Error: --layout and --output are required\n";
        print_usage(argv[0]);
        return 1;
    }
    
    // Validate read length
    if (read_min <= 0 || read_max < read_min) {
        std::cerr << "Error: Invalid read length range\n";
        return 1;
    }
    
    // Validate RC fraction
    if (rc_frac < 0.0 || rc_frac > 1.0) {
        std::cerr << "Error: --rc-frac must be in [0.0, 1.0]\n";
        return 1;
    }
    
    // Validate format
    bool output_fastq = true;
    if (format == "fasta" || format == "fa") {
        output_fastq = false;
    } else if (format != "fastq" && format != "fq") {
        std::cerr << "Error: format must be fastq or fasta\n";
        return 1;
    }
    
    // Validate mode selection
    if (num_cells > 0) {
        // Cells mode
        if (reads_per_cell_fixed == 0 && !reads_per_cell_range.has_value()) {
            std::cerr << "Error: --cells requires --reads-per-cell or --reads-per-cell-range\n";
            return 1;
        }
        if (reads_per_cell_fixed > 0 && reads_per_cell_range.has_value()) {
            std::cerr << "Error: Use either --reads-per-cell OR --reads-per-cell-range\n";
            return 1;
        }
    } else {
        // Reads mode
        if (num_reads == 0) {
            std::cerr << "Error: Specify --num-reads or --cells\n";
            return 1;
        }
    }
    
    // Run generator
    try {
        run_synthetic_generator(layout_path, output_prefix, output_fastq, verbose,
                               seed, read_min, read_max, rc_all, rc_frac,
                               num_reads, num_cells, reads_per_cell_fixed, 
                               reads_per_cell_range);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
}