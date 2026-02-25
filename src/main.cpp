// main.cpp
#include "include/rad/rad_headers.h"

// ============================================================================
// USAGE FUNCTIONS
// ============================================================================

static void print_main_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <command> [options]\n\n"
              << "Commands:\n"
              << "  prep              Prepare read layout and/or position map\n"
              << "  demux             Run full demultiplexing pipeline\n"
              << "  reformat          Reformat fastx.gz outputs for downstream applications\n"
              << "  list              List available read layouts and/or whitelist kits\n"
              << "  modify            Add/remove read layout or whitelist-kit mappings\n"
              << "  help              Show help for a specific command\n\n"
              << "Use '" << prog << " <command> --help' for command-specific options\n";
}

static void usage_demux(const char* prog) {
    std::cerr
      << "Usage: rad " << prog << " -l LAYOUT -q FASTQ [options]\n\n"
      << "Required:\n"
      << "  -l, --layout                      layout key (five_prime, three_prime, splitseq) or path\n"
      << "  -q, --fastq                       input FASTQ file\n\n"
      << "Optional:\n"
      << "  -k, --kit                         use this kit's default whitelist\n"
      << "  -g, --global_whitelist            path to custom global whitelist CSV\n"
      << "  -c, --custom_whitelist            path to custom whitelist CSV\n"
      << "  -R, --bc_correction_mode          'offensive' (default) or 'defensive'\n"
      << "      --joint-bc-mode               'default' (permissive) or 'strict' (all barcodes must pass)\n"
      << "  -M, --whitelist_mutation          mutation space for whitelist (default: 2)\n"
      << "  -m, --generated_mutation          mutation space for generated barcodes (default: 2)\n"
      << "  -n, --max_reads                   maximum number of reads (default: all)\n"
      << "  -z, --chunk_size                  chunk size for processing (default: 5000)\n"
      << "  -o, --output                      filename prefix (default: output)\n"
      << "  -d, --dir                         output directory (default: current)\n"
      << "  -F, --log-file                    save output to log file\n"
      << "  -w, --write_dbg                   write debug files (.sig, .csv, metrics)\n"
      << "  -b, --bc_split                    split reads into per-barcode FASTQs\n"
      << "  -t, --threads                     number of threads (default: 1)\n"
      << "  -v, --verbose                     verbose mode\n"
      << "  -D, --max_verbose                 maximum verbosity (debug)\n"
      << "  -h, --help                        show this help\n";
}

static void usage_reformat(const char* prog) {
    std::cerr
      << "Usage: rad " << prog << " -q INPUT [options]\n\n"
      << "Options:\n"
      << "  -q, --fastq           input FASTQ/FASTA (.fq/.fa/.gz)\n"
      << "  -o, --outdir          output directory for per-barcode fastqs (required if --split-bc)\n"
      << "      --split-bc        split reads into per-barcode .fq.gz files (CB:Z tag)\n"
      << "      --reformat-header collapse headers to QNAME_CB_UB (default underscore)\n"
      << "      --reformat-delim  delimiter for collapsed headers (single char; use \" \" for space)\n"
      << "      --coordinate[=MODE] rewrite headers using coordinate mode (default MODE: vizHD-v1)\n"
      << "      --bin-size INT    target bin size in microns for coordinate mode (default: 2)\n"
      << "  -t, --threads         worker threads (default: 2)\n"
      << "  -v, --verbose         verbose logging\n"
      << "  -h, --help            show this help\n\n"
      << "If both --split-bc and --reformat-header are set, headers are collapsed before splitting.\n"
      << "If only --reformat-header is set, the input file is streamed and rewritten in place.\n";
}

static void usage_prep(const char* prog) {
    std::cerr
      << "Usage: rad " << prog << " -l LAYOUT [options]\n\n"
      << "Generate read layout structure and/or position map with misalignment statistics.\n\n"
      << "Required:\n"
      << "  -l, --layout                      layout key or path to layout CSV\n\n"
      << "Mode flags (at least one required):\n"
      << "      --read-layout                 generate and display read layout structure\n"
      << "      --position-map                generate position map with misalignment stats\n"
      << "                                    (requires -q/--fastq)\n\n"
      << "Optional:\n"
      << "  -q, --fastq                       input FASTQ file (required for --position-map)\n"
      << "  -o, --output                      output base path (saves _layout.csv, _position_map.csv)\n"
      << "  -n, --max_reads                   max reads for misalignment sampling (default: 50000)\n"
      << "  -t, --threads                     number of threads (default: 1)\n"
      << "  -v, --verbose                     verbose mode\n"
      << "  -D, --max_verbose                 maximum verbosity (debug)\n"
      << "  -h, --help                        show this help\n\n"
      << "Examples:\n"
      << "  rad " << prog << " -l five_prime --read-layout\n"
      << "  rad " << prog << " -l my_layout.csv --position-map -q reads.fq.gz -o output\n"
      << "  rad " << prog << " -l five_prime --read-layout --position-map -q reads.fq.gz -o output\n";
}

static void usage_list(const char* prog) {
    std::cerr
      << "Usage: rad " << prog << " [--layouts] [--whitelists]\n\n"
      << "List registered layout keys and whitelist kit keys with resolved paths.\n\n"
      << "Options:\n"
      << "      --layouts                    list read layout keys and paths\n"
      << "      --whitelists                 list whitelist kit keys and paths\n"
      << "  -h, --help                       show this help\n\n"
      << "If neither --layouts nor --whitelists is provided, both are listed.\n";
}

static void usage_modify(const char* prog) {
    std::cerr
      << "Usage: rad " << prog
      << " (--layout KEY | --whitelist KIT) (--add PATH | --remove)\n\n"
      << "Modify registered mappings for layout keys or whitelist kits.\n\n"
      << "Options:\n"
      << "      --layout KEY                 target read layout key\n"
      << "      --whitelist KIT              target whitelist kit key\n"
      << "      --add PATH                   add or update mapping to PATH\n"
      << "      --remove                     remove mapping\n"
      << "  -h, --help                       show this help\n";
}

// =============================================================================
// HELPER: tee_buff for logging to both console and file, JIC want -D output
// =============================================================================

class tee_buff : public std::streambuf {
private:
    std::streambuf* sb1;
    std::streambuf* sb2;
    
public:
    tee_buff(std::streambuf* sb1, std::streambuf* sb2) 
        : sb1(sb1), sb2(sb2) {}
    
protected:
    virtual int overflow(int c) override {
        if (c == EOF) {
            return !EOF;
        } else {
            int const r1 = sb1->sputc(c);
            int const r2 = sb2->sputc(c);
            return (r1 == EOF || r2 == EOF) ? EOF : c;
        }
    }
    
    virtual int sync() override {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }
};

// ============================================================================
// COMMAND: prep (combined layout + position_map)
// ============================================================================

int cmd_prep(int argc, char* argv[]) {
    std::string layout_key, fastq_path, output_base;
    bool verbose = false, max_verbose = false;
    bool do_read_layout = false, do_position_map = false;
    int nthreads = 1;
    size_t max_reads = 50000;
    
    const char* optstring = "l:q:o:n:t:vDh";
    struct option longopts[] = {
        {"layout",       required_argument, nullptr, 'l'},
        {"fastq",        required_argument, nullptr, 'q'},
        {"output",       required_argument, nullptr, 'o'},
        {"max_reads",    required_argument, nullptr, 'n'},
        {"threads",      required_argument, nullptr, 't'},
        {"verbose",      no_argument,       nullptr, 'v'},
        {"max_verbose",  no_argument,       nullptr, 'D'},
        {"help",         no_argument,       nullptr, 'h'},
        {"read-layout",  no_argument,       nullptr,  1 },
        {"position-map", no_argument,       nullptr,  2 },
        {nullptr, 0, nullptr, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'l': layout_key = optarg; break;
            case 'q': fastq_path = optarg; break;
            case 'o': output_base = optarg; break;
            case 'n': max_reads = std::stoull(optarg); break;
            case 't': nthreads = std::stoi(optarg); break;
            case 'v': verbose = true; break;
            case 'D': max_verbose = true; verbose = true; break;
            case 'h': usage_prep(argv[0]); return 0;
            case 1:   do_read_layout = true; break;
            case 2:   do_position_map = true; break;
            default:  usage_prep(argv[0]); return 1;
        }
    }
    
    // Validate required args
    if (layout_key.empty()) {
        std::cerr << "[ERROR] --layout is required\n\n";
        usage_prep(argv[0]);
        return 1;
    }
    
    if (!do_read_layout && !do_position_map) {
        std::cerr << "[ERROR] must specify --read-layout and/or --position-map\n\n";
        usage_prep(argv[0]);
        return 1;
    }
    
    if (do_position_map && fastq_path.empty()) {
        std::cerr << "[ERROR] --fastq is required when using --position-map\n\n";
        usage_prep(argv[0]);
        return 1;
    }
    
    if (do_position_map && output_base.empty()) {
        std::cerr << "[ERROR] --output is required when using --position-map\n\n";
        usage_prep(argv[0]);
        return 1;
    }
    
    try {
        auto start = std::chrono::steady_clock::now();
        
        ReadLayout read_layout;
        std::string layout_csv;
        
        // Determine if it's a built-in layout or custom path
        if (!config_utils::check_if_custom_rl(layout_key)) {
            layout_csv = config_utils::get_read_layout(layout_key);
            if (verbose) std::cout << "[prep] Using built-in layout: " << layout_key << "\n";
        } else {
            layout_csv = layout_key;
            if (verbose) std::cout << "[prep] Using custom layout: " << layout_key << "\n";
        }
        
        // Prep layout
        if (verbose) std::cout << "[prep] Preparing layout...\n";
        read_layout.prep_new_layout(layout_csv, max_verbose);
        
        // Display layout if requested
        if (do_read_layout) {
            std::cout << "\n========================================\n";
            std::cout << "         LAYOUT STRUCTURE               \n";
            std::cout << "========================================\n\n";
            
            read_layout.display_read_layout();
            
            std::cout << "\n========================================\n";
            std::cout << "Layout contains " << read_layout.size() << " elements\n";
            std::cout << "========================================\n";
        }
        
        // Generate position map if requested
        if (do_position_map) {
            // Save initial layout
            read_layout.write_to_csv(output_base, "layout");
            
            // Run misalignment analysis
            if (verbose) std::cout << "\n[prep] Computing misalignment statistics from FASTQ...\n";
            Misalignment_Setup mis(read_layout);
            mis.generate_misalignment_data(fastq_path, read_layout, nthreads, max_reads);
            
            // Generate position mapping
            if (verbose) std::cout << "[prep] Generating position mapping...\n";
            read_layout.generate_position_mapping();
            
            // Display if max verbose
            if (verbose) {
                std::cout << "\n";
                read_layout.display_read_layout();
            }
            
            // Write updated files
            if (verbose) std::cout << "[prep] Writing output files...\n";
            read_layout.write_to_csv(output_base, "both");
            
            auto elapsed = std::chrono::steady_clock::now() - start;
            std::cout << "\n[prep] Complete! Time: "
                      << std::chrono::duration_cast<std::chrono::seconds>(elapsed).count()
                      << " seconds\n";
            
            std::cout << "\n[prep] Output files:\n"
                      << "  " << output_base << "_layout.csv\n"
                      << "  " << output_base << "_position_map.csv\n";
        } else if (!output_base.empty()) {
            // Just save layout if output specified but no position map
            read_layout.write_to_csv(output_base, "layout");
            std::cout << "\n[prep] Saved to: " << output_base << "_layout.csv\n";
        }
        
        return 0;
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << "\n";
        return 1;
    }
}

// ============================================================================
// COMMAND: reformat
// ============================================================================

static std::optional<std::string> extract_cb(const std::string& id, const std::string& comment) {
    const std::string tag = "CB:Z:";
    const std::string* src = &comment;
    auto pos = comment.find(tag);
    if (pos == std::string::npos) {
        pos = id.find(tag);
        if (pos == std::string::npos) return std::nullopt;
        src = &id;
    }
    pos += tag.size();
    auto end = src->find_first_of(" \t", pos);
    if (end == std::string::npos) {
        return src->substr(pos);
    }
    return src->substr(pos, end - pos);
}

static std::optional<std::string> extract_ub(const std::string& comment) {
    const std::string tag = "UB:Z:";
    auto pos = comment.find(tag);
    if (pos == std::string::npos) return std::nullopt;
    pos += tag.size();
    auto end = comment.find_first_of(" \t", pos);
    if (end == std::string::npos) {
        return comment.substr(pos);
    }
    return comment.substr(pos, end - pos);
}

static std::string collapse_id(const std::string& qname,
                               const std::optional<std::string>& cb, const std::optional<std::string>& ub,
                               char sep = '_') {
    std::string out = qname;
    if (cb && !cb->empty()) {
        out.push_back(sep);
        out.append(*cb);
    }
    if (ub && !ub->empty()) {
        out.push_back(sep);
        out.append(*ub);
    }
    return out;
}

static std::optional<std::tuple<std::string, std::string, std::string>>
parse_collapsed_id(const std::string& id, char sep) {
    // Parse from the right so QNAME can safely contain the delimiter.
    auto p2 = id.rfind(sep);
    if (p2 == std::string::npos) return std::nullopt;
    auto p1 = id.rfind(sep, p2 == 0 ? std::string::npos : p2 - 1);
    if (p1 == std::string::npos) return std::nullopt;

    std::string qname = id.substr(0, p1);
    std::string cb = id.substr(p1 + 1, p2 - p1 - 1);
    std::string ub = id.substr(p2 + 1);
    if (qname.empty() || cb.empty() || ub.empty()) return std::nullopt;
    return std::make_tuple(std::move(qname), std::move(cb), std::move(ub));
}

static std::string make_visium_hd_comment(const std::string& bc, const std::string& cb,
                                          const std::string& ub, const std::string& sample_tag) {
    std::string comment;
    comment.reserve(bc.size() + cb.size() + ub.size() + sample_tag.size() + 40);
    comment += "BC:Z:";
    comment += bc;
    comment += "\t";
    comment += "CB:Z:";
    comment += cb;
    comment += "\tUB:Z:";
    comment += ub;
    comment += "\tRG:Z:";
    comment += sample_tag;
    return comment;
}

static std::string to_lower_copy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::string zpad(int v, int width) {
    if (v < 0) v = 0;
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << v;
    return oss.str();
}

static bool parse_positive_int(const std::string& s, int& out) {
    if (s.empty()) return false;
    try {
        size_t consumed = 0;
        int value = std::stoi(s, &consumed);
        if (consumed != s.size() || value <= 0) return false;
        out = value;
        return true;
    } catch (...) {
        return false;
    }
}

static bool parse_legacy_bin_token(const std::string& token, int& out_bin_size_um) {
    // Accept either integer microns (e.g. "8") or legacy token "s008um".
    if (parse_positive_int(token, out_bin_size_um)) return true;
    if (token.size() >= 4 && (token[0] == 's' || token[0] == 'S')) {
        auto lower = to_lower_copy(token);
        if (lower.size() > 3 && lower.substr(lower.size() - 2) == "um") {
            std::string num = lower.substr(1, lower.size() - 3);
            return parse_positive_int(num, out_bin_size_um);
        }
    }
    return false;
}

static std::vector<std::string> split_string(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= s.size()) {
        size_t pos = s.find(sep, start);
        if (pos == std::string::npos) pos = s.size();
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
        if (pos == s.size()) break;
    }
    return out;
}

struct vizhd_axis_reference {
    std::unordered_map<std::string, int> x_index_by_bc;
    std::unordered_map<std::string, int> y_index_by_bc;
    std::vector<int> x_lengths;
    std::vector<int> y_lengths;
};

static std::vector<std::string> load_barcode_csv_first_column(const std::string& path) {
    auto lines = streaming_utils::import_text(path);
    std::vector<std::string> out;
    out.reserve(lines.size());
    bool first = true;
    for (auto& line : lines) {
        auto trimmed = seq_utils::trim(line);
        if (trimmed.empty()) continue;
        if (first) {
            first = false;
            if (to_lower_copy(trimmed) == "whitelist_bcs") continue;
        }
        auto comma = trimmed.find(',');
        if (comma != std::string::npos) trimmed = trimmed.substr(0, comma);
        trimmed = seq_utils::trim(trimmed);
        if (!trimmed.empty()) out.push_back(trimmed);
    }
    return out;
}

static vizhd_axis_reference load_vizhd_v1_axis_reference() {
    vizhd_axis_reference ref;

    const std::string x_path = whitelist_utils::kit_to_path("visium_hd_bc1");
    const std::string y_path = whitelist_utils::kit_to_path("visium_hd_bc2");
    auto x_bcs = load_barcode_csv_first_column(x_path);
    auto y_bcs = load_barcode_csv_first_column(y_path);

    ref.x_index_by_bc.reserve(x_bcs.size());
    ref.y_index_by_bc.reserve(y_bcs.size());

    for (size_t i = 0; i < x_bcs.size(); ++i) {
        ref.x_index_by_bc.emplace(x_bcs[i], static_cast<int>(i));
        ref.x_lengths.push_back(static_cast<int>(x_bcs[i].size()));
    }
    for (size_t i = 0; i < y_bcs.size(); ++i) {
        ref.y_index_by_bc.emplace(y_bcs[i], static_cast<int>(i));
        ref.y_lengths.push_back(static_cast<int>(y_bcs[i].size()));
    }

    std::sort(ref.x_lengths.begin(), ref.x_lengths.end());
    ref.x_lengths.erase(std::unique(ref.x_lengths.begin(), ref.x_lengths.end()), ref.x_lengths.end());
    std::sort(ref.y_lengths.begin(), ref.y_lengths.end());
    ref.y_lengths.erase(std::unique(ref.y_lengths.begin(), ref.y_lengths.end()), ref.y_lengths.end());
    return ref;
}

static std::optional<std::pair<std::string, std::string>>
resolve_xy_barcodes_from_cb(const std::string& cb_tag, const vizhd_axis_reference& ref) {
    if (cb_tag.empty()) return std::nullopt;

    auto direct = split_string(cb_tag, '-');
    if (direct.size() >= 2) {
        const std::string& a = direct[0];
        const std::string& b = direct[1];
        bool a_is_x = ref.x_index_by_bc.find(a) != ref.x_index_by_bc.end();
        bool b_is_y = ref.y_index_by_bc.find(b) != ref.y_index_by_bc.end();
        if (a_is_x && b_is_y) return std::make_pair(a, b);

        bool b_is_x = ref.x_index_by_bc.find(b) != ref.x_index_by_bc.end();
        bool a_is_y = ref.y_index_by_bc.find(a) != ref.y_index_by_bc.end();
        if (b_is_x && a_is_y) return std::make_pair(b, a);
    }

    // Fallback for concatenated CB tags with no delimiter.
    std::vector<std::pair<std::string, std::string>> candidates;
    for (int lx : ref.x_lengths) {
        for (int ly : ref.y_lengths) {
            if (static_cast<int>(cb_tag.size()) != lx + ly) continue;

            std::string xs = cb_tag.substr(0, static_cast<size_t>(lx));
            std::string ys = cb_tag.substr(static_cast<size_t>(lx), static_cast<size_t>(ly));
            if (ref.x_index_by_bc.find(xs) != ref.x_index_by_bc.end() &&
                ref.y_index_by_bc.find(ys) != ref.y_index_by_bc.end()) {
                candidates.emplace_back(xs, ys);
            }

            ys = cb_tag.substr(0, static_cast<size_t>(ly));
            xs = cb_tag.substr(static_cast<size_t>(ly), static_cast<size_t>(lx));
            if (ref.x_index_by_bc.find(xs) != ref.x_index_by_bc.end() &&
                ref.y_index_by_bc.find(ys) != ref.y_index_by_bc.end()) {
                candidates.emplace_back(xs, ys);
            }
        }
    }

    if (candidates.size() == 1) return candidates.front();
    return std::nullopt;
}

static int project_axis_for_bin_size(int coord_2um, int bin_size_um) {
    if (bin_size_um <= 0) return coord_2um;
    long long scaled = static_cast<long long>(coord_2um) * 2LL;
    return static_cast<int>(scaled / static_cast<long long>(bin_size_um));
}

static std::string make_bin_rg_tag(int bin_size_um) {
    std::ostringstream oss;
    oss << "s" << std::setw(3) << std::setfill('0') << bin_size_um << "um";
    return oss.str();
}

static std::string make_spatial_mapping_id(int bin_size_um, int x_bin, int y_bin) {
    std::ostringstream oss;
    oss << "s_"
        << std::setw(3) << std::setfill('0') << bin_size_um
        << "um_"
        << zpad(x_bin, 5)
        << "_"
        << zpad(y_bin, 5)
        << "-1";
    return oss.str();
}

static std::string make_coordinate_read_id(const std::string& base_id, int x_bin, int y_bin) {
    std::string id = base_id;
    id += "_";
    id += zpad(x_bin, 5);
    id += "_";
    id += zpad(y_bin, 5);
    return id;
}

int cmd_reformat(int argc, char* argv[]) {
    std::string fastq_path, outdir;
    int nthreads = 2;
    bool verbose = false;
    bool do_split = false;
    bool do_reformat = false;
    bool do_coordinate = false;
    std::string coordinate_mode = "vizHD-v1";
    int bin_size_um = 2;
    char reformat_sep = '_';

    const char* optstring = "q:o:t:vh";
    struct option longopts[] = {
        {"fastq",            required_argument, nullptr, 'q'},
        {"outdir",           required_argument, nullptr, 'o'},
        {"threads",          required_argument, nullptr, 't'},
        {"verbose",          no_argument,       nullptr, 'v'},
        {"help",             no_argument,       nullptr, 'h'},
        {"split-bc",         no_argument,       nullptr,  1 },
        {"reformat-header",  no_argument,       nullptr,  2 },
        {"reformat-delim",   required_argument, nullptr,  3 },
        {"coordinate",       optional_argument, nullptr,  4 },
        {"bin-size",         required_argument, nullptr,  5 },
        // Backward-compatible aliases:
        {"visium-hd-header", no_argument,       nullptr,  6 },
        {"visium-hd-sample", required_argument, nullptr,  7 },
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'q': fastq_path = optarg; break;
            case 'o': outdir = optarg; break;
            case 't': nthreads = std::max(1, std::stoi(optarg)); break;
            case 'v': verbose = true; break;
            case 1: do_split = true; break;
            case 2: do_reformat = true; break;
            case 3: {
                std::string arg = optarg ? optarg : "";
                if (arg == "\\t") {
                    reformat_sep = '\t';
                } else if (arg == "\\n") {
                    reformat_sep = '\n';
                } else if (arg == "\\r") {
                    reformat_sep = '\r';
                } else if (arg.size() == 1) {
                    reformat_sep = arg[0];
                } else {
                    std::cerr << "[ERROR] --reformat-delim must be a single character (or \\t/\\n/\\r)\n\n";
                    usage_reformat(argv[0]);
                    return 1;
                }
                do_reformat = true;
                break;
            }
            case 4: {
                do_coordinate = true;
                do_reformat = true;
                if (optarg && *optarg) {
                    coordinate_mode = optarg;
                } else if (optind < argc && argv[optind] && argv[optind][0] != '-') {
                    coordinate_mode = argv[optind];
                    ++optind;
                }
                break;
            }
            case 5: {
                int parsed = 0;
                if (!parse_positive_int(optarg ? std::string(optarg) : "", parsed)) {
                    std::cerr << "[ERROR] --bin-size must be a positive integer (microns)\n\n";
                    usage_reformat(argv[0]);
                    return 1;
                }
                bin_size_um = parsed;
                do_coordinate = true;
                do_reformat = true;
                break;
            }
            case 6: {
                do_coordinate = true;
                do_reformat = true;
                if (coordinate_mode.empty()) coordinate_mode = "vizHD-v1";
                break;
            }
            case 7: {
                int parsed = 0;
                std::string arg = optarg ? optarg : "";
                if (!parse_legacy_bin_token(arg, parsed)) {
                    std::cerr << "[ERROR] --visium-hd-sample expects legacy sXXXum token or integer bin size\n\n";
                    usage_reformat(argv[0]);
                    return 1;
                }
                bin_size_um = parsed;
                do_coordinate = true;
                do_reformat = true;
                break;
            }
            case 'h': usage_reformat(argv[0]); return 0;
            default:  usage_reformat(argv[0]); return 1;
        }
    }

    if (!do_split && !do_reformat) {
        std::cerr << "[ERROR] must specify --split-bc and/or --reformat-header\n\n";
        usage_reformat(argv[0]);
        return 1;
    }
    if (fastq_path.empty()) {
        std::cerr << "[ERROR] --fastq is required\n\n";
        usage_reformat(argv[0]);
        return 1;
    }
    if (do_split && outdir.empty()) {
        std::cerr << "[ERROR] --outdir is required when --split-bc is set\n\n";
        usage_reformat(argv[0]);
        return 1;
    }

    boost::filesystem::path outdir_path(outdir);
    if (do_split && !boost::filesystem::exists(outdir_path)) {
        boost::filesystem::create_directories(outdir_path);
    }

    std::optional<vizhd_axis_reference> spatial_ref;
    if (do_coordinate) {
        std::string mode_lc = to_lower_copy(coordinate_mode);
        if (!(mode_lc == "vizhd-v1" || mode_lc == "visium-hd-v1" || mode_lc == "visium_hd_v1")) {
            std::cerr << "[ERROR] Unsupported --coordinate mode: " << coordinate_mode
                      << " (supported: vizHD-v1)\n";
            return 1;
        }
        try {
            spatial_ref = load_vizhd_v1_axis_reference();
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Failed to load coordinate reference for mode " << coordinate_mode
                      << ": " << e.what() << "\n";
            return 1;
        }
        if (verbose) {
            std::cout << "[reformat] coordinate mode: " << coordinate_mode
                      << ", bin-size: " << bin_size_um << "um"
                      << ", x barcodes: " << spatial_ref->x_index_by_bc.size()
                      << ", y barcodes: " << spatial_ref->y_index_by_bc.size() << "\n";
        }
    }

    struct ref_writer {
        std::mutex mtx;
        FILE* pigz_fp = nullptr;
        gzFile gz = nullptr;
        std::ofstream out;
        bool using_pigz = false;
        bool use_gz = false;
    };
    std::unordered_map<std::string, std::shared_ptr<ref_writer>> writers;
    std::mutex writers_mtx;
    pigz_writing pigz_out;

    auto get_writer = [&](const std::string& bc) -> ref_writer& {
        std::lock_guard<std::mutex> lock(writers_mtx);
        auto it = writers.find(bc);
        if (it != writers.end()) return *(it->second);

        boost::filesystem::path p = outdir_path / (bc + ".fq.gz");
        auto w = std::make_shared<ref_writer>();
        w->pigz_fp = pigz_out.open_pigz_pipe(p.string(), nthreads);
        if (w->pigz_fp) {
            w->using_pigz = true;
        } else {
            w->gz = gzopen(p.string().c_str(), "wb");
            if (!w->gz) throw std::runtime_error("Failed to open output " + p.string());
            gzbuffer(w->gz, 1 << 20); // 1 MiB buffer
            w->use_gz = true;
        }
        auto [ins_it, _] = writers.emplace(bc, std::move(w));
        return *(ins_it->second);
    };
    auto open_single_writer = [&](const std::string& path) -> std::unique_ptr<ref_writer> {
        auto w = std::make_unique<ref_writer>();
        bool is_gz = path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
        w->use_gz = is_gz;
        if (is_gz) {
            w->pigz_fp = pigz_out.open_pigz_pipe(path, nthreads);
            if (w->pigz_fp) {
                w->using_pigz = true;
                return w;
            }
            w->gz = gzopen(path.c_str(), "wb");
            if (!w->gz) throw std::runtime_error("Failed to open gzip output: " + path);
            gzbuffer(w->gz, 1 << 20);
        } else {
            w->out.open(path, std::ios::out | std::ios::binary);
            if (!w->out) throw std::runtime_error("Failed to open output: " + path);
        }
        return w;
    };

    auto close_writer = [&](ref_writer& w) {
        if (w.using_pigz && w.pigz_fp) {
            pigz_out.close_pigz_pipe(w.pigz_fp);
        } else if (w.gz) {
            gzclose(w.gz);
        } else if (w.out.is_open()) {
            w.out.close();
        }
    };

    auto write_record = [&](ref_writer& w, const read_streaming::sequence& r, bool use_lock) {
        std::string rec;
        if (r.is_fastq) {
            rec.reserve(r.id.size() + r.comment.size() + r.seq.size() + r.qual.size() + 16);
            rec += '@'; rec += r.id;
            if (!r.comment.empty()) { rec += '\t'; rec += r.comment; }
            rec += '\n'; rec += r.seq; rec += "\n+\n"; rec += r.qual; rec += '\n';
        } else {
            rec.reserve(r.id.size() + r.comment.size() + r.seq.size() + 8);
            rec += '>'; rec += r.id;
            if (!r.comment.empty()) { rec += '\t'; rec += r.comment; }
            rec += '\n'; rec += r.seq; rec += '\n';
        }
        auto do_write = [&]() {
            if (w.using_pigz) {
                pigz_out.pigz_write(w.pigz_fp, rec.data(), rec.size());
            } else if (w.gz) {
                gzwrite(w.gz, rec.data(), static_cast<unsigned int>(rec.size()));
            } else {
                w.out.write(rec.data(), rec.size());
            }
        };
        if (use_lock) {
            std::lock_guard<std::mutex> lock(w.mtx);
            do_write();
        } else {
            do_write();
        }
    };

    try {
        const size_t progress_interval = 500000;
        std::atomic<size_t> total{0};

        chunk_streaming<read_streaming::sequence, 
        std::function<void(std::vector<read_streaming::sequence>&, const std::string&)>> cs(5000);

        // set up single-writer if reformat only
        std::unique_ptr<ref_writer> single_writer;
        std::string single_tmp;
        if (do_reformat && !do_split) {
            boost::filesystem::path in(fastq_path);
            const std::string in_name = in.filename().string();
            const bool in_is_gz = in_name.size() >= 3 &&
                                  in_name.substr(in_name.size() - 3) == ".gz";
            const std::string tmp_name = in_name + ".tmp" + (in_is_gz ? ".gz" : "");
            single_tmp = (in.parent_path() / tmp_name).string();
            single_writer = open_single_writer(single_tmp);
        }

        auto chunk_func = [&](std::vector<read_streaming::sequence>& chunk, const std::string&) {
            for (auto& r : chunk) {
                auto cb = extract_cb(r.id, r.comment);
                std::optional<std::string> ub = extract_ub(r.comment);
                std::string qname = r.id;

                // If ID was previously collapsed, restore qname and recover CB/UB from the ID.
                if (auto parsed = parse_collapsed_id(r.id, reformat_sep); parsed.has_value()) {
                    auto [parsed_qname, parsed_cb, parsed_ub] = *parsed;
                    qname = parsed_qname;
                    cb = parsed_cb;
                    ub = parsed_ub;
                    r.id = qname;
                }

                if (do_coordinate) {
                    if (cb.has_value() && ub.has_value() && !cb->empty() && !ub->empty()) {
                        auto xy = resolve_xy_barcodes_from_cb(*cb, *spatial_ref);
                        if (xy.has_value()) {
                            const auto& x_bc = xy->first;
                            const auto& y_bc = xy->second;
                            auto xit = spatial_ref->x_index_by_bc.find(x_bc);
                            auto yit = spatial_ref->y_index_by_bc.find(y_bc);
                            if (xit != spatial_ref->x_index_by_bc.end() && yit != spatial_ref->y_index_by_bc.end()) {
                                int x_bin = project_axis_for_bin_size(xit->second, bin_size_um);
                                int y_bin = project_axis_for_bin_size(yit->second, bin_size_um);
                                r.id = make_coordinate_read_id(qname, x_bin, y_bin);
                                std::string spatial_bc = make_spatial_mapping_id(bin_size_um, x_bin, y_bin);
                                r.comment = make_visium_hd_comment(spatial_bc, *cb, *ub, make_bin_rg_tag(bin_size_um));
                            }
                        } else if (verbose) {
                            std::cout << "[reformat] skipped coordinate mapping for read with unresolvable CB pair: "
                                      << *cb << "\n";
                        }
                    } else if (verbose) {
                        std::cout << "[reformat] skipped coordinate mapping for read without parsable CB/UB: "
                                  << r.id << "\n";
                    }
                } else if (do_reformat) {
                    r.id = collapse_id(r.id, cb, ub, reformat_sep);
                    r.comment.clear(); // drop tags after collapsing
                }

                if (do_split) {
                    if (!cb || cb->empty()) continue;
                    ref_writer& w = get_writer(*cb);
                    write_record(w, r, true);
                } else {
                    write_record(*single_writer, r, false);
                }

                size_t now = total.fetch_add(1, std::memory_order_relaxed) + 1;
                if (verbose && now % progress_interval == 0) {
                    std::cout << "[reformat] processed " << now << " reads\n";
                }
            }
        };

        cs.process_chunks(fastq_path, chunk_func, nthreads, -1);

        if (do_split) {
            for (auto& kv : writers) {
                auto& w = *(kv.second);
                close_writer(w);
            }
        } else if (single_writer) {
            close_writer(*single_writer);
            // replace original
            boost::filesystem::path orig(fastq_path);
            boost::filesystem::path tmp(single_tmp);
            boost::filesystem::remove(orig);
            boost::filesystem::rename(tmp, orig);
        }

        if (verbose) {
            std::cout << "[reformat] wrote " << total.load() << " reads\n";
            if (do_split) {
                std::cout << "[reformat] barcodes: " << writers.size() << "\n";
            }
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[reformat][ERROR] " << ex.what() << "\n";
        for (auto& kv : writers) {
            auto& w = *(kv.second);
            close_writer(w);
        }
        return 1;
    }
}

// ============================================================================
// COMMAND: demux
// ============================================================================

int cmd_demux(int argc, char* argv[]) {
    std::string layout_key, fastq_path, custom_kit, global_whitelist_path, custom_whitelist_path, output_prefix, output_dir, log_file;
    bool verbose = false, max_verbose = false, split_bc = false, write_debug = false;
    std::optional<int> wl_mut, gen_mut;
    int nthreads = 1;
    size_t max_reads = -1;
    size_t chunk_size = 5000;
    std::string bc_corr_mode = "offensive";
    std::string joint_bc_mode = "default";

    const char* optstring = "l:q:k:g:c:R:M:m:n:z:o:d:F:wbt:vDh";
    struct option longopts[] = {
        {"layout",                      required_argument, nullptr, 'l'},
        {"fastq",                       required_argument, nullptr, 'q'},
        {"kit",                         required_argument, nullptr, 'k'},
        {"global_whitelist",            required_argument, nullptr, 'g'},
        {"custom_whitelist",            required_argument, nullptr, 'c'},
        {"bc_correction_mode",          required_argument, nullptr, 'R'},
        {"whitelist_mutation",          required_argument, nullptr, 'M'},
        {"generated_mutation",          required_argument, nullptr, 'm'},
        {"max_reads",                   required_argument, nullptr, 'n'},
        {"chunk_size",                  required_argument, nullptr, 'z'},
        {"output",                      required_argument, nullptr, 'o'},
        {"dir",                         required_argument, nullptr, 'd'},
        {"log-file",                    required_argument, nullptr, 'F'},
        {"write_dbg",                   no_argument,       nullptr, 'w'},
        {"bc_split",                    no_argument,       nullptr, 'b'},
        {"joint-bc-mode",               required_argument, nullptr,  1 },
        {"threads",                     required_argument, nullptr, 't'},
        {"verbose",                     no_argument,       nullptr, 'v'},
        {"max_verbose",                 no_argument,       nullptr, 'D'},
        {"help",                        no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
          case 'l': layout_key             = optarg;            break;
          case 'q': fastq_path             = optarg;            break;
          case 'k': custom_kit             = optarg;            break;
          case 'g': global_whitelist_path  = optarg;            break;
          case 'c': custom_whitelist_path  = optarg;            break;
          case 'R': bc_corr_mode           = optarg;            break;
          case 'M': wl_mut                 = std::stoi(optarg); break;
          case 'm': gen_mut                = std::stoi(optarg); break;
          case 'n': max_reads              = std::stoull(optarg); break;
          case 'z': chunk_size             = std::stoull(optarg); break;
          case 'o': output_prefix          = optarg;            break;
          case 'd': output_dir             = optarg;            break;
          case 'F': log_file               = optarg;            break;
          case 'w': write_debug            = true;              break;
          case 'b': split_bc               = true;              break;
          case 1: {
              joint_bc_mode = optarg ? optarg : "default";
              std::transform(joint_bc_mode.begin(), joint_bc_mode.end(), joint_bc_mode.begin(),
                             [](unsigned char ch){ return static_cast<char>(std::tolower(ch)); });
              if (joint_bc_mode == "auto") joint_bc_mode = "default";
              if (joint_bc_mode != "default" && joint_bc_mode != "strict") {
                  std::cerr << "[ERROR] --joint-bc-mode must be one of: default, strict\n";
                  return 1;
              }
              break;
          }
          case 't': nthreads               = std::stoi(optarg); break;
          case 'v': verbose                = true;              break;
          case 'D': max_verbose            = true;              break;
          case 'h': usage_demux(argv[0]);                       return 0;
          default:  usage_demux(argv[0]);                       return 1;
        }
    }

    if (layout_key.empty() || fastq_path.empty()) {
        std::cerr << "[ERROR] --layout and --fastq are required\n\n";
        usage_demux(argv[0]);
        return 1;
    }

    if (max_verbose) verbose = true;

    if (chunk_size == 0) {
        std::cerr << "[ERROR] chunk_size must be greater than 0\n";
        return 1;
    }

    // Setup logging
    std::ofstream log_stream;
    std::unique_ptr<tee_buff> tee_buf;
    std::streambuf* cout_backup = nullptr;
    std::streambuf* cerr_backup = nullptr;

    if (!log_file.empty()) {
        boost::filesystem::path log_path(log_file);
        if (log_path.has_parent_path()) {
            boost::filesystem::create_directories(log_path.parent_path());
        }
        
        log_stream.open(log_file);
        if (!log_stream.is_open()) {
            std::cerr << "[ERROR] Cannot open log file: " << log_file << "\n";
            return 1;
        }
        
        cout_backup = std::cout.rdbuf();
        cerr_backup = std::cerr.rdbuf();
        if (max_verbose) {
            // In max verbose mode, send all stdout/stderr to the log file to avoid flooding the terminal.
            std::cout.rdbuf(log_stream.rdbuf());
            std::cerr.rdbuf(log_stream.rdbuf());
            std::cerr << "[main] Max verbose output redirected to log file: " << log_file << "\n";
        } else {
            tee_buf = std::make_unique<tee_buff>(std::cout.rdbuf(), log_stream.rdbuf());
            std::cout.rdbuf(tee_buf.get());
            std::cerr << "[main] Logging output to: " << log_file << "\n";
        }
    }

    // Output paths
    boost::filesystem::path outdir = output_dir.empty()
        ? boost::filesystem::current_path()
        : boost::filesystem::path(output_dir);

    auto ensure_dir = [](const boost::filesystem::path& dir_path) -> bool {
        if (dir_path.empty()) return true;
        if (boost::filesystem::exists(dir_path)) return true;
        return boost::filesystem::create_directories(dir_path);
    };

    boost::filesystem::path outbase;
    if (output_prefix.empty()) {
        if (!ensure_dir(outdir)) {
            std::cerr << "[ERROR] Cannot create output directory: " << outdir.string() << "\n";
            return 1;
        }
        outbase = outdir / "output";
    } else {
        boost::filesystem::path prefix_path(output_prefix);
        if (prefix_path.is_absolute()) {
            outbase = prefix_path;
        } else {
            if (!ensure_dir(outdir)) {
                std::cerr << "[ERROR] Cannot create output directory: " << outdir.string() << "\n";
                return 1;
            }
            outbase = outdir / prefix_path;
        }
    }

    boost::filesystem::path outbase_parent = outbase.parent_path();
    if (!outbase_parent.empty() && !ensure_dir(outbase_parent)) {
        std::cerr << "[ERROR] Cannot create output directory: " << outbase_parent.string() << "\n";
        return 1;
    }

    std::string base = outbase.filename().string();
    std::string outdir_display = outbase_parent.empty() ? outdir.string() : outbase_parent.string();

    if (verbose) {
        std::cout << "============= Configuration =============\n"
                  << "  Layout key/path              : " << layout_key               << "\n"
                  << "  FASTA/Q input                : " << fastq_path               << "\n"
                  << "  custom kit                   : " << (custom_kit.empty() ? "[none]" : custom_kit) << "\n"
                  << "  global whitelist             : " << (global_whitelist_path.empty()? "[none]" : global_whitelist_path) << "\n"
                  << "  custom whitelist             : " << (custom_whitelist_path.empty()? "[none]" : custom_whitelist_path) << "\n"
                  << "  barcode correction mode      : " << bc_corr_mode             << "\n"
                  << "  joint barcode mode           : " << joint_bc_mode            << "\n"
                  << "  whitelist-generated mutations: " << (wl_mut ? std::to_string(*wl_mut) : "default(2)") << "\n"
                  << "  read-generated mutations     : " << (gen_mut ? std::to_string(*gen_mut) : "default(2)") << "\n"
                  << "  max reads                    : " << (max_reads == -1 ? "all" : std::to_string(max_reads)) << "\n"
                  << "  chunk size                   : " << chunk_size               << "\n"
                  << "  output directory             : " << outdir_display           << "\n"
                  << "  write debug files            : " << std::boolalpha << write_debug << "\n"
                  << "  filename base                : " << base                     << "\n"
                  << "  log file                     : " << (log_file.empty() ? "[none]" : log_file) << "\n"
                  << "  barcode split                : " << std::boolalpha << split_bc << "\n"
                  << "  threads                      : " << nthreads                 << "\n"
                  << "  verbose                      : " << std::boolalpha << verbose << "\n"
                  << "  max_threads avail.           : " << omp_get_max_threads()    << "\n\n";
    }

    std::cout << "[main] Starting main processing...\n";
    auto main_start = std::chrono::steady_clock::now();

    try {
        // ReadLayout setup
        ReadLayout read_layout;
        std::string layout_csv;
        if (!config_utils::check_if_custom_rl(layout_key)) {
            layout_csv = config_utils::get_read_layout(layout_key);
        } else {
            layout_csv = layout_key;
        }

        // Import or generate layout
        bool have_layout = boost::filesystem::exists(outbase.string() + "_layout.csv");
        bool have_map = boost::filesystem::exists(outbase.string() + "_position_map.csv");

        if (have_layout) {
            if (verbose) std::cout << "[import_read_layout] from " << outbase.string() + "_layout.csv\n";
            read_layout.import_read_layout(outbase.string() + "_layout.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[prep_new_layout] Generating read layout...\n";
            read_layout.prep_new_layout(layout_csv, verbose);
        }
        auto rl_gen_elapsed = std::chrono::steady_clock::now() - main_start;
        std::cout << "[main] Read layout generation time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(rl_gen_elapsed).count()
                  << " ms\n";

        memory_utils::get_rss();

        // Import or generate position map
        if (have_map) {
            if (verbose) std::cout << "[pos_map] Importing position map...\n";
            read_layout.import_position_map(outbase.string() + "_position_map.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[pos_map] Generating position map...\n";
            read_layout.generate_position_mapping();
        }
        auto pos_map_elapsed = std::chrono::steady_clock::now() - main_start;
        std::cout << "[main] Position map time: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(pos_map_elapsed).count()
                  << " ms\n";

        // Misalignment if needed
        if (!have_layout || !have_map) {
            read_layout.write_to_csv(outbase.string(), "layout");
            if (verbose) std::cout << "[misalignment_stats] Computing misalignment...\n";
            Misalignment_Setup mis(read_layout);
            mis.generate_misalignment_data(fastq_path, read_layout, nthreads);
            auto mis_time = std::chrono::steady_clock::now() - main_start;
            std::cout << "[main] Misalignment time: "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(mis_time).count()
                      << " ms\n";
            if (verbose) std::cout << "[main] Writing updated layout & map...\n";
            if (max_verbose) read_layout.display_read_layout();
            read_layout.write_to_csv(outbase.string(), "both");
        }

        memory_utils::get_rss();

        // Whitelist loading
        std::string kit_or_wl = custom_kit.empty() ? global_whitelist_path : custom_kit;
        std::optional<std::string> whitelist_path;
        if (!kit_or_wl.empty() && !custom_whitelist_path.empty()) {
            whitelist_path = kit_or_wl + ":" + custom_whitelist_path;
        }

        if (whitelist_path) {
            if (verbose) std::cout << "[main] Loading custom kit & whitelist...\n";
            read_layout.load_wl(whitelist_path.value(), wl_mut, verbose, nthreads);
        } else if (!custom_kit.empty()) {
            if (verbose) std::cout << "[main] Loading kit whitelist...\n";
            read_layout.load_wl(std::nullopt, wl_mut, verbose, nthreads);
        } else {
            if (verbose) std::cout << "[main] Loading all whitelists...\n";
            read_layout.load_wl(std::nullopt, wl_mut, verbose, nthreads);
        }

        // Demultiplex
        auto sigalign_start = std::chrono::steady_clock::now();
        std::cout << "[sigalign] Running sigalign with chunk size " << chunk_size << "...\n";
        memory_utils::get_rss();
        
        size_t max_reads_param = (max_reads == -1) ? -1 : max_reads;
        SigString::sigalign(fastq_path, read_layout, outbase.string(), gen_mut, max_verbose, nthreads, chunk_size, max_reads_param, write_debug, bc_corr_mode, joint_bc_mode);
        
        auto sig_time = std::chrono::steady_clock::now() - sigalign_start;
        std::cout << "[sigalign] End-to-end wall time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(sig_time).count()
                  << " s\n";

        // Save whitelist
        if (verbose) std::cout << "[main] Saving whitelist summary...\n";
        read_layout.save_wl(outbase.string() + "_whitelist.csv", false);

        auto final_elapsed = std::chrono::steady_clock::now() - main_start;
        std::cout << "[main] Total code duration: "
                  << std::chrono::duration_cast<std::chrono::seconds>(final_elapsed).count()
                  << " s\n";

        memory_utils::get_rss();

        if (split_bc) {
            if (verbose) std::cout << "[main] Splitting output by barcode...\n";
            auto split_start = std::chrono::steady_clock::now();
            try {
                //io_bc_split_utils::split_fastqas(outbase.string(), verbose, nthreads);
                auto split_time = std::chrono::steady_clock::now() - split_start;
                std::cout << "[main] Barcode splitting time: "
                          << std::chrono::duration_cast<std::chrono::seconds>(split_time).count()
                          << " s\n";
            } catch (const std::exception& ex) {
                std::cerr << "[WARNING] Barcode splitting failed: " << ex.what() << std::endl;
            }
        }
        
        // Restore cout if logging
        if (cout_backup) {
            std::cout.rdbuf(cout_backup);
        }
        if (cerr_backup) {
            std::cerr.rdbuf(cerr_backup);
        }
        if (log_stream.is_open()) {
            log_stream.close();
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what()
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n";
        
        // Restore cout on error too
        if (cout_backup) {
            std::cout.rdbuf(cout_backup);
        }
        if (cerr_backup) {
            std::cerr.rdbuf(cerr_backup);
        }
        if (log_stream.is_open()) {
            log_stream.close();
        }
        
        return 1;
    }
    return 0;
}

// ============================================================================
// COMMAND: list
// ============================================================================

int cmd_list(int argc, char* argv[]) {
    bool show_layouts = false;
    bool show_whitelists = false;

    const char* optstring = "h";
    struct option longopts[] = {
        {"layouts", no_argument, nullptr, 1},
        {"whitelists", no_argument, nullptr, 2},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 1: show_layouts = true; break;
            case 2: show_whitelists = true; break;
            case 'h': usage_list(argv[0]); return 0;
            default: usage_list(argv[0]); return 1;
        }
    }

    if (!show_layouts && !show_whitelists) {
        show_layouts = true;
        show_whitelists = true;
    }

    if (show_layouts) {
        std::cout << "[layouts]\n";
        std::vector<std::pair<std::string, std::string>> rows;
        const auto& layouts = config_utils::list_read_layouts();
        rows.reserve(layouts.size());
        for (const auto& kv : layouts) rows.emplace_back(kv.first, kv.second);
        std::sort(rows.begin(), rows.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        for (const auto& kv : rows) {
            try {
                std::string resolved = config_utils::get_read_layout(kv.first);
                std::cout << kv.first << "\t" << resolved << "\n";
            } catch (const std::exception& e) {
                std::cout << kv.first << "\t[unresolved] " << kv.second
                          << " (" << e.what() << ")\n";
            }
        }
        if (show_whitelists) std::cout << "\n";
    }

    if (show_whitelists) {
        std::cout << "[whitelists]\n";
        std::vector<std::pair<std::string, std::string>> rows;
        const auto& whitelists = whitelist_utils::list_whitelist_paths();
        rows.reserve(whitelists.size());
        for (const auto& kv : whitelists) rows.emplace_back(kv.first, kv.second);
        std::sort(rows.begin(), rows.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        for (const auto& kv : rows) {
            try {
                std::string resolved = whitelist_utils::kit_to_path(kv.first);
                std::cout << kv.first << "\t" << resolved << "\n";
            } catch (const std::exception& e) {
                std::cout << kv.first << "\t[unresolved] " << kv.second
                          << " (" << e.what() << ")\n";
            }
        }
    }

    return 0;
}

// ============================================================================
// COMMAND: modify
// ============================================================================

int cmd_modify(int argc, char* argv[]) {
    std::string layout_key;
    std::string whitelist_key;
    std::string add_path;
    bool do_remove = false;

    const char* optstring = "h";
    struct option longopts[] = {
        {"layout", required_argument, nullptr, 1},
        {"whitelist", required_argument, nullptr, 2},
        {"add", required_argument, nullptr, 3},
        {"remove", no_argument, nullptr, 4},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 1: layout_key = optarg; break;
            case 2: whitelist_key = optarg; break;
            case 3: add_path = optarg; break;
            case 4: do_remove = true; break;
            case 'h': usage_modify(argv[0]); return 0;
            default: usage_modify(argv[0]); return 1;
        }
    }

    const bool target_layout = !layout_key.empty();
    const bool target_whitelist = !whitelist_key.empty();
    if (target_layout == target_whitelist) {
        std::cerr << "[ERROR] Specify exactly one of --layout or --whitelist\n\n";
        usage_modify(argv[0]);
        return 1;
    }

    const bool do_add = !add_path.empty();
    if (do_add == do_remove) {
        std::cerr << "[ERROR] Specify exactly one of --add or --remove\n\n";
        usage_modify(argv[0]);
        return 1;
    }

    if (target_layout) {
        if (do_add) {
            config_utils::save_read_layout(layout_key, add_path);
            std::cout << "[modify] layout[" << layout_key << "] = " << add_path << "\n";
            try {
                std::string resolved = config_utils::get_read_layout(layout_key);
                std::cout << "[modify] resolved path: " << resolved << "\n";
            } catch (const std::exception& e) {
                std::cout << "[modify] warning: layout is currently unresolved (" << e.what() << ")\n";
            }
        } else {
            bool ok = config_utils::remove_read_layout(layout_key);
            if (!ok) {
                std::cerr << "[modify] layout key not found: " << layout_key << "\n";
                return 1;
            }
            std::cout << "[modify] removed layout key: " << layout_key << "\n";
        }
    } else {
        if (do_add) {
            whitelist_utils::set_whitelist_path(whitelist_key, add_path);
            std::cout << "[modify] whitelist[" << whitelist_key << "] = " << add_path << "\n";
            try {
                std::string resolved = whitelist_utils::kit_to_path(whitelist_key);
                std::cout << "[modify] resolved path: " << resolved << "\n";
            } catch (const std::exception& e) {
                std::cout << "[modify] warning: whitelist is currently unresolved (" << e.what() << ")\n";
            }
        } else {
            bool ok = whitelist_utils::remove_whitelist_path(whitelist_key);
            if (!ok) {
                std::cerr << "[modify] whitelist key not found: " << whitelist_key << "\n";
                return 1;
            }
            std::cout << "[modify] removed whitelist key: " << whitelist_key << "\n";
        }
    }

    return 0;
}

// ============================================================================
//  COMMAND: help (i need somebody (help (not just anybody)))
// ============================================================================

int cmd_help(int argc, char* argv[]) {
    if (argc < 2) {
        print_main_usage("rad");
        return 0;
    }
    
    std::string topic = argv[1];
    if (topic == "demux") {
        usage_demux("demux");
    } else if (topic == "prep") {
        usage_prep("prep");
    } else if (topic == "reformat") {
        usage_reformat("reformat");
    } else if (topic == "list") {
        usage_list("list");
    } else if (topic == "modify") {
        usage_modify("modify");
    } else {
        std::cerr << "Unknown command: " << topic << "\n\n";
        print_main_usage("rad");
        return 1;
    }
    return 0;
}

// ============================================================================
// MAIN: Command dispatcher
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_main_usage(argv[0]);
        return 1;
    }

    std::string command = argv[1];

    // Dispatch to subcommands
    if (command == "demux") {
        return cmd_demux(argc - 1, argv + 1);
    } else if (command == "prep") {
        return cmd_prep(argc - 1, argv + 1);
    } else if (command == "reformat") {
        return cmd_reformat(argc - 1, argv + 1);
    } else if (command == "list") {
        return cmd_list(argc - 1, argv + 1);
    } else if (command == "modify") {
        return cmd_modify(argc - 1, argv + 1);
    } else if (command == "help") {
        return cmd_help(argc - 1, argv + 1);
    } else if (command == "--help" || command == "-h") {
        print_main_usage(argv[0]);
        return 0;
    } else {
        std::cerr << "Unknown command: " << command << "\n\n";
        print_main_usage(argv[0]);
        return 1;
    }
}
