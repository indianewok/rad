// main.cpp
#include "include/rad/rad_headers.h"

// ============================================================================
// USAGE FUNCTIONS
// ============================================================================

static void print_main_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <command> [options]\n\n"
              << "Commands:\n"
              << "  demux       Run full demultiplexing pipeline\n"
              << "  prep_layout      Generate and display read layout structure\n"
              << "  position_map      Generate position map with misalignment statistics\n"
              << "  help        Show help for a specific command\n\n"
              << "Use '" << prog << " <command> --help' for command-specific options\n";
}

static void usage_demux(const char* prog) {
    std::cerr
      << "Usage: " << prog << " demux -l LAYOUT -q FASTQ [options]\n\n"
      << "Required:\n"
      << "  -l, --layout                      layout key (five_prime, three_prime, splitseq) or path\n"
      << "  -q, --fastq                       input FASTQ file\n\n"
      << "Optional:\n"
      << "  -k, --kit                         use this kit's default whitelist\n"
      << "  -g, --global_whitelist            path to custom global whitelist CSV\n"
      << "  -c, --custom_whitelist            path to custom whitelist CSV\n"
      << "  -R, --bc_correction_mode          'offensive' (default) or 'defensive'\n"
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

static void usage_layout(const char* prog) {
    std::cerr
      << "Usage: " << prog << " layout -l LAYOUT [options]\n\n"
      << "Generate and display read layout structure from CSV\n\n"
      << "Required:\n"
      << "  -l, --layout                      layout key or path to layout CSV\n\n"
      << "Optional:\n"
      << "  -o, --output                      output base path (saves to _layout.csv)\n"
      << "  -v, --verbose                     verbose mode\n"
      << "  -h, --help                        show this help\n";
}

static void usage_position_map(const char* prog) {
    std::cerr
      << "Usage: " << prog << " position_map -l LAYOUT -q FASTQ -o OUTPUT [options]\n\n"
      << "Generate position map with misalignment statistics\n"
      << "This command runs misalignment analysis and generates both\n"
      << "layout and position_map CSV files.\n\n"
      << "Required:\n"
      << "  -l, --layout                      layout key or path to layout CSV\n"
      << "  -q, --fastq                       input FASTQ file for misalignment stats\n"
      << "  -o, --output                      output base path\n\n"
      << "Optional:\n"
      << "  -n, --max_reads                   max reads for sampling (default: 50000)\n"
      << "  -t, --threads                     number of threads (default: 1)\n"
      << "  -v, --verbose                     verbose mode\n"
      << "  -D, --max_verbose                 maximum verbosity (debug)\n"
      << "  -h, --help                        show this help\n";
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
// COMMAND: layout
// ============================================================================

int cmd_layout(int argc, char* argv[]) {
    std::string layout_key;
    std::string output_base;
    bool verbose = false;
    
    const char* optstring = "l:o:vh";
    struct option longopts[] = {
        {"layout",  required_argument, nullptr, 'l'},
        {"output",  required_argument, nullptr, 'o'},
        {"verbose", no_argument,       nullptr, 'v'},
        {"help",    no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'l': layout_key = optarg; break;
            case 'o': output_base = optarg; break;
            case 'v': verbose = true; break;
            case 'h': usage_layout(argv[0]); return 0;
            default: usage_layout(argv[0]); return 1;
        }
    }
    
    if (layout_key.empty()) {
        std::cerr << "[ERROR] --layout is required\n\n";
        usage_layout(argv[0]);
        return 1;
    }
    
    try {
        std::cout << "[layout] Generating layout structure...\n\n";
        
        ReadLayout read_layout;
        std::string layout_csv;
        
        // Determine if it's a built-in layout or custom path
        if (!config_utils::check_if_custom_rl(layout_key)) {
            layout_csv = config_utils::get_read_layout(layout_key);
            std::cout << "[layout] Using built-in layout: " << layout_key << "\n";
        } else {
            layout_csv = layout_key;
            std::cout << "[layout] Using custom layout: " << layout_key << "\n";
        }
        
        // Prep layout with verbosity
        read_layout.prep_new_layout(layout_csv, verbose);
        
        std::cout << "\n========================================\n";
        std::cout << "         LAYOUT STRUCTURE               \n";
        std::cout << "========================================\n\n";
        
        // Display the layout
        read_layout.display_read_layout();
        
        std::cout << "\n========================================\n";
        std::cout << "Layout contains " << read_layout.size() << " elements\n";
        std::cout << "========================================\n";
        
        // Optionally save to file
        if (!output_base.empty()) {
            read_layout.write_to_csv(output_base, "layout");
            std::cout << "\n[layout] Saved to: " << output_base << "_layout.csv\n";
        }
        
        return 0;
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << "\n";
        return 1;
    }
}

// ============================================================================
// COMMAND: position_map
// ============================================================================

int cmd_position_map(int argc, char* argv[]) {
    std::string layout_key, fastq_path, output_base;
    bool verbose = false, max_verbose = false;
    int nthreads = 1;
    size_t max_reads = 50000;
    
    const char* optstring = "l:q:o:n:t:vDh";
    struct option longopts[] = {
        {"layout",      required_argument, nullptr, 'l'},
        {"fastq",       required_argument, nullptr, 'q'},
        {"output",      required_argument, nullptr, 'o'},
        {"max_reads",   required_argument, nullptr, 'n'},
        {"threads",     required_argument, nullptr, 't'},
        {"verbose",     no_argument,       nullptr, 'v'},
        {"max_verbose", no_argument,       nullptr, 'D'},
        {"help",        no_argument,       nullptr, 'h'},
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
            case 'h': usage_position_map(argv[0]); return 0;
            default: usage_position_map(argv[0]); return 1;
        }
    }
    
    if (layout_key.empty() || fastq_path.empty() || output_base.empty()) {
        std::cerr << "[ERROR] --layout, --fastq, and --output are required\n\n";
        usage_position_map(argv[0]);
        return 1;
    }
    
    try {
        std::cout << "[position_map] Generating position map with misalignment statistics...\n";
        auto start = std::chrono::steady_clock::now();
        
        ReadLayout read_layout;
        std::string layout_csv;
        
        // Load layout
        if (!config_utils::check_if_custom_rl(layout_key)) {
            layout_csv = config_utils::get_read_layout(layout_key);
            if (verbose) std::cout << "[position_map] Using built-in layout: " << layout_key << "\n";
        } else {
            layout_csv = layout_key;
            if (verbose) std::cout << "[position_map] Using custom layout: " << layout_key << "\n";
        }
        
        // Prep layout
        if (verbose) std::cout << "[position_map] Preparing layout...\n";
        read_layout.prep_new_layout(layout_csv, max_verbose);
        
        // Save initial layout
        read_layout.write_to_csv(output_base, "layout");
        
        // Run misalignment analysis
        if (verbose) std::cout << "[position_map] Computing misalignment statistics from FASTQ...\n";
        Misalignment_Setup mis(read_layout);
        mis.generate_misalignment_data(fastq_path, read_layout, nthreads, max_reads);
        
        // Generate position mapping
        if (verbose) std::cout << "[position_map] Generating position mapping...\n";
        read_layout.generate_position_mapping();
        
        // Display if max verbose
        if (max_verbose) {
            std::cout << "\n";
            read_layout.display_read_layout();
        }
        
        // Write updated files
        if (verbose) std::cout << "[position_map] Writing output files...\n";
        read_layout.write_to_csv(output_base, "both");
        
        auto elapsed = std::chrono::steady_clock::now() - start;
        std::cout << "\n[position_map] Complete! Time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(elapsed).count()
                  << " seconds\n";
        
        std::cout << "\n[position_map] Output files:\n"
                  << "  " << output_base << "_layout.csv\n"
                  << "  " << output_base << "_position_map.csv\n";
        
        return 0;
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << "\n";
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

    if (!boost::filesystem::exists(outdir)) {
        if (!boost::filesystem::create_directories(outdir)) {
            std::cerr << "[ERROR] Cannot create output directory: " << outdir.string() << "\n";
            return 1;
        }
    }

    std::string base = output_prefix.empty() ? "output" : output_prefix;
    boost::filesystem::path outbase = outdir / base;

    if (verbose) {
        std::cout << "============= Configuration =============\n"
                  << "  Layout key/path              : " << layout_key               << "\n"
                  << "  FASTA/Q input                : " << fastq_path               << "\n"
                  << "  custom kit                   : " << (custom_kit.empty() ? "[none]" : custom_kit) << "\n"
                  << "  global whitelist             : " << (global_whitelist_path.empty()? "[none]" : global_whitelist_path) << "\n"
                  << "  custom whitelist             : " << (custom_whitelist_path.empty()? "[none]" : custom_whitelist_path) << "\n"
                  << "  barcode correction mode      : " << bc_corr_mode             << "\n"
                  << "  whitelist-generated mutations: " << (wl_mut ? std::to_string(*wl_mut) : "default(2)") << "\n"
                  << "  read-generated mutations     : " << (gen_mut ? std::to_string(*gen_mut) : "default(2)") << "\n"
                  << "  max reads                    : " << (max_reads == -1 ? "all" : std::to_string(max_reads)) << "\n"
                  << "  chunk size                   : " << chunk_size               << "\n"
                  << "  output directory             : " << outdir.string()          << "\n"
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
            if (verbose) std::cout << "[prep_read_layout] Generating read layout...\n";
            read_layout.prep_new_layout(layout_csv, max_verbose);
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
        SigString::sigalign(fastq_path, read_layout, outbase.string(), gen_mut, max_verbose, nthreads, chunk_size, max_reads_param, write_debug, bc_corr_mode);
        
        auto sig_time = std::chrono::steady_clock::now() - sigalign_start;
        std::cout << "[sigalign] Time: "
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
//  help (i need somebody (help (not just anybody)))
// ============================================================================

int cmd_help(int argc, char* argv[]) {
    if (argc < 2) {
        print_main_usage("rad");
        return 0;
    }
    
    std::string topic = argv[1];
    if (topic == "demux") {
        usage_demux("rad");
    } else if (topic == "layout") {
        usage_layout("rad");
    } else if (topic == "position_map") {
        usage_position_map("rad");
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
    } else if (command == "prep_layout") {
        return cmd_layout(argc - 1, argv + 1);
    } else if (command == "position_map") {
        return cmd_position_map(argc - 1, argv + 1);
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
