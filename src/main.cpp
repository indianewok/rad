// main.cpp
#include "include/rad/rad_headers.h"

static void usage(const char* prog) {
    std::cerr
      << "Usage: " << prog
      << " -l LAYOUT -q FASTQ [-k KIT] [-g GLOBAL_WL] [-c CUSTOM_WL] [-o PREFIX] [-d DIR] [-b] [-t THREADS] [-v] [-D] [-h]\n"
      << "  -l, --layout                      layout key (five_prime, three_prime, splitseq) or a custom path\n"
      << "  -q, --fastq                       input FASTQ file\n"
      << "  -k, --kit                         use this kit's default whitelist\n"
      << "  -g, --global_whitelist            path to a single custom global whitelist CSV(i.e, all barcodes that could potentially appear in this dataset)\n"
      << "  -c, --custom_whitelist            path to a single custom whitelist CSV (i.e, barcodes that you've already found in short-read matched data)\n"
      << "  -R, --bc_correction_mode          either 'offensive' (default) or 'defensive', changes the behavior of barcode correction\n"
      << "  -M, --whitelist_mutation          mutation space of each barcode in the user-passed whitelist (default: 2)\n"
      << "  -m, --generated_mutation          mutation space of each barcode in the user-passed whitelist (default: 2)\n"
      << "  -n, --max_reads                   maximum number of reads to process (default: all)\n"
      << "  -z, --chunk_size                  chunk size for processing reads (default: 5000)\n"
      << "  -o, --output                      filename prefix  (default: output)\n"
      << "  -d, --dir                         output directory (default: current directory)\n"
      << "  -w, --write_dbg                   writes all debug .sig files, metrics.tsv, and .csv (debug only)\n"
      << "  -b, --bc_split                    write reads into per-barcode FASTQ files\n"
      << "  -t, --threads                     number of threads (default: 1)\n"
      << "  -v, --verbose                     verbose mode\n"
      << "  -D, --max_verbose                 max verbose level (debug only)\n"
      << "  -h, --help                        prints this menu\n";
}

int main(int argc, char* argv[]) {
    std::string layout_key, fastq_path, custom_kit, global_whitelist_path, custom_whitelist_path, output_prefix, output_dir;
    bool verbose = false, max_verbose = false, split_bc = false, write_debug = false;
    std::optional<int> wl_mut, gen_mut;
    int nthreads = 1;
    size_t max_reads = -1;  // -1 means unlimited
    size_t chunk_size = 5000;  // Default chunk size
    std::string bc_corr_mode = "offensive"; // default mode

    const char* optstring = "l:q:k:g:c:R:M:m:S:s:n:z:o:d:w:bt:vDh";
    struct option longopts[] = {
        {"layout",                      required_argument, nullptr, 'l'},
        {"fastq",                       required_argument, nullptr, 'q'},
        {"kit",                         required_argument, nullptr, 'k'},
        {"global_whitelist",            required_argument, nullptr, 'g'},
        {"custom_whitelist",            required_argument, nullptr, 'c'},
        {"bc_corr_mode",                required_argument, nullptr, 'R'},
        {"whitelist_mutation",          required_argument, nullptr, 'M'},
        {"generated_mutation",          required_argument, nullptr, 'm'},
        {"max_reads",                   required_argument, nullptr, 'n'},
        {"chunk_size",                  required_argument, nullptr, 'z'},
        {"output",                      required_argument, nullptr, 'o'},
        {"dir",                         required_argument, nullptr, 'd'},
        {"write_dbg",                   no_argument,       nullptr, 'w'},
        {"bc_split",                    no_argument,       nullptr, 'b'},
        {"threads",                     required_argument, nullptr, 't'},
        {"verbose",                     no_argument,       nullptr, 'v'},
        {"max_verbose",                 no_argument,       nullptr, 'D'},
        {"help",                        no_argument,       nullptr, 'h'},
        {nullptr, 0,nullptr, 0}
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
          case 'n': max_reads            = std::stoull(optarg); break;
          case 'z': chunk_size           = std::stoull(optarg); break;
          case 'o': output_prefix          = optarg;            break;
          case 'd': output_dir             = optarg;            break;
          case 'w': write_debug            = true;              break;
          case 'b': split_bc               = true;              break;
          case 't': nthreads               = std::stoi(optarg); break;
          case 'v': verbose                = true;              break;
          case 'D': max_verbose            = true;              break;
          case 'h': usage(argv[0]);                             return 0;
          default:  usage(argv[0]);                             return 1;
        }
    }

    if (layout_key.empty() || fastq_path.empty()) {
        std::cerr << "[ERROR] --layout and --fastq are required\n\n";
        usage(argv[0]);
        return 1;
    }

    if (max_verbose){
        verbose = true;
    }

    // Validate chunk_size
    if (chunk_size == 0) {
        std::cerr << "[ERROR] chunk_size must be greater than 0\n";
        return 1;
    }
    
    //
    // ─── OUTPUT PATHS ────────────────────────────────────────────────────────────
    //

    // Default output directory = current working directory
    boost::filesystem::path outdir = output_dir.empty()
        ? boost::filesystem::current_path()
        : boost::filesystem::path(output_dir);

    // Create it if necessary
    if (!boost::filesystem::exists(outdir)) {
        if (!boost::filesystem::create_directories(outdir)) {
            std::cerr << "[ERROR] Cannot create output directory: "
                      << outdir.string() << "\n";
            return 1;
        }
    }

    // Filename base
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

        // Demultiplex with chunk_size parameter
        auto sigalign_start = std::chrono::steady_clock::now();
        std::cout << "[sigalign] Running sigalign with chunk size " << chunk_size << "...\n";
        memory_utils::get_rss();
        
        // if max reads is -1, it's unlimited--if not, then default to the max number of reads that are user-specified
        size_t max_reads_param = (max_reads == -1) ? -1 : max_reads;
        // main function
        SigString::sigalign(fastq_path, read_layout, outbase.string(), gen_mut, max_verbose, nthreads, chunk_size, max_reads_param, write_debug, bc_corr_mode);
        
        auto sig_time = std::chrono::steady_clock::now() - sigalign_start;
        std::cout << "[sigalign] Time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(sig_time).count()
                  << " s\n";

        // Save whitelist summary
        if (verbose) std::cout << "[main] Saving whitelist summary...\n";
        read_layout.save_stats_wl(outbase.string() + "_whitelist_stats.csv", false);
        read_layout.save_wl(outbase.string() + "_whitelist.csv", false);

        auto final_elapsed = std::chrono::steady_clock::now() - main_start;
        std::cout << "[main] Completed in "
                  << std::chrono::duration_cast<std::chrono::seconds>(final_elapsed).count()
                  << " s\n";

        memory_utils::get_rss();

        if (split_bc) {
            if (verbose) std::cout << "[main] Splitting output by barcode...\n";
            auto split_start = std::chrono::steady_clock::now();
            try {
                io_bc_split_utils::split_fastqas(outbase.string(), verbose, nthreads);
                auto split_time = std::chrono::steady_clock::now() - split_start;
                std::cout << "[main] Barcode splitting time: "
                          << std::chrono::duration_cast<std::chrono::seconds>(split_time).count()
                          << " s\n";
            } catch (const std::exception& ex) {
                std::cerr << "[WARNING] Barcode splitting failed: " << ex.what() << std::endl;
            }
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what()
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n";
        return 1;
    }
    return 0;
}