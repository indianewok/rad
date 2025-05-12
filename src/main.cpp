// main.cpp
#include "include/rad/rad_headers.h"
/*
static void usage(const char* prog) {
    std::cerr
      << "Usage: " << prog
      << " -l LAYOUT -q FASTQ [-k KIT] [-w PATH] [-o PREFIX] [-d DIR] [-t THREADS] [-v VERBOSE] [-V MAX_VERBOSE]\n\n"
      << "  -l, --layout            layout key (five_prime|three_prime|splitseq) or a custom path\n"
      << "  -q, --fastq             input FASTQ file\n"
      << "  -k, --kit               use this kit's default whitelist\n"
      << "  -g, --global_whitelist  path to a single global whitelist CSV (all possible barcodes, overrides read layout)\n"
      << "  -c, --custom_whitelist  path to a single whitelist CSV (selected barcodes within the global whitelist, overrides read layout)\n"
      << "  -m, --mutation          mutation space of each barcode (look for barcodes up to m edit distance away)\n"
      << "  -s, --shift             shift space of each barcode (shifts barcodes s bases to the left or right)\n"
      << "  -o, --output            filename prefix  (default: [output])\n"
      << "  -d, --dir               output directory (default: [$HOME/Desktop])\n"
      << "  -b, --bc_split          write reads into per-barcode FASTQ files\n"
      << "  -t, --threads           number of threads (default: [1])\n"
      << "  -v, --verbose           verbose mode\n"
      << "  -D, --max_verbose       max verbose level: for debug only, generates a lot of text\n"
      << "  -h, --help              prints this menu\n";
}

int main(int argc, char* argv[]) {
    std::string layout_key, fastq_path, kit, global_whitelist_path, custom_whitelist_path, output_prefix, output_dir;
    std::optional<int> shift, mut;
    int nthreads = 1;
    bool verbose, max_verbose = false;
    bool split_bc = false;

    const char* optstring = "l:q:g:c:k:m:s:o:d:t:vhD";
    struct option longopts[] = {
        {"layout",    required_argument, nullptr, 'l'},
        {"fastq",     required_argument, nullptr, 'q'},
        {"kit",       required_argument, nullptr, 'k'},
        {"global_whitelist", required_argument, nullptr, 'g'},
        {"custom_whitelist", required_argument, nullptr, 'c'},
        {"mutation",  required_argument, nullptr, 'm'},
        {"shift",     required_argument, nullptr, 's'},
        {"output",    required_argument, nullptr, 'o'},
        {"dir",       required_argument, nullptr, 'd'},
        {"bc_split",  no_argument,       nullptr, 'b'},
        {"threads",   required_argument, nullptr, 't'},
        {"verbose",   no_argument,       nullptr, 'v'},
        {"max_verbose",no_argument,      nullptr, 'D'},
        {"help",      no_argument,       nullptr, 'h'},
        {nullptr,     0,                 nullptr,  0 }
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
          case 'l': layout_key     = optarg;            break;
          case 'q': fastq_path     = optarg;            break;
          case 'k': kit            = optarg;            break;
          case 'g': global_whitelist_path = optarg;     break;
          case 'c': custom_whitelist_path = optarg;     break;
          case 'm': mut            = std::stoi(optarg); break;
          case 's': shift          = std::stoi(optarg); break;
          case 'o': output_prefix  = optarg;            break;
          case 'b': split_bc       = true;              break;
          case 'd': output_dir     = optarg;            break;
          case 't': nthreads       = std::stoi(optarg); break;
          case 'v': verbose        = true;              break;
          case 'D': max_verbose   = true;               break;
          case 'h': usage(argv[0]); return 0;
          default:  usage(argv[0]); return 1;
        }
    }

    if (layout_key.empty() || fastq_path.empty()) {
        std::cerr << "[ERROR] --layout and --fastq are required\n\n";
        usage(argv[0]);
        return 1;
    }

    // the hodgepodge mess that comes with trying to build all the output directories, etc., etc.
    // need to compartmentalize split_bc into its own function for sure, but not sure how to do it re: streaming output
    std::filesystem::path outdir;
    std::string base = output_prefix.empty() ? "output" : output_prefix;

    std::filesystem::path outbase = outdir / base;
    if (!output_dir.empty()) {
        outdir = output_dir;
    } else if (auto h = std::getenv("HOME")) {
        outdir = std::filesystem::path(h) / "Desktop";
    } else {
        outdir = std::filesystem::current_path();
    }
    std::error_code ec;
    std::filesystem::create_directories(outdir, ec);
    if (ec && ec != std::errc::file_exists) {
        std::cerr << "[ERROR] Cannot create " << outdir << ": " << ec.message() << "\n";
        return 1;
    }

    //this doesn't go anywhere yet--it just directs the creation of the folder bc_split
    if (split_bc) {
        std::unordered_map<std::string,std::ofstream> bc_streams;
        std::string bc_dir = outdir / (base + "bc_split");
        std::error_code ec;
        std::filesystem::create_directories(bc_dir, ec);
        if (ec && ec != std::errc::file_exists) {
            throw std::runtime_error("Cannot create barcode folder: " + bc_dir);
        }
    }
    
    // setting up all of the verbosity to run
    if(max_verbose) verbose = true;

    if (verbose) {
        std::cout << "=== Configuration ===\n"
                  << " Layout key/path    : " << layout_key    << "\n"
                  << " FASTQ input        : " << fastq_path    << "\n"
                  << " kit                : " << (kit.empty() ? "[none specified]" : kit) << "\n"
                  << " global whitelist   : " << (global_whitelist_path.empty() ? "[none specified]" : global_whitelist_path) << "\n"
                  << " custom whitelist   : " << (custom_whitelist_path.empty() ? "[none specified]" : custom_whitelist_path) << "\n"
                  << " output directory   : " << outdir        << "\n"
                  << " filename base      : " << base          << "\n"
                  << " threads            : " << nthreads   << "\n"
                  << " verbose            : " << std::boolalpha << verbose << "\n"
                  << " max_verbose        : " << std::boolalpha << max_verbose << "\n"
                  << " max threads avail. : " << omp_get_max_threads()
                  << "\n\n";
    }

    static const auto main_start = std::chrono::steady_clock::now();
    std::cout << "[main] Starting main processing...\n";

    try {
        // setting up read layout
       ReadLayout read_layout;
       std::string layout_csv;
       if(!config_utils::check_if_custom_rl(layout_key))
         layout_csv = config_utils::get_read_layout(layout_key);
        else {
            layout_csv = layout_key;
        }

        // import or generate layout + map
        bool have_layout = std::filesystem::exists(outbase.string() + "_layout.csv");
        bool have_map = std::filesystem::exists(outbase.string() + "_position_map.csv");

        // do we have a layout?
        if (have_layout) {
            if (verbose) std::cout << "[import_read_layout] Importing read layout from " << outbase.string() + "_layout.csv" << "\n";
            read_layout.import_read_layout(outbase.string() + "_layout.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[prep_read_layout] Generating read layout...\n";
            read_layout.prep_new_layout(layout_csv, max_verbose);
        }

        //timing read layout generation
        auto rl_gen_time = std::chrono::steady_clock::now();
        auto rl_gen_elapsed = rl_gen_time - main_start;
        std::cout << "[main] Read layout generation time: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(rl_gen_elapsed).count()  << " ms\n";

        // do we have a position map? 
        if (have_map) {
            if (verbose) std::cout << "[pos_map] Importing position map from " << outbase.string() + "_position_map.csv" << "\n";
            read_layout.import_position_map(outbase.string() + "_position_map.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[pos_map] Generating position map...\n";
            read_layout.generate_position_mapping();
        }

        // setting up timing for position mapping
        auto pos_map_time = std::chrono::steady_clock::now();
        auto pos_map_elapsed = pos_map_time - main_start;
        std::cout << "[main] Position map generation time: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(pos_map_elapsed).count() << " ms\n";

        //  compute misalignment if we regenerated either
        if (!have_layout || !have_map) {
            read_layout.write_to_csv(outbase.string(), "layout");
            if (verbose) std::cout << "[misalignment_stats] Computing misalignment...";
            Misalignment_Setup mis(read_layout);
            mis.generate_misalignment_data(fastq_path, read_layout, nthreads);

            // timing misalignment threshold
            auto misalignment_time = std::chrono::steady_clock::now();
            auto mis_time_elapsed = misalignment_time - main_start;
            std::cout << "[main] Misalignment threshold generation time: " 
                      << std::chrono::duration_cast<std::chrono::milliseconds>(mis_time_elapsed).count() << " ms\n";
    
            if(verbose) std::cout << "[main] Writing generated layout and filled-out map...\n";
            if(max_verbose) read_layout.display_read_layout();
            read_layout.write_to_csv(outbase.string(), "both");
        }

        // wl loading
        std::string kit_or_wl = kit.empty() ? global_whitelist_path : kit;
        std::optional<std::string> whitelist_path;

        // wl kit 
        if(!kit_or_wl.empty() && !custom_whitelist_path.empty()) {
           whitelist_path = kit_or_wl + ":" + custom_whitelist_path;
        }
        if (whitelist_path.has_value()) {
            if (verbose) std::cout << "[main] Loading custom kit & whitelist: " << whitelist_path.value() << "\n";
            read_layout.load_wl(whitelist_path.value(), shift, mut, verbose, nthreads);
        }
        else if (!kit.empty()) {
            if (verbose) std::cout << "[main] Loading kit whitelist: " << kit << "\n";
            read_layout.load_wl(whitelist_path, shift, mut, max_verbose, nthreads);
        }
        else {
            if (verbose) std::cout << "[main] Loading all whitelists from read layout...\n";
            read_layout.load_wl(whitelist_path, shift, mut, max_verbose, nthreads);
        }

        // master demuxn fxn
        auto sigalign_start = std::chrono::steady_clock::now();
        std::cout << "[sigalign] Running sigalign...\n";
        SigString::sigalign(fastq_path, read_layout, outbase.string(), verbose=max_verbose, nthreads, max_reads=-1);
        auto sigalign_stop = std::chrono::steady_clock::now();
        auto sig_time = sigalign_stop - sigalign_start;
        std::cout << "[sigalign] Total demultiplexing time: " 
                  << std::chrono::duration_cast<std::chrono::seconds>(sig_time).count() << " seconds\n";

        // save whitelist summary
        if (verbose) std::cout << "[main] Saving whitelist summary to " << outbase.string() + "_whitelist.csv" << "\n";
        read_layout.save_wl(outbase.string() + "_whitelist.csv", verbose, full=false);

        auto final_time = std::chrono::steady_clock::now();
        auto final_elapsed = final_time - main_start;
        std::cout << "[main] Processing completed!\n";
        std::cout << "[main] Total run time: " 
                  << std::chrono::duration_cast<std::chrono::seconds>(final_elapsed).count() << " seconds\n";
        
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << " for " << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
        return 1;
    }

    return 0;
}
*/

static void usage(const char* prog) {
    std::cerr
      << "Usage: " << prog
      << " -l LAYOUT -q FASTQ [-k KIT] [-g GLOBAL_WL] [-c CUSTOM_WL] [-o PREFIX] [-d DIR] [-b] [-t THREADS] [-v] [-D] [-h]\n"
      << "  -l, --layout            layout key (five_prime, three_prime, splitseq) or a custom path\n"
      << "  -q, --fastq             input FASTQ file\n"
      << "  -k, --kit               use this kit's default whitelist\n"
      << "  -g, --global_whitelist  path to a single custom global whitelist (i.e, all barcodes that could potentially appear in this dataset) CSV\n"
      << "  -c, --custom_whitelist  path to a single custom whitelist CSV (i.e, barcodes that you've already found in short-read matched data)\n"
      << "  -m, --mutation          mutation space of each barcode\n"
      << "  -s, --shift             shift space of each barcode\n"
      << "  -o, --output            filename prefix  (default: output)\n"
      << "  -d, --dir               output directory (default: current directory)\n"
      << "  -b, --bc_split          write reads into per-barcode FASTQ files\n"
      << "  -t, --threads           number of threads (default: 1)\n"
      << "  -v, --verbose           verbose mode\n"
      << "  -D, --max_verbose       max verbose level (debug only)\n"
      << "  -h, --help              prints this menu\n";
}

int main(int argc, char* argv[]) {
    std::string layout_key, fastq_path, custom_kit, global_whitelist_path, custom_whitelist_path, output_prefix, output_dir;
    bool verbose = false, max_verbose = false, split_bc = false;
    std::optional<int> shift, mut;
    int nthreads = 1;

    const char* optstring = "l:q:k:g:c:m:s:o:d:bt:vDh";
    struct option longopts[] = {
        {"layout",            required_argument, nullptr, 'l'},
        {"fastq",             required_argument, nullptr, 'q'},
        {"kit",               required_argument, nullptr, 'k'},
        {"global_whitelist",  required_argument, nullptr, 'g'},
        {"custom_whitelist",  required_argument, nullptr, 'c'},
        {"mutation",          required_argument, nullptr, 'm'},
        {"shift",             required_argument, nullptr, 's'},
        {"output",            required_argument, nullptr, 'o'},
        {"dir",               required_argument, nullptr, 'd'},
        {"bc_split",          no_argument,       nullptr, 'b'},
        {"threads",           required_argument, nullptr, 't'},
        {"verbose",           no_argument,       nullptr, 'v'},
        {"max_verbose",       no_argument,       nullptr, 'D'},
        {"help",              no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
          case 'l': layout_key             = optarg;            break;
          case 'q': fastq_path             = optarg;            break;
          case 'k': custom_kit                    = optarg;            break;
          case 'g': global_whitelist_path  = optarg;            break;
          case 'c': custom_whitelist_path  = optarg;            break;
          case 'm': mut   = std::stoi(optarg);                  break;
          case 's': shift = std::stoi(optarg);                  break;
          case 'o': output_prefix          = optarg;            break;
          case 'd': output_dir             = optarg;            break;
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
    if (max_verbose) verbose = true;

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

    // Per-barcode split folder
    boost::filesystem::path bc_dir;
    if (split_bc) {
        bc_dir = outdir / (base + "_bc_split");
        if (!boost::filesystem::exists(bc_dir) && !boost::filesystem::create_directories(bc_dir)) {
            std::cerr << "[ERROR] Cannot create barcode-split directory: "
                      << bc_dir.string() << "\n";
            return 1;
        }
    }

    if (verbose) {
        std::cout << "============= Configuration =============\n"
                  << "  Layout key/path    : " << layout_key               << "\n"
                  << "  FASTQ input        : " << fastq_path               << "\n"
                  << "  custom kit         : " << (custom_kit.empty() ? "[none]" : custom_kit) << "\n"
                  << "  global whitelist   : " << (global_whitelist_path.empty()? "[none]" : global_whitelist_path) << "\n"
                  << "  custom whitelist   : " << (custom_whitelist_path.empty()? "[none]" : custom_whitelist_path) << "\n"
                  << "  output directory   : " << outdir.string()          << "\n"
                  << "  filename base      : " << base                     << "\n"
                  << "  threads            : " << nthreads                 << "\n"
                  << "  verbose            : " << std::boolalpha <<verbose  << "\n"
                  << "  max_threads avail. : " << omp_get_max_threads()    << "\n\n";
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
        bool have_map    = boost::filesystem::exists(outbase.string() + "_position_map.csv");

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

        // Whitelist loading
        std::string kit_or_wl = custom_kit.empty() ? global_whitelist_path : custom_kit;
        std::optional<std::string> whitelist_path;
        if (!kit_or_wl.empty() && !custom_whitelist_path.empty()) {
            whitelist_path = kit_or_wl + ":" + custom_whitelist_path;
        }
        if (whitelist_path) {
            if (verbose) std::cout << "[main] Loading custom kit & whitelist...\n";
            read_layout.load_wl(whitelist_path.value(), shift, mut, verbose, nthreads);
        } else if (!custom_kit.empty()) {
            if (verbose) std::cout << "[main] Loading kit whitelist...\n";
            read_layout.load_wl(std::nullopt, shift, mut, verbose, nthreads);
        } else {
            if (verbose) std::cout << "[main] Loading all whitelists...\n";
            read_layout.load_wl(std::nullopt, shift, mut, verbose, nthreads);
        }

        // Demultiplex
        auto sigalign_start = std::chrono::steady_clock::now();
        if (verbose) std::cout << "[sigalign] Running sigalign...\n";
        SigString::sigalign(fastq_path, read_layout, outbase.string(),
                            max_verbose, nthreads, -1);
        auto sig_time = std::chrono::steady_clock::now() - sigalign_start;
        std::cout << "[sigalign] Time: "
                  << std::chrono::duration_cast<std::chrono::seconds>(sig_time).count()
                  << " s\n";

        // Save whitelist summary
        if (verbose) std::cout << "[main] Saving whitelist summary...\n";
        read_layout.save_wl(outbase.string() + "_whitelist.csv", false);

        auto final_elapsed = std::chrono::steady_clock::now() - main_start;
        std::cout << "[main] Completed in "
                  << std::chrono::duration_cast<std::chrono::seconds>(final_elapsed).count()
                  << " s\n";
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what()
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n";
        return 1;
    }

    return 0;
}
