// main.cpp
#include "include/rad/rad_headers.h"

static void usage(const char* prog) {
    std::cerr
      << "Usage: " << prog
      << " -l LAYOUT -q FASTQ [-k KIT] [-w PATH] [-o PREFIX] [-d DIR] [-t THREADS] [-v VERBOSE] [-V MAX_VERBOSE]\n\n"
      << "  -l, --layout     layout key (five_prime|three_prime|splitseq) or a custom path\n"
      << "  -q, --fastq      input FASTQ file\n"
      << "  -k, --kit        use this kit's default whitelist (optional)\n"
      << "  -w, --whitelist  path to a single whitelist CSV (overrides kits)\n"
      << "  -m, --mutation   mutation space of each barcode (look for barcodes up to x edit distance away)\n"
      << "  -s, --shift      shift space of each barcode (shift barcodes x bases to the left or right)\n"
      << "  -o, --output     filename prefix (default: output)\n"
      << "  -d, --dir        output directory (default: $HOME/Desktop)\n"
      << "  -b, --bc_split   write reads into per-barcode FASTQ files\n"
      << "  -t, --threads    number of threads (default: 1)\n"
      << "  -v, --verbose    verbose mode\n"
      << "  -V, --max_verbose max verbose level (debug only, generates a lot of text)\n"
      << "  -h, --help        prints this menu and nothing else\n";
}

int main(int argc, char* argv[]) {
    std::string layout_key, fastq_path, kit, whitelist_path, output_prefix, output_dir;
    std::optional<int> shift, mut;
    int num_threads = 1;
    bool verbose, max_verbose = false;
    bool split_bc = false;

    const char* optstring = "l:q:k:w:m:s:o:d:t:vhw";
    struct option longopts[] = {
        {"layout",    required_argument, nullptr, 'l'},
        {"fastq",     required_argument, nullptr, 'q'},
        {"kit",       required_argument, nullptr, 'k'},
        {"whitelist", required_argument, nullptr, 'w'},
        {"mutation",  required_argument, nullptr, 'm'},
        {"shift",     required_argument, nullptr, 's'},
        {"output",    required_argument, nullptr, 'o'},
        {"dir",       required_argument, nullptr, 'd'},
        {"bc_split",  no_argument,       nullptr, 'b'},
        {"threads",   required_argument, nullptr, 't'},
        {"verbose",   no_argument,       nullptr, 'v'},
        {"max_verbose",no_argument,      nullptr, 'V'},
        {"help",      no_argument,       nullptr, 'h'},
        {nullptr,     0,                 nullptr,  0 }
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
          case 'l': layout_key     = optarg;            break;
          case 'q': fastq_path     = optarg;            break;
          case 'k': kit            = optarg;            break;
          case 'w': whitelist_path = optarg;            break;
          case 'm': mut            = std::stoi(optarg); break;
          case 's': shift          = std::stoi(optarg); break;
          case 'o': output_prefix  = optarg;            break;
          case 'b': split_bc       = true;              break;
          case 'd': output_dir     = optarg;            break;
          case 't': num_threads    = std::stoi(optarg); break;
          case 'v': verbose        = true;              break;
          case 'V': max_verbose    = true;              break;
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
                  << " layout key/path    : " << layout_key    << "\n"
                  << " FASTQ input        : " << fastq_path    << "\n"
                  << " kit                : " << (kit.empty() ? "[none specified]" : kit) << "\n"
                  << " whitelist          : " << (whitelist_path.empty() ? "[none specified]" : whitelist_path) << "\n"
                  << " output dir         : " << outdir        << "\n"
                  << " filename base      : " << base          << "\n"
                  << " threads            : " << num_threads   << "\n"
                  << " verbose            : " << std::boolalpha << verbose << "\n"
                  << " omp max thr        : " << omp_get_max_threads()
                  << "\n\n";
    }


    std::cout << "[main] Starting main processing...\n";

    try {
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

        if (have_layout) {
            if (verbose) std::cout << "[main] Importing read layout from " << outbase.string() + "_layout.csv" << "\n";
            read_layout.import_read_layout(outbase.string() + "_layout.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[main] Generating read layout...\n";
            read_layout.prep_new_layout(layout_csv, max_verbose);
        }

        if (have_map) {
            if (verbose) std::cout << "[main] Importing position map from " << outbase.string() + "_position_map.csv" << "\n";
            read_layout.import_position_map(outbase.string() + "_position_map.csv", max_verbose);
        } else {
            if (verbose) std::cout << "[main] Generating position map...\n";
            read_layout.generate_position_mapping();
        }

        //  compute misalignment if we regenerated either
        if (!have_layout || !have_map) {
            read_layout.write_to_csv(outbase.string(), "layout");
            if (verbose) std::cout << "[main] Computing misalignment...\n";
            Misalignment_Setup mis(read_layout);
            mis.generate_misalignment_data(fastq_path, read_layout, num_threads);

            if (verbose) std::cout << "[main] Writing layout and map...\n";
            read_layout.write_to_csv(outbase.string(), "both");
        }

        // wl loading
        if (!whitelist_path.empty()) {
            if (verbose) std::cout << "[main] Loading single whitelist: " << whitelist_path << "\n";
            read_layout.load_wl(shift, mut, max_verbose);
        }
        else if (!kit.empty()) {
            if (verbose) std::cout << "[main] Loading kit whitelist: " << kit << "\n";
            read_layout.load_wl(shift, mut, max_verbose);
        }
        else {
            if (verbose) std::cout << "[main] Loading all whitelists...\n";
            read_layout.load_wl(shift, mut, max_verbose);
        }

        // master demuxn fxn 
        if (verbose) std::cout << "[main] Running sigalign...\n";
        SigString::sigalign(fastq_path, read_layout, outbase.string(), /*verbose=*/max_verbose, num_threads, /*max_reads=*/-1);

        // save whitelist summary
        if (verbose) std::cout << "[main] Saving whitelist summary to " << outbase.string() + "_whitelist.csv" << "\n";
        read_layout.save_wl(outbase.string() + "_whitelist.csv", /*full=*/max_verbose);
    }
    catch (const std::exception &ex) {
        std::cerr << "[ERROR] " << ex.what() << "for " << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
        return 1;
    }

    return 0;
}
