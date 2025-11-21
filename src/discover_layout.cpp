#include "include/rad/rad_headers.h"

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --parts parts.csv --fastq reads.fastq --output PREFIX [options]\n"
              << "\n"
              << "Required arguments:\n"
              << "  --parts, -p    CSV with at least 'id' and 'seq' columns (static parts list)\n"
              << "  --fastq, -q    FASTQ/FASTA file to sample reads from\n"
              << "  --output, -o   Prefix for generated <prefix>_layout.csv\n"
              << "\n"
              << "Options:\n"
              << "  --max_reads, -n    Maximum reads to sample (default: 25000)\n"
              << "  --chunk_size, -c   Reads per streaming chunk (default: 2000)\n"
              << "  --threads, -t      Number of threads for streaming (default: 1)\n"
              << "  --verbose, -v      Verbose logging\n"
              << "  --help, -h         Show this message\n";
}

int main(int argc, char* argv[]) {
    std::string parts_csv;
    std::string fastq_path;
    std::string output_prefix;
    size_t max_reads = 25000;
    size_t chunk_size = 2000;
    int threads = 1;
    bool verbose = false;

    const char* optstring = "p:q:o:n:c:t:vh";
    struct option longopts[] = {
        {"parts", required_argument, nullptr, 'p'},
        {"fastq", required_argument, nullptr, 'q'},
        {"output", required_argument, nullptr, 'o'},
        {"max_reads", required_argument, nullptr, 'n'},
        {"chunk_size", required_argument, nullptr, 'c'},
        {"threads", required_argument, nullptr, 't'},
        {"verbose", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'p': parts_csv = optarg; break;
            case 'q': fastq_path = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 'n': max_reads = std::stoull(optarg); break;
            case 'c': chunk_size = std::stoull(optarg); break;
            case 't': threads = std::stoi(optarg); break;
            case 'v': verbose = true; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (parts_csv.empty() || fastq_path.empty() || output_prefix.empty()) {
        std::cerr << "[ERROR] --parts, --fastq, and --output are required\n";
        usage(argv[0]);
        return 1;
    }

    try {
        ReadLayout rl;
        rl.discover_layout(parts_csv, fastq_path, output_prefix, max_reads, chunk_size, threads, verbose);
        if (verbose) {
            std::cout << "[generate_layout_from_parts] Wrote layout to '" << output_prefix << "_layout.csv'\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "[ERROR] " << ex.what() << "\n";
        return 1;
    }

    return 0;
}