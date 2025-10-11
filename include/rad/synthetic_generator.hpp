#include "rad_headers.h"

/**
 * @brief Simple synthetic read generator - assembles clean reads from layout + barcodes
 */
class SimpleSyntheticGenerator {
private:
    std::mt19937 rng;
    std::uniform_int_distribution<int> base_dist;
    const std::vector<char> bases = {'A', 'C', 'G', 'T'};

public:
    SimpleSyntheticGenerator(int seed = 42) : rng(seed), base_dist(0, 3) {}

    /**
     * @brief Generate synthetic reads from layout and barcode list
     */
    void generate_reads(const std::string& layout_csv_path,
                       const std::string& barcodes_csv_path, 
                       const std::string& output_prefix,
                       bool output_fastq = true,
                       bool verbose = false) {
        
        if (verbose) {
            std::cout << "[SimpleSynthetic] Loading layout from: " << layout_csv_path << std::endl;
            std::cout << "[SimpleSynthetic] Loading barcodes from: " << barcodes_csv_path << std::endl;
            std::cout << "[SimpleSynthetic] Output format: " << (output_fastq ? "FASTQ" : "FASTA") << std::endl;
        }

        // Load the read layout
        ReadLayout layout;
        layout.prep_new_layout(layout_csv_path, verbose);
        
        // Load barcodes
        auto barcodes = load_barcodes(barcodes_csv_path, verbose);
        
        if (barcodes.empty()) {
            throw std::runtime_error("No barcodes loaded from " + barcodes_csv_path);
        }

        // Generate reads
        generate_synthetic_reads(layout, barcodes, output_prefix, output_fastq, verbose);
    }

private:
    /**
     * @brief Load barcodes from CSV file
     */
    std::vector<std::string> load_barcodes(const std::string& csv_path, bool verbose) {
        std::vector<std::string> barcodes;
        
        try {
            auto lines = streaming_utils::import_text(csv_path);
            if (lines.empty()) {
                throw std::runtime_error("Empty barcode file");
            }

            // Skip header if it contains "barcode" or "Barcode"
            size_t start_idx = 0;
            if (!lines.empty() && (lines[0].find("barcode") != std::string::npos || 
                                  lines[0].find("Barcode") != std::string::npos)) {
                start_idx = 1;
            }

            for (size_t i = start_idx; i < lines.size(); ++i) {
                const auto& line = lines[i];
                if (line.empty()) continue;

                // Extract first column (barcode)
                auto comma_pos = line.find(',');
                std::string barcode = (comma_pos != std::string::npos) ? 
                                    line.substr(0, comma_pos) : line;
                
                // Trim whitespace
                barcode = seq_utils::trim(barcode);
                
                if (!barcode.empty()) {
                    barcodes.push_back(barcode);
                }
            }

            if (verbose) {
                std::cout << "[SimpleSynthetic] Loaded " << barcodes.size() 
                          << " barcodes" << std::endl;
            }

        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to load barcodes: " + std::string(e.what()));
        }

        return barcodes;
    }

    /**
     * @brief Structure to hold barcode and UMI information for a read
     */
    struct ReadBarcodeInfo {
        std::string barcode;
        std::string umi;
    };

    /**
     * @brief Generate synthetic reads for each barcode
     */
    void generate_synthetic_reads(const ReadLayout& layout,
                                 const std::vector<std::string>& barcodes,
                                 const std::string& output_prefix,
                                 bool output_fastq,
                                 bool verbose) {
        
        std::string file_extension = output_fastq ? ".fq" : ".fa";
        std::string reads_path = output_prefix + file_extension;
        std::string metadata_path = output_prefix + "_metadata.csv";
        
        std::ofstream reads_file(reads_path);
        std::ofstream metadata_file(metadata_path);
        
        if (!reads_file || !metadata_file) {
            throw std::runtime_error("Failed to open output files");
        }

        // Write metadata file header
        metadata_file << "read_id,original_barcode,original_umi\n";

        if (verbose) {
            std::cout << "[SimpleSynthetic] Generating " << barcodes.size() 
                      << " synthetic reads..." << std::endl;
        }

        // Generate one read per barcode
        for (size_t i = 0; i < barcodes.size(); ++i) {
            const std::string& barcode = barcodes[i];
            std::string read_id = "synthetic_" + std::to_string(i + 1);
            
            auto [sequence, barcode_info] = assemble_read(layout, barcode, read_id);
            
            // Write sequence file (FASTQ or FASTA)
            if (output_fastq) {
                // Generate perfect quality scores (all high quality)
                std::string quality_scores(sequence.length(), '5'); // Phred 20
                
                reads_file << "@" << read_id << "\n"
                          << sequence << "\n"
                          << "+\n" 
                          << quality_scores << "\n";
            } else {
                // FASTA format
                reads_file << ">" << read_id << "\n"
                          << sequence << "\n";
            }
            
            // Write metadata
            metadata_file << read_id << "," 
                         << barcode_info.barcode << ","
                         << barcode_info.umi << "\n";

            if (verbose && (i + 1) % 10000 == 0) {
                std::cout << "[SimpleSynthetic] Generated " << (i + 1) 
                          << " reads..." << std::endl;
            }
        }

        if (verbose) {
            std::cout << "[SimpleSynthetic] Generation complete!" << std::endl;
            std::cout << "[SimpleSynthetic] Files written:" << std::endl;
            std::cout << "  Reads: " << reads_path << std::endl;
            std::cout << "  Metadata: " << metadata_path << std::endl;
        }
    }

    /**
     * @brief Assemble a single read according to the layout
     */
    std::pair<std::string, ReadBarcodeInfo> 
    assemble_read(const ReadLayout& layout, const std::string& barcode, const std::string& read_id) {
        
        std::string full_sequence;
        ReadBarcodeInfo barcode_info;
        
        // Process elements in order - we'll handle both forward and reverse if present
        std::vector<const ReadElement*> ordered_elements;
        for (const auto& elem : layout.by_order()) {
            ordered_elements.push_back(&elem);
        }

        for (const auto* elem : ordered_elements) {
            std::string elem_sequence = generate_element_sequence(*elem, barcode, barcode_info);
            
            // Add to full sequence
            full_sequence += elem_sequence;
        }

        return {full_sequence, barcode_info};
    }

    /**
     * @brief Generate sequence for a single element
     */
    std::string generate_element_sequence(const ReadElement& elem, const std::string& barcode, ReadBarcodeInfo& barcode_info) {
        
        if (elem.global_class == "start" || elem.global_class == "stop") {
            return ""; // These are just positional markers
        }
        
        if (elem.type == "static") {
            if (elem.global_class == "poly_tail") {
                // Generate poly-A or poly-T based on the sequence pattern
                if (elem.seq.find('A') != std::string::npos || elem.class_id.find("poly_a") != std::string::npos) {
                    return std::string(elem.expected_length.value_or(16), 'A');
                } else {
                    return std::string(elem.expected_length.value_or(16), 'T');
                }
            } else {
                // Regular static sequence
                return elem.seq;
            }
        }
        
        if (elem.global_class == "barcode") {
            std::string barcode_seq;
            // Use the provided barcode, potentially reverse complement for reverse elements
            if (elem.direction == "reverse" && elem.class_id.find("rc_") == 0) {
                barcode_seq = seq_utils::revcomp(barcode);
            } else {
                barcode_seq = barcode;
            }
            
            // Store the original barcode (only store once, prefer forward direction)
            if (barcode_info.barcode.empty() || elem.direction == "forward") {
                barcode_info.barcode = barcode;
            }
            
            return barcode_seq;
        }
        
        if (elem.global_class == "umi") {
            // Generate random UMI sequence
            int length = elem.expected_length.value_or(10);
            std::string umi_seq = generate_random_sequence(length);
            
            // Store the original UMI (only store once, prefer forward direction)
            if (barcode_info.umi.empty() || elem.direction == "forward") {
                barcode_info.umi = umi_seq;
            }
            
            return umi_seq;
        }
        
        if (elem.global_class == "read") {
            // Generate random read sequence
            int length = elem.expected_length.value_or(400);
            return generate_random_sequence(length);
        }
        
        // Default: generate random sequence for variable elements
        int length = elem.expected_length.value_or(20);
        return generate_random_sequence(length);
    }

    /**
     * @brief Generate random DNA sequence
     */
    std::string generate_random_sequence(int length, double gc_content = 0.5) {
        std::string sequence;
        sequence.reserve(length);
        
        //use base_dist for equal probability, could be enhanced for GC content
        for (int i = 0; i < length; ++i) {
            sequence += bases[base_dist(rng)];
        }
        
        return sequence;
    }
};

/**
 * @brief Command-line interface for simple synthetic generation
 */
static void usage_synthetic(const char* prog) {
    std::cerr
      << "Usage: " << prog << " [options]\n"
      << "Generate clean synthetic reads from layout specification and barcode list\n\n"
      << "Required arguments:\n"
      << "  -l, --layout FILE             Read layout CSV file\n"
      << "  -b, --barcodes FILE           Barcode CSV file (one barcode per read)\n"
      << "  -o, --output PREFIX           Output file prefix\n\n"
      << "Optional arguments:\n"
      << "  -s, --seed SEED               Random seed for UMIs/reads (default: 42)\n"
      << "  -f, --format FORMAT           Output format: fastq or fasta (default: fastq)\n"
      << "  -v, --verbose                 Verbose output\n"
      << "  -h, --help                    Show this help message\n\n"
      << "Example:\n"
      << "  " << prog << " -l three_prime_layout.csv -b my_barcodes.csv -o synthetic_data -f fasta -v\n\n"
      << "Output files:\n"
      << "  PREFIX.fq/.fa         - Synthetic reads (FASTQ/FASTA format)\n"
      << "  PREFIX_metadata.csv   - Metadata with read_id, original_barcode, original_umi\n\n";
}

int main_synthetic_simple(int argc, char* argv[]) {
    std::string layout_path, barcodes_path, output_prefix, format = "fastq";
    int seed = 42;
    bool verbose = false;
    
    const char* optstring = "l:b:o:s:f:vh";
    struct option longopts[] = {
        {"layout",    required_argument, nullptr, 'l'},
        {"barcodes",  required_argument, nullptr, 'b'},
        {"output",    required_argument, nullptr, 'o'},
        {"seed",      required_argument, nullptr, 's'},
        {"format",    required_argument, nullptr, 'f'},
        {"verbose",   no_argument,       nullptr, 'v'},
        {"help",      no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, optstring, longopts, nullptr)) != -1) {
        switch (c) {
            case 'l': layout_path = optarg; break;
            case 'b': barcodes_path = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 's': seed = std::stoi(optarg); break;
            case 'f': format = optarg; break;
            case 'v': verbose = true; break;
            case 'h': usage_synthetic(argv[0]); return 0;
            default:
                std::cerr << "Use -h for help\n";
                return 1;
        }
    }
    
    if (layout_path.empty() || barcodes_path.empty() || output_prefix.empty()) {
        std::cerr << "Error: Layout file, barcodes file, and output prefix are all required\n";
        usage_synthetic(argv[0]);
        return 1;
    }
    
    // Validate format
    bool output_fastq = true;
    if (format == "fasta" || format == "fa") {
        output_fastq = false;
    } else if (format != "fastq" && format != "fq") {
        std::cerr << "Error: Format must be 'fastq' or 'fasta'\n";
        return 1;
    }
    
    try {
        SimpleSyntheticGenerator generator(seed);
        generator.generate_reads(layout_path, barcodes_path, output_prefix, output_fastq, verbose);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
}