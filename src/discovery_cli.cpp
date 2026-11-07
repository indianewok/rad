// Included by main.cpp: RAD's existing header implementation uses one main TU.
#include "include/rad/discovery.hpp"

#include <limits>
#include <stdexcept>

namespace rad_discovery {
namespace cli_detail {

void usage() {
    std::cerr
        << "Usage: rad discover --wl -i INPUT -p PRIMER -n LENGTH -o PREFIX [options]\n\n"
        << "Identify a nonspatial whitelist from a known flanking primer or masked pattern.\n"
        << "Reuse scan-wl raw extraction, count exact read hits independently in every\n"
        << "eligible catalog, and return the highest count. Equal top counts are ties.\n"
        << "Spatial catalog identification is unsupported.\n\n"
        << "  -i, --input PATH           FASTQ/FASTA input (.gz supported)\n"
        << "  -o, --output PREFIX        output prefix (existing outputs are not overwritten)\n"
        << "      --layout               write inferred layout and parts only\n"
        << "      --wl                   identify whitelist only; requires --pattern\n"
        << "                              supplying both flags requests both outputs\n"
        << "      --parts CSV            additional static parts (id,seq columns)\n"
        << "      --sample-reads N       prefix sample limit in reads (default: 10000)\n"
        << "      --sample-bases N       prefix sample limit in bases (default: 50000000)\n"
        << "  -t, --threads N            worker threads (default: 1)\n"
        << "  -p, --pattern SEQUENCE     known flanking primer or masked pattern (--wl only)\n"
        << "                              N captures; n matches without capture\n"
        << "  -n, --barcode-length N     barcode length for an unmasked pattern (1-64)\n"
        << "  -h, --help                 show this help\n\n"
        << "Writes PREFIX_whitelist_match.csv (winner, ties or no hits),\n"
        << "PREFIX_whitelists.csv (full ranking), and PREFIX_discover.log.\n"
        << "No minimum distinct-barcode count, fraction or catalog-size penalty is applied.\n\n"
        << "For compatibility, no mode flags requests both layout and whitelist using\n"
        << "inferred extraction recipes; --layout alone writes layout and parts only.\n"
        << "Parts output is PREFIX_parts.csv. A complete layout is PREFIX_layout.csv; an unresolved\n"
        << "structure is PREFIX_layout.partial.csv and is not a runnable layout.\n"
        << "Exit codes: 0 each requested whitelist recipe has a unique nonzero winner\n"
        << "and any requested layout is complete; 2 tied/no-hit match or incomplete layout;\n"
        << "1 invalid arguments, input or processing error.\n";
}

std::size_t positive_size(const char* value, const std::string& name) {
    const std::string text = value ? value : "";
    if (text.empty() || !std::all_of(text.begin(), text.end(), [](unsigned char c) {
            return c >= '0' && c <= '9';
        })) {
        throw std::invalid_argument(name + " must be a positive integer");
    }
    std::size_t used = 0;
    unsigned long long number = 0;
    try { number = std::stoull(text, &used); }
    catch (const std::exception&) {
        throw std::invalid_argument(name + " is outside the supported range");
    }
    if (used != text.size() || number == 0 ||
        number > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(name + " must be a positive integer in range");
    }
    return static_cast<std::size_t>(number);
}

int positive_int(const char* value, const std::string& name, int maximum) {
    const auto number = positive_size(value, name);
    if (number > static_cast<std::size_t>(maximum)) {
        throw std::invalid_argument(name + " exceeds " + std::to_string(maximum));
    }
    return static_cast<int>(number);
}

std::string csv_string(const std::string& value) {
    std::string result = "\"";
    for (char c : value) {
        if (c == '"') result += '"';
        result += c;
    }
    return result + '"';
}

std::string join_aliases(const std::vector<std::string>& aliases) {
    std::string result;
    for (const auto& alias : aliases) {
        if (!result.empty()) result += ';';
        result += alias;
    }
    return result;
}

void ensure_new_output(const std::string& path) {
    namespace fs = boost::filesystem;
    if (fs::exists(path) || fs::is_symlink(fs::symlink_status(path))) {
        throw std::runtime_error("output already exists: " + path + "; choose a new prefix");
    }
    const auto parent = fs::path(path).parent_path();
    if (!parent.empty() && !fs::is_directory(parent)) {
        throw std::runtime_error("output directory does not exist: " + parent.string());
    }
}

std::ofstream output_file(const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write output: " + path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    return out;
}

struct TemporarySample {
    boost::filesystem::path directory;
    TemporarySample() {
        namespace fs = boost::filesystem;
        const auto parent = fs::temp_directory_path();
        for (int attempt = 0; attempt < 10; ++attempt) {
            auto candidate = parent / fs::unique_path("rad-discover-%%%%-%%%%-%%%%");
            if (fs::create_directory(candidate)) {
                directory = std::move(candidate);
                try {
                    fs::permissions(directory, fs::owner_all);
                } catch (...) {
                    boost::system::error_code ignored;
                    fs::remove_all(directory, ignored);
                    directory.clear();
                    throw;
                }
                return;
            }
        }
        throw std::runtime_error("cannot create temporary discovery sample directory");
    }
    ~TemporarySample() {
        if (!directory.empty()) {
            boost::system::error_code ignored;
            boost::filesystem::remove_all(directory, ignored);
        }
    }
    TemporarySample(const TemporarySample&) = delete;
    TemporarySample& operator=(const TemporarySample&) = delete;
};

std::vector<Read> sample_reads(const Options& options, std::size_t& bases) {
    // zlib keeps prefix sampling bounded without a decompressor helper process.
    file_streaming files(options.input, options.threads, false);
    read_streaming reader(files);
    std::vector<Read> reads;
    reads.reserve(std::min<std::size_t>(options.sample_reads, 10000));
    bases = 0;
    while (reads.size() < options.sample_reads && bases < options.sample_bases) {
        auto record = reader.next_sequence();
        if (!record) break;
        if (record->seq.size() > options.sample_bases - bases) break;
        if (record->seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("sampled read exceeds the supported alignment length");
        }
        std::transform(record->seq.begin(), record->seq.end(), record->seq.begin(), [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
        bases += record->seq.size();
        reads.push_back({std::move(record->id), std::move(record->seq), std::move(record->qual)});
    }
    if (reader.read_failed()) {
        throw std::runtime_error("malformed sequence input in " + reader.read_error_path() +
                                 " (reader error " + std::to_string(reader.read_error_code()) + ")");
    }
    if (reads.empty()) {
        throw std::runtime_error("sample contains no nonempty reads within the base limit; check input or increase --sample-bases");
    }
    return reads;
}

std::string write_sample(const TemporarySample& temporary, const std::vector<Read>& reads) {
    const bool fastq = std::all_of(reads.begin(), reads.end(), [](const Read& read) {
        return read.quality.size() == read.sequence.size();
    });
    const auto path = (temporary.directory / (fastq ? "sample.fastq" : "sample.fasta")).string();
    auto out = output_file(path);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        const auto& read = reads[i];
        // Sample-local identifiers avoid header ambiguities; the cached Read retains the original ID.
        if (fastq) out << "@sample_" << i << '\n' << read.sequence << "\n+\n" << read.quality << '\n';
        else out << ">sample_" << i << '\n' << read.sequence << '\n';
    }
    out.close();
    return path;
}

void write_parts(const std::string& path, const LayoutInference& inference) {
    auto out = output_file(path);
    out << "id,sequence,reads,forward,reverse\n";
    for (const auto& part : inference.parts) {
        out << csv_string(part.id) << ',' << csv_string(part.sequence) << ','
            << part.reads << ',' << part.forward << ',' << part.reverse << '\n';
    }
    out.close();
}

void write_matches(const std::string& path, const std::vector<CatalogMatch>& matches) {
    auto out = output_file(path);
    out << "recipe_id,catalog,aliases,length,catalog_size,matched_reads,matched_distinct,extracted_reads,observed_distinct,sampled_reads\n";
    for (const auto& match : matches) {
        out << csv_string(match.recipe_id) << ',' << csv_string(match.path) << ','
            << csv_string(join_aliases(match.aliases)) << ',' << match.length << ','
            << match.catalog_size << ',' << match.reads << ',' << match.distinct << ','
            << match.extracted << ',' << match.observed_distinct << ',' << match.sampled_reads << '\n';
    }
    out.close();
}

void write_choices(const std::string& path, const std::vector<WhitelistChoice>& choices,
                   const std::vector<CatalogMatch>& matches) {
    auto out = output_file(path);
    out << "status,recipe_id,catalog,aliases,matched_reads,extracted_reads,matched_distinct\n";
    if (choices.empty()) out << "no_hits,\"\",\"\",\"\",0,,0\n";
    for (const auto& choice : choices) {
        if (choice.winners.empty()) {
            // If multiple geometries for this recipe have different extraction
            // totals, no single denominator exists. An absent eligible catalog
            // also leaves this denominator unknown rather than establishing 0.
            std::optional<std::uint64_t> extracted;
            bool ambiguous = false;
            for (const auto& match : matches) {
                if (match.recipe_id != choice.recipe_id) continue;
                if (!extracted) extracted = match.extracted;
                else if (*extracted != match.extracted) ambiguous = true;
            }
            out << "no_hits," << csv_string(choice.recipe_id) << ",\"\",\"\",0,";
            if (extracted && !ambiguous) out << *extracted;
            out << ",0\n";
            continue;
        }
        for (const auto& match : choice.winners) {
            out << choice.status << ',' << csv_string(choice.recipe_id) << ','
                << csv_string(match.path) << ',' << csv_string(join_aliases(match.aliases)) << ','
                << match.reads << ',' << match.extracted << ',' << match.distinct << '\n';
        }
    }
    out.close();
}

Recipe explicit_recipe(Options& options, bool length_given) {
    std::size_t capture = 0, literal = 0;
    for (char& c : options.pattern) {
        if (c == 'N') { ++capture; continue; }
        if (c == 'n') continue;
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            throw std::invalid_argument("--pattern accepts A, C, G, T, N (capture), and n (match-only)");
        }
        ++literal;
    }
    if (literal == 0) throw std::invalid_argument("--pattern needs at least one literal anchor base");
    if (capture > 64) throw std::invalid_argument("--pattern captures more than 64 barcode bases");
    if (capture > 0) {
        if (length_given && options.barcode_length != static_cast<int>(capture)) {
            throw std::invalid_argument("--barcode-length must match the uppercase N capture count");
        }
        options.barcode_length = static_cast<int>(capture);
    } else if (!length_given) {
        throw std::invalid_argument("an unmasked --pattern requires --barcode-length");
    }
    return {"explicit_pattern", options.pattern, options.barcode_length, 0, 0};
}

} // namespace cli_detail

int run(int argc, char** argv) {
    using namespace cli_detail;
    Options options;
    bool layout_flag = false, wl_flag = false, length_given = false, pattern_given = false;
    std::ofstream log;
    try {
        const option long_options[] = {
            {"input", required_argument, nullptr, 'i'},
            {"output", required_argument, nullptr, 'o'},
            {"threads", required_argument, nullptr, 't'},
            {"layout", no_argument, nullptr, 1000},
            {"wl", no_argument, nullptr, 1001},
            {"sample-reads", required_argument, nullptr, 1002},
            {"sample-bases", required_argument, nullptr, 1003},
            {"parts", required_argument, nullptr, 1004},
            {"pattern", required_argument, nullptr, 'p'},
            {"barcode-length", required_argument, nullptr, 'n'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}
        };
        optind = 1;
        opterr = 0;
        int choice;
        while ((choice = getopt_long(argc, argv, "i:o:t:p:n:h", long_options, nullptr)) != -1) {
            switch (choice) {
                case 'i': options.input = optarg; break;
                case 'o': options.output = optarg; break;
                case 't': options.threads = positive_int(optarg, "--threads", std::numeric_limits<int>::max()); break;
                case 'p': options.pattern = optarg; pattern_given = true; break;
                case 'n': options.barcode_length = positive_int(optarg, "--barcode-length", 64); length_given = true; break;
                case 'h': usage(); return 0;
                case 1000: layout_flag = true; break;
                case 1001: wl_flag = true; break;
                case 1002: options.sample_reads = positive_size(optarg, "--sample-reads"); break;
                case 1003: options.sample_bases = positive_size(optarg, "--sample-bases"); break;
                case 1004: options.parts_csv = optarg; break;
                default: throw std::invalid_argument("unknown option or missing value; see 'rad discover --help'");
            }
        }
        if (optind != argc) throw std::invalid_argument("unexpected positional argument: " + std::string(argv[optind]));
        if (layout_flag || wl_flag) { options.layout = layout_flag; options.wl = wl_flag; }
        if (options.input.empty() || options.output.empty()) {
            throw std::invalid_argument("--input and --output are required");
        }
        if (!boost::filesystem::is_regular_file(options.input)) {
            throw std::invalid_argument("input is not a readable sequence file: " + options.input);
        }
        if (!options.parts_csv.empty() && !boost::filesystem::is_regular_file(options.parts_csv)) {
            throw std::invalid_argument("parts file does not exist: " + options.parts_csv);
        }
        if (pattern_given && options.pattern.empty()) {
            throw std::invalid_argument("--pattern must not be empty");
        }
        if (length_given && options.pattern.empty()) {
            throw std::invalid_argument("--barcode-length requires --pattern");
        }
        if (!options.pattern.empty() && (options.layout || !options.wl)) {
            throw std::invalid_argument("--pattern is available with --wl only; omit it to infer structure and extraction recipes");
        }
        if (options.wl && !options.layout && options.pattern.empty()) {
            throw std::invalid_argument("--wl requires --pattern with a known flanking primer or masked capture pattern");
        }
        std::vector<Recipe> manual_recipes;
        if (!options.pattern.empty()) manual_recipes.push_back(explicit_recipe(options, length_given));
        const bool infer = manual_recipes.empty();
        const std::string log_path = options.output + "_discover.log";
        ensure_new_output(log_path);
        if (infer) ensure_new_output(options.output + "_parts.csv");
        if (options.wl) {
            ensure_new_output(options.output + "_whitelists.csv");
            ensure_new_output(options.output + "_whitelist_match.csv");
        }
        if (options.layout) {
            ensure_new_output(options.output + "_layout.csv");
            ensure_new_output(options.output + "_layout.partial.csv");
        }
        log = output_file(log_path);
        log << "command=discover\ninput=" << options.input << "\nmode_layout=" << options.layout
            << "\nmode_wl=" << options.wl << "\nsampling=prefix_whole_reads\n"
            << "sample_reads_limit=" << options.sample_reads << "\nsample_bases_limit=" << options.sample_bases
            << "\nthreads=" << options.threads << '\n';
        if (options.wl) {
            log << "whitelist_scope=nonspatial_catalogs_only\n"
                << "whitelist_ranking=highest_independent_exact_read_hits\n"
                << "note=Spatial catalog identification is unsupported. No distinct-count or catalog-size penalty is applied.\n";
        }
        std::size_t bases = 0;
        auto reads = sample_reads(options, bases);
        log << "sample_reads=" << reads.size() << "\nsample_bases=" << bases << '\n';
        std::cout << "[discover] Sampled " << reads.size() << " reads (" << bases << " bases).\n";
        TemporarySample temporary;
        const auto sample_path = write_sample(temporary, reads);
        LayoutInference inference;
        if (infer) {
            inference = infer_layout(reads, options);
            write_parts(options.output + "_parts.csv", inference);
            for (const auto& note : inference.notes) log << "note=" << note << '\n';
            log << "layout_complete=" << inference.complete << "\nobserved_parts=" << inference.parts.size() << '\n';
        } else {
            inference.recipes = std::move(manual_recipes);
            log << "extraction=explicit_pattern\npattern=" << options.pattern
                << "\nbarcode_length=" << options.barcode_length << '\n';
        }
        log << "extraction_recipes=" << inference.recipes.size() << '\n';
        std::vector<CatalogMatch> matches;
        if (options.wl) {
            if (!inference.recipes.empty()) matches = identify_whitelists(sample_path, inference.recipes, options);
            write_matches(options.output + "_whitelists.csv", matches);
        }
        const auto exact = std::count_if(matches.begin(), matches.end(), [](const CatalogMatch& match) {
            return match.reads > 0 && match.distinct > 0;
        });
        std::vector<WhitelistChoice> choices;
        bool whitelist_resolved = true;
        if (options.wl) {
            choices = select_whitelists(matches, inference.recipes);
            write_choices(options.output + "_whitelist_match.csv", choices, matches);
            whitelist_resolved = !choices.empty() && std::all_of(choices.begin(), choices.end(), [](const WhitelistChoice& choice) {
                return choice.status == "matched";
            });
            const bool any_tie = std::any_of(choices.begin(), choices.end(), [](const WhitelistChoice& choice) {
                return choice.status == "tie";
            });
            log << "catalog_rows=" << matches.size() << "\ncatalog_rows_with_exact_matches=" << exact
                << "\nwhitelist_status=" << (whitelist_resolved ? "matched" : any_tie ? "tie" : "no_hits") << '\n';
            log << "note=Winner selection uses exact matching read count only. Catalog aliases do not distinguish chemistry versions.\n";
            for (const auto& choice : choices) {
                log << "whitelist_status." << choice.recipe_id << '=' << choice.status << '\n'
                    << "top_read_hits." << choice.recipe_id << '=' << choice.matched_reads << '\n';
                if (choices.size() == 1) log << "top_read_hits=" << choice.matched_reads << '\n';
                if (choice.status == "matched") {
                    const auto& winner = choice.winners.front();
                    log << "selected_whitelist." << choice.recipe_id << '=' << winner.path << '\n';
                    if (choices.size() == 1) log << "selected_whitelist=" << winner.path << '\n';
                    std::cout << "[discover] Selected whitelist for " << choice.recipe_id << ": " << winner.path
                              << " (" << join_aliases(winner.aliases) << "); " << winner.reads << " exact read hits.\n";
                } else if (choice.status == "tie") {
                    log << "selected_whitelist." << choice.recipe_id << "=\n";
                    if (choices.size() == 1) log << "selected_whitelist=\n";
                    std::cout << "[discover] Tied whitelist catalogs for " << choice.recipe_id << ": "
                              << choice.matched_reads << " exact read hits each.\n";
                    for (const auto& winner : choice.winners) {
                        log << "tied_whitelist." << choice.recipe_id << '=' << winner.path << '\n';
                        std::cout << "[discover]   " << winner.path << " (" << join_aliases(winner.aliases) << ")\n";
                    }
                } else {
                    log << "selected_whitelist." << choice.recipe_id << "=\n";
                    if (choices.size() == 1) log << "selected_whitelist=\n";
                    std::cout << "[discover] No exact read hits in eligible nonspatial whitelist catalogs for "
                              << choice.recipe_id << ".\n";
                }
            }
            if (choices.empty()) std::cout << "[discover] No whitelist extraction recipe is available.\n";
        }
        if (options.layout) {
            const auto path = options.output + (inference.complete ? "_layout.csv" : "_layout.partial.csv");
            write_layout(path, inference, matches);
            log << "layout_output=" << path << '\n';
            std::cout << "[discover] " << (inference.complete ? "Layout: " : "Partial structure (not runnable): ") << path << '\n';
        }
        const bool incomplete = (options.layout && !inference.complete) || (options.wl && !whitelist_resolved);
        log << "status=" << (incomplete ? "incomplete" : "complete") << '\n';
        log.close();
        return incomplete ? 2 : 0;
    } catch (const std::exception& error) {
        if (log.is_open()) {
            try { log << "status=error\nerror=" << error.what() << '\n'; log.close(); }
            catch (...) {}
        }
        std::cerr << "[discover] Error: " << error.what() << '\n';
        return 1;
    }
}

} // namespace rad_discovery
