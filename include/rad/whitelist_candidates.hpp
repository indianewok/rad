#pragma once

// Reference-guided candidates are a superset, not abundance-selected calls.
namespace rad_whitelist_candidates {

struct Options {
    bool enabled = false;
    std::string layout, position_map;
    int mismatches = 2;
    bool mismatches_explicit = false, concat_hmm = false;
};

struct temporary_outputs {
    std::vector<boost::filesystem::path> paths;
    ~temporary_outputs() {
        for (const auto& path : paths) {
            boost::system::error_code ec;
            boost::filesystem::remove(path, ec);
        }
    }
};

inline bool parse_nonnegative_int(const char* text, int& value) {
    if (!text || !*text) return false;
    const std::string token(text);
    if (!std::all_of(token.begin(), token.end(),
            [](unsigned char c) { return std::isdigit(c); })) return false;
    try {
        size_t consumed = 0;
        const auto parsed = std::stoull(token, &consumed);
        if (consumed != token.size() || parsed > static_cast<unsigned long long>(std::numeric_limits<int>::max())) return false;
        value = static_cast<int>(parsed);
        return true;
    } catch (...) { return false; }
}

inline std::string unpack(whitelist_relevance::packed_barcode barcode) {
    std::string sequence(barcode.length, 'A');
    for (size_t i = sequence.size(); i > 0; --i) {
        sequence[i - 1] = "ACTG"[barcode.bits & 3U];
        barcode.bits >>= 2;
    }
    return sequence;
}

inline bool parse_catalog_token(const std::string& token, uint16_t numeric_length,
        whitelist_relevance::packed_barcode& barcode) {
    if (whitelist_relevance::pack(token, barcode)) return true;
    if (token.empty() || numeric_length == 0 || numeric_length > 32) return false;
    const size_t first = token.front() == '-' ? 1 : 0;
    if (first == token.size() || !std::all_of(token.begin() + first, token.end(),
            [](unsigned char c) { return std::isdigit(c); })) return false;
    try {
        size_t consumed = 0;
        if (first && numeric_length != 32) return false;
        const uint64_t bits = first
            ? static_cast<uint64_t>(std::stoll(token, &consumed))
            : std::stoull(token, &consumed);
        if (consumed != token.size() || (numeric_length < 32 &&
                bits >= (uint64_t{1} << (2 * numeric_length)))) return false;
        barcode = {bits, numeric_length};
        return true;
    } catch (...) { return false; }
}

inline int run(const Options& options, const std::string& input,
        const std::string& prefix, const std::string& catalog_spec,
        int max_reads, int num_threads, int chunk_size, bool verbose) {
    const auto started = std::chrono::steady_clock::now();
    try {
        if (input.empty() || prefix.empty() || catalog_spec.empty() || options.layout.empty() || options.position_map.empty())
            throw std::runtime_error("--catalog-candidates requires -i, -o, -w, --layout and --position-map");
        if (max_reads < 0 || num_threads < 0 || chunk_size <= 0 || options.mismatches < 0)
            throw std::runtime_error("candidate generation requires nonnegative read/thread/mismatch limits and a positive chunk size");
        const auto specs = whitelist_utils::parse_whitelist_specs(catalog_spec);
        if (specs.size() != 1) throw std::runtime_error("candidate generation requires one catalogue without dual-whitelist syntax");
        const std::string catalog_path = whitelist_utils::kit_to_path(specs.front());
        const std::vector<std::string> sources{input, options.layout, options.position_map, catalog_path};
        for (const auto& source : sources)
            if (!boost::filesystem::is_regular_file(source))
                throw std::runtime_error("candidate generation needs a regular, replayable file: " + source);
        const std::string text_path = prefix + "_candidates.txt", log_path = prefix + "_catalog_candidates.log";
        for (const auto& output : {text_path, log_path}) {
            if (boost::filesystem::is_directory(output))
                throw std::runtime_error("candidate output is a directory: " + output);
            for (const auto& source : sources) {
                boost::system::error_code ec;
                if (boost::filesystem::equivalent(output, source, ec) && !ec)
                    throw std::runtime_error("candidate output would overwrite an input: " + output);
            }
        }
        if (!ReadLayout::is_generated_cache(options.layout))
            throw std::runtime_error("--layout must be a prepared RAD layout CSV; use rad prep first");
        ReadLayout layout;
        layout.import_read_layout(options.layout, verbose, options.layout);
        layout.import_position_map(options.position_map, verbose);
        std::set<std::string> classes;
        std::set<uint16_t> lengths;
        for (const auto& element : layout.layout) {
            if (element.global_class != "barcode" || element.type != "variable") continue;
            classes.insert(seq_utils::remove_rc(element.class_id));
            std::string flags = element.flags;
            std::transform(flags.begin(), flags.end(), flags.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (flags.find("joint_barcode") != std::string::npos)
                throw std::runtime_error("catalogue candidate generation supports non-joint barcode layouts only");
            if (!element.expected_length || *element.expected_length <= 0 || *element.expected_length > 32 ||
                    std::any_of(element.length_candidates.begin(), element.length_candidates.end(),
                        [](int length) { return length <= 0 || length > 32; }))
                throw std::runtime_error("catalogue candidate generation supports barcode lengths 1..32 only");
            lengths.insert(static_cast<uint16_t>(*element.expected_length));
            if (layout.position_map.find(element.class_id) == layout.position_map.end())
                throw std::runtime_error("prepared position map is missing barcode coordinates: " + element.class_id);
        }
        if (classes.size() != 1 || lengths.size() != 1)
            throw std::runtime_error("catalogue candidate generation requires one barcode class with a consistent expected length");
        if (options.concat_hmm) {
            if (layout.concat_params.empty())
                throw std::runtime_error("--concat-hmm requires frozen hmm_* parameters in --position-map; candidate generation does not calibrate");
            if (!layout.build_concat_model(verbose))
                throw std::runtime_error("cannot build the candidate-generation HMM from the prepared layout and position map");
        }
        const int threads = num_threads > 0 ? num_threads : std::max(1, omp_get_max_threads());
        const auto query_started = std::chrono::steady_clock::now();
        std::cout << "[scan-wl candidates] Collecting pre-correction queries without abundance or spatial filtering\n";
        auto observations = SigString::observe_whitelist_queries<whitelist_relevance::observation_map>(
            input, layout, threads, static_cast<size_t>(chunk_size), static_cast<size_t>(max_reads), 0);
        const double query_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - query_started).count();
        auto& observed = observations[*classes.begin()];
        if (observed.unsupported_long_queries)
            throw std::runtime_error("observed barcode queries longer than 32 bases are unsupported for candidate generation");
        whitelist_relevance::candidate_index index(observed, options.mismatches);
        whitelist_relevance::barcode_set retained;
        uint64_t rows = 0, valid_rows = 0, skipped_rows = 0;
        const auto catalog_started = std::chrono::steady_clock::now();
        // zlib transparently reads plain files too, and exposes truncation/CRC
        // errors that an istream wrapper can otherwise disguise as normal EOF.
        std::unique_ptr<gzFile_s, decltype(&gzclose)> catalog(
            gzopen(catalog_path.c_str(), "rb"), &gzclose);
        if (!catalog) throw std::runtime_error("cannot open catalogue: " + catalog_path);
        std::cout << "[scan-wl candidates] Streaming complete catalogue " << catalog_path << '\n';
        auto process_line = [&](const std::string& line) {
            ++rows;
            const auto delimiter = line.find_first_of(",\t");
            const std::string token = seq_utils::trim(delimiter == std::string::npos ? line : line.substr(0, delimiter));
            whitelist_relevance::packed_barcode barcode;
            if (!parse_catalog_token(token, *lengths.begin(), barcode)) { ++skipped_rows; return; }
            ++valid_rows;
            if (index.contains(barcode.bits, barcode.length)) retained.insert(barcode);
        };
        std::array<char, 65536> buffer;
        std::string line;
        while (gzgets(catalog.get(), buffer.data(), static_cast<int>(buffer.size()))) {
            line += buffer.data();
            if (!line.empty() && line.back() == '\n') {
                process_line(line);
                line.clear();
            }
        }
        int gzip_error = Z_OK;
        const char* gzip_message = gzerror(catalog.get(), &gzip_error);
        if (gzip_error != Z_OK && gzip_error != Z_STREAM_END)
            throw std::runtime_error("error while reading catalogue " + catalog_path + ": " + gzip_message);
        if (!line.empty()) process_line(line);
        if (gzclose(catalog.release()) != Z_OK)
            throw std::runtime_error("error while closing catalogue: " + catalog_path);
        if (!valid_rows) throw std::runtime_error("catalogue contains no valid barcode identities: " + catalog_path);
        const double catalog_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - catalog_started).count();
        std::vector<std::string> candidates;
        candidates.reserve(retained.size());
        for (const auto& barcode : retained) candidates.push_back(unpack(barcode));
        std::sort(candidates.begin(), candidates.end());
        if (!ensure_output_directory(prefix)) return 1;
        temporary_outputs temporary;
        temporary.paths.push_back(boost::filesystem::unique_path(text_path + ".tmp-%%%%-%%%%-%%%%"));
        temporary.paths.push_back(boost::filesystem::unique_path(log_path + ".tmp-%%%%-%%%%-%%%%"));
        std::ofstream text(temporary.paths[0].string());
        if (!text) throw std::runtime_error("cannot write candidates: " + text_path);
        for (const auto& candidate : candidates) text << candidate << '\n';
        text.close();
        if (!text) throw std::runtime_error("error while writing candidates: " + text_path);
        std::ofstream log(temporary.paths[1].string());
        if (!log) throw std::runtime_error("cannot write candidate log: " + log_path);
        log << std::setprecision(9)
            << "stage=whitelist_candidate_generation\ninput=" << input
            << "\nlayout=" << options.layout << "\nposition_map=" << options.position_map
            << "\ncatalogue=" << catalog_path << "\ncandidates=" << text_path
            << "\nrequested_max_reads=" << max_reads << "\nthreads=" << threads
            << "\nconcat_hmm_requested=" << options.concat_hmm << "\ncandidate_mismatches=" << options.mismatches
            << "\nquery_observations=" << observed.observations_seen
            << "\nraw_orientation_keys=" << observed.raw.size()
            << "\nexact_window_orientation_keys=" << observed.exact.size()
            << "\ncatalogue_rows=" << rows << "\ncatalogue_valid_rows=" << valid_rows
            << "\ncatalogue_skipped_rows=" << skipped_rows << "\nunique_candidate_identities=" << candidates.size()
            << "\nabundance_selection_applied=false\nspatial_selection_applied=false"
            << "\nquery_seconds=" << query_seconds << "\ncatalogue_seconds=" << catalog_seconds
            << "\nwall_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() << '\n';
        log.close();
        if (!log) throw std::runtime_error("error while writing candidate log: " + log_path);
        // Both files are complete before either final path is published.
        // Publish the candidate list last; each rename stays on the same FS.
        boost::filesystem::rename(temporary.paths[1], log_path);
        boost::filesystem::rename(temporary.paths[0], text_path);
        std::cout << "[scan-wl candidates] Wrote " << candidates.size() << " unique candidate identities to " << text_path << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

} // namespace rad_whitelist_candidates
