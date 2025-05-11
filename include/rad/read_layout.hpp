#pragma once
#include "rad_headers.h"

/**
 * @brief Represents a position in the ReadLayout. Each string from PositionStrings is processed in this manner.
 */
struct ParsedPosition {
    std::string ref_id; // Reference element ID this position is based on
    bool is_start; // True if this is a start position, false for stop
    int offset; // Offset from the reference position
    std::string add_flags; // Additional flags (e.g., "terminal_linked", "skipped")

    std::string to_string() const {
        std::stringstream ss;
        ss << ref_id << "|"
           << (is_start ? "start" : "stop")
           << (offset >= 0 ? "+" : "") << offset;
        if (!add_flags.empty()) {
            ss << "|" << add_flags;
        }
        return ss.str();
    }

    static ParsedPosition from_string(const std::string& pos_str, bool verbose) {
        ParsedPosition pos;
        if(pos_str.empty()){
            std::cout << "[from_string] Empty position string\n"; 
            return pos;
        } 

        // First split by pipes to see what we're dealing with
        std::vector<std::string> parts;
        std::string temp;
        std::istringstream pipe_stream(pos_str);
        if(verbose){
            std::cout << "[from_string] Parsing position string: '" << pos_str << "'\n";
        }
        assert(pipe_stream.good());
        while(std::getline(pipe_stream, temp, '|')) {
            if(verbose){
                std::cout << "[from_string] Flags: eof =" << pipe_stream.eof()  << " fail =" << pipe_stream.fail()  << " bad =" << pipe_stream.bad() << "\n";
            }
            parts.push_back(temp);
            if(verbose){
                std::cout << "[from_string] Found part: '" << temp << "'\n";
            }
        }

        // parts[0] = ref_id
        // parts[1] = start/stop+offset or start/stop-offset
        // parts[2] = flags (if exists)
        
        if(parts.size() >= 1) {
            pos.ref_id = parts[0];
            if(verbose) std::cout << "[from_string] Set ref_id: " << pos.ref_id << "\n";      
        }

        if(parts.size() >= 2) {
            // Handle the start/stop and offset part
            std::string pos_and_offset = parts[1];
            pos.is_start = (pos_and_offset.find("start") != std::string::npos);
            std::string is_start = pos.is_start ? "yes" : "no"; 
            if(verbose) std::cout << "[from_string] Is start: " << is_start << "\n";

            // Find the offset
            size_t plus_pos = pos_and_offset.find('+');
            size_t minus_pos = pos_and_offset.find('-');
            
            if(plus_pos != std::string::npos) {
                std::string offset_str = pos_and_offset.substr(plus_pos + 1);
                pos.offset = std::stoi(offset_str);
                if(verbose) std::cout << "[from_string] Found positive offset: " << pos.offset << "\n";
            } else if(minus_pos != std::string::npos) {
                std::string offset_str = pos_and_offset.substr(minus_pos + 1);
                pos.offset = -std::stoi(offset_str);
                if(verbose) std::cout << "[from_string] Found negative offset: " << pos.offset << "\n";
            }
        }

        if(parts.size() >= 3) {
            pos.add_flags = parts[2];
            if(verbose) std::cout << "[from_string] Set flags: " << pos.add_flags << "\n";
        }

        return pos;
    }

};

/**
 * @brief Groups primary and secondary positions for an element. These are the stored breakdowns of the string for easy access.
 */
struct ReferencePositions {
    ParsedPosition primary_start;      //< Primary starting position
    ParsedPosition primary_stop;       //< Primary stopping position
    ParsedPosition secondary_start;    //< Backup starting position
    ParsedPosition secondary_stop;     //< Backup stopping position
};

/**
 * @brief These are alignment positions for start and stop info in non-variable regions.
 */
struct AlignmentPositions {
    std::pair<double, double> start_stats;  // mean, variance
    std::pair<double, double> stop_stats;   // mean, variance
};

// Multi-index container tags
struct id_tag {};               // Tag for accessing by ID
struct length_tag {};           // Tag for accessing by length
struct type_tag {};             // Tag for accessing by type
struct order_tag {};            // Tag for accessing by order
struct direction_tag {};        // Tag for accessing by direction
struct global_class_tag {};     // Tag for accessing by global class
struct dir_order_tag {};        // Tag for accessing by order...and direction!


/**
 * @brief Core element structure for read processing.
 * @brief This structure is used to store information about each element in the read layout.
 * @param class_id Unique identifier for the element class
 * @param seq Raw sequence data
 * @param masked_seq Masked version of the sequence (if static)
 * @param expected_length Expected sequence length if known
 * @param type Element type (e.g., "static", "variable")
 * @param order Position in the layout
 * @param direction Orientation ("forward" or "reverse")
 * @param global_class Global classification
 * @param whitelist_path Path to whitelist file (if applicable)
 * @param misalignment_threshold Threshold for misalignment detection
 * @param aligned_positions Alignment position storage
 * @param misaligned_positions Alignment position storage
 * @param ref_pos Position information
 * 
 */
struct ReadElement {
    std::string class_id;        // Unique identifier for the element class
    std::string seq;             // Raw sequence data
    std::string masked_seq;      // Masked version of the sequence (if static)
    std::optional<int> expected_length;  // Expected sequence length if known
    std::string type;            // Element type (e.g., "static", "variable")
    int order;                   // Position in the layout
    std::string direction;       // Orientation ("forward" or "reverse")
    std::string global_class;    // Global classification
    std::string whitelist_path;     //path to whitelist_file (if applicable)
    std::optional<std::tuple<int, int, int>> misalignment_threshold;  // Threshold for misalignment detection
    std::optional<AlignmentPositions> aligned_positions;  // Alignment position storage
    std::optional<AlignmentPositions> misaligned_positions;  // Alignment position storage
    std::optional<ReferencePositions> ref_pos;  // Position information

    ReadElement(
        const std::string& class_id,
        const std::string& seq,
        const std::string& masked_seq,
        const std::optional<int>& expected_length,
        const std::string& type,
        int order,
        const std::string& direction,
        const std::string& global_class,
        const std::string& whitelist_path = "",
        const std::optional<std::tuple<int, int, int>>& misalignment_threshold = std::nullopt,
        const std::optional<AlignmentPositions>& aligned_positions = std::nullopt,
        const std::optional<AlignmentPositions>& misaligned_positions = std::nullopt,
        const std::optional<ReferencePositions>& ref_pos = std::nullopt
    ) : class_id(std::move(class_id)),
        seq(std::move(seq)),
        masked_seq(std::move(masked_seq)),
        expected_length(expected_length),
        type(std::move(type)),
        order(order),
        direction(std::move(direction)),
        global_class(std::move(global_class)),
        whitelist_path(std::move(whitelist_path)),
        aligned_positions(aligned_positions),
        misaligned_positions(misaligned_positions),
        misalignment_threshold(misalignment_threshold),
        ref_pos(ref_pos) {}
};

/**
 * @brief Multi-index container for read layout elements
 */
typedef boost::multi_index::multi_index_container<
    ReadElement,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<id_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::class_id>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<length_tag>,
            boost::multi_index::member<ReadElement, std::optional<int>, &ReadElement::expected_length>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<type_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::type>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<order_tag>,
            boost::multi_index::member<ReadElement, int, &ReadElement::order>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<direction_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::direction>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<global_class_tag>,
            boost::multi_index::member<ReadElement, std::string, &ReadElement::global_class>
        >,
        boost::multi_index::ordered_non_unique<
            boost::multi_index::tag<dir_order_tag>,
            boost::multi_index::composite_key<
                ReadElement,
                boost::multi_index::member<ReadElement, std::string, &ReadElement::direction>,
                boost::multi_index::member<ReadElement, int, &ReadElement::order>
            >
        >
    >
> Layout_Struct;

/**
 * @brief Class for managing ReadLayout operations
 */
class ReadLayout {
public:
    whitelist wl_map;
    Layout_Struct layout;
    std::unordered_map<std::string, ReferencePositions> position_map;

    //A bunch of access-based methods for the multi-index container
    auto& by_id() { return layout.get<id_tag>(); }
    auto& by_length() { return layout.get<length_tag>(); }
    auto& by_type() { return layout.get<type_tag>(); }
    auto& by_order() { return layout.get<order_tag>(); }
    auto& by_direction() { return layout.get<direction_tag>(); }
    auto& by_global_class() { return layout.get<global_class_tag>(); }
    auto& by_dir_order() { return layout.get<dir_order_tag>(); }

    // Const versions for read-only access
    const auto& by_id() const { return layout.get<id_tag>(); }
    const auto& by_length() const { return layout.get<length_tag>(); }
    const auto& by_type() const { return layout.get<type_tag>(); }
    const auto& by_order() const { return layout.get<order_tag>(); }
    const auto& by_direction() const { return layout.get<direction_tag>(); }
    const auto& by_global_class() const { return layout.get<global_class_tag>(); }
    const auto& by_dir_order() const { return layout.get<dir_order_tag>(); }

    // get read layout mode
    std::string get_rl_mode(std::string input_file) const {
        std::string mode;
        std::ifstream fin(input_file);
        std::string title;
        std::getline(fin, title);           // e.g. "Read Layout:R1"
        std::stringstream ss(title);

        std::string drop;
        std::getline(ss, drop, ':');        // drop “Read Layout”
        if (!std::getline(ss, mode, ':')) { // try to read the part after the first colon
            mode.clear();
        }

        // now trim whitespace and any double-quotes
        const char* trim_chars = " \r\n\t\",";
        auto start = mode.find_first_not_of(trim_chars);
        if (start == std::string::npos) {
            mode.clear();
        } else {
            auto end = mode.find_last_not_of(trim_chars);
            mode = mode.substr(start, end - start + 1);
        }
        return mode;
    }

    // Import data from CSV
    void prep_new_layout(const std::string& input_file, bool verbose) {
        std::unordered_map<std::string, int> class_id_counts;
        std::unordered_map<std::string, int> class_id_total_counts;
        bool build_forward_only, build_reverse_only = false;

        auto log_elem = [&](const ReadElement &e){
            std::cout << "[element] order=" << e.order
                      << " id=" << e.class_id
                      << " type=" << e.type
                      << " dir=" << e.direction
                      << " class=" << e.global_class
                      << " seq=\"" << e.seq << "\""
                      << " exp_len=" << (e.expected_length? std::to_string(*e.expected_length) : "none")
                      << "\n";
        };
       
        std::string mode = get_rl_mode(input_file);
        // set flags
        build_forward_only = (mode == "R1");
        build_reverse_only = (mode == "R2");
        if(verbose){
            std::cout << "[prep_new_layout] Mode: '" << mode << "'\n";
            std::cout << "[prep_new_layout] Build mode: " << (build_forward_only ? "forward"  : build_reverse_only  ? "reverse" : "both")<< "\n";    
        }
        
        // Configure CSV reader
        csv::CSVFormat format;
        format.delimiter(',').header_row(1);
        csv::CSVReader reader(input_file, format);

        // First pass: Collect total counts of class_ids
        std::vector<csv::CSVRow> rows;
        for (auto& row : reader) {
            rows.push_back(row);
            std::string class_id = row["class"].is_null() || row["class"].get<std::string>().empty()
                                    ? row["id"].get<std::string>() : row["class"].get<std::string>();
            class_id_total_counts[class_id]++;
        }

        // Initialize order counter and insert seq_start
        int order_counter = 1;
        ReadElement seq_start("seq_start", "", "", 0, "static", order_counter++, "forward", "start");
        //if(!build_reverse_only) 
        layout.insert(seq_start);

        // Process forward elements
        for (auto& row : rows) {
            std::string id = row["id"].is_null() ? "" : row["id"].get<std::string>();
            std::string seq = row["seq"].is_null() ? "" : row["seq"].get<std::string>();
            std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            std::string type = row["type"].is_null() ? "variable" : row["type"].get<std::string>();
            std::string global_class = row["class"].is_null() ? "" : row["class"].get<std::string>();
            std::string direction = "forward";
            std::string class_id = row["class"].is_null() ? id : row["class"].get<std::string>();
            std::string whitelist_path = row["whitelist"].is_null() ? "" : row["whitelist"].get<std::string>();

            // Fill empty class, class_id, or expected length
            std::optional<int> expected_length;
            if (!row["expected_length"].is_null()) {
                expected_length = row["expected_length"].get<int>();
            }
            if (global_class.empty()) {
                global_class = class_id;
            }
            if (class_id.empty()) {
                class_id = global_class;
            }

            //populate poly-tail cases
            if (class_id == "poly_t") {
                expected_length = expected_length.value_or(12);
                seq = "T{" + std::to_string(*expected_length) + ",}+";
                global_class = "poly_tail";
            } else if (class_id == "poly_a") {
                expected_length = expected_length.value_or(12);
                seq = "A{" + std::to_string(*expected_length) + ",}+";
                global_class = "poly_tail";
            }

            if (class_id_total_counts[class_id] > 1) {
                int &cnt = class_id_counts[class_id];
                ++cnt;
                std::string original = class_id;
                class_id += "_" + std::to_string(cnt);
                if (verbose) {
                    std::cout << "[prep_new_layout] Duplicate class " << original << "' found (total = " 
                              << class_id_total_counts[original] << "), renaming to " << class_id << "\n";
                }
            }

            if (!seq.empty() && !expected_length.has_value()) {
                expected_length = static_cast<int>(seq.length());
            }
            std::string masked_seq = (type == "static" && global_class != "poly_tail") ? "MASKED" : "";

            ReadElement elem(class_id, seq, masked_seq, expected_length, type, order_counter++, direction, global_class, whitelist_path);
            //if(!build_reverse_only) 
            layout.insert(elem);
            if (verbose) {
                std::cout << "[read_layout] Inserted element: ";
                log_elem(elem);
            }
        }
        // Insert seq_stop
        ReadElement seq_stop("seq_stop", "", "", 0, "static", order_counter++, "forward", "stop");
        //if(!build_reverse_only) 
        layout.insert(seq_stop);
        if(build_forward_only){
            std::cout << "[prep_new_layout] Forward-only (R1) layout generated.\n";
            return;
        }
        // Generate reverse complement entries
        std::vector<ReadElement> reverse_elements;
        // Determine the starting point for reverse complement orders
        int reverse_order_start = order_counter;
        const auto& ordered_index = by_order();

        for (const auto& elem : ordered_index) {
            ReadElement reverse_elem = elem;
            reverse_elem.direction = "reverse";
            if (elem.class_id == "poly_t") {
                reverse_elem.global_class = "poly_tail";
                reverse_elem.seq = "A{" + std::to_string(elem.expected_length.value_or(12)) + ",}+";
                reverse_elem.class_id = "poly_a";
            } else if (elem.class_id == "poly_a") {
                reverse_elem.global_class = "poly_tail";
                reverse_elem.seq = "T{" + std::to_string(elem.expected_length.value_or(12)) + ",}+";
                reverse_elem.class_id = "poly_t";
            } else {
                reverse_elem.seq = seq_utils::revcomp(elem.seq);
                reverse_elem.class_id = "rc_" + elem.class_id;
            }
            reverse_elements.push_back(reverse_elem);
        }
        // Organize reverse complement rows
        auto start_it = std::stable_partition(reverse_elements.begin(), reverse_elements.end(),
            [](const ReadElement& elem) { return elem.global_class == "start"; });
        auto stop_it = std::stable_partition(start_it, reverse_elements.end(),
            [](const ReadElement& elem) { return elem.global_class != "stop"; });

        // Reverse the middle portion
        std::reverse(start_it, stop_it);
        // Assign order to reverse complement elements and add them to layout
        for (auto& reverse_elem : reverse_elements) {
            reverse_elem.order = reverse_order_start++;
            layout.insert(reverse_elem);
            if (verbose) {
                std::cout << "[prep_new_layout] Inserted element: ";
                log_elem(reverse_elem);
            }
        }
        if (build_reverse_only) {
            // Get the index keyed on direction:
            auto& dir_index = layout.get<direction_tag>();
            // Find the sub-range where direction == "forward"
            auto range = dir_index.equal_range("forward");
            // Erase them
            //dir_index.erase(range.first, range.second);
            std::cout << "[prep_new_layout] Removed all forward elements; reverse-only layout ready. New layout :\n";
            for(auto& elem : layout) {
                log_elem(elem);
            }
            return;
        }
    }

    // Export layout to CSV
    void write_to_csv(const std::string& base_path = "", const std::string& which_file = "both") const {
        std::string path_prefix = base_path.empty() ?
            std::string(std::getenv("HOME")) + "/Desktop/" : base_path;
        
        if(which_file == "layout" || which_file == "both") {
            // Basic layout information
            std::ofstream layout_file(path_prefix + "_layout.csv");
            layout_file << "id,seq,masked_seq,expected_length,type,class,direction,class_id,whitelist,order\n";
            for (const auto& element : by_order()) {
                layout_file << "\"" << element.class_id << "\","
                    << "\"" << element.seq << "\","
                    << (element.masked_seq.empty() ? "" : element.masked_seq) << ","
                    << (element.expected_length.has_value() ? std::to_string(*element.expected_length) : "") << ","
                    << element.type << ","
                    << element.global_class << ","
                    << element.direction << ","
                    << element.class_id << ","
                    << (element.whitelist_path.empty() ? "" : element.whitelist_path) << ","
                    << element.order << "\n";
            }
        }

        if(which_file == "pos_map" || which_file == "both") {
            // Alignment and threshold information
            std::ofstream align_file(path_prefix + "_position_map.csv");
            align_file << "id,primary_start,primary_stop,secondary_start,secondary_stop,"
                << "misalign_lower,misalign_mean,misalign_upper,"
                << "align_start,var_start,align_stop,var_stop,"
                << "misalign_start,mvar_start,misalign_stop,mvar_stop\n";
            
            for (const auto& element : by_order()) {
                align_file << element.class_id;
                
                // Position information
                auto pos_it = position_map.find(element.class_id);
                if (pos_it != position_map.end()) {
                    align_file << "," << pos_it->second.primary_start.to_string()
                            << "," << pos_it->second.primary_stop.to_string()
                            << "," << pos_it->second.secondary_start.to_string()
                            << "," << pos_it->second.secondary_stop.to_string();
                } else {
                    // Four commas for four empty fields
                    align_file << ",,,,";  
                }

                // Misalignment thresholds and positions
                if (element.misalignment_threshold) {
                    const auto& [lower, mean, upper] = *element.misalignment_threshold;
                    align_file << "," << lower << "," << mean << "," << upper;
                    
                    if (element.aligned_positions) {
                        const auto& [start_mean, start_var] = element.aligned_positions->start_stats;
                        const auto& [stop_mean, stop_var] = element.aligned_positions->stop_stats;
                        align_file << "," << start_mean << "," << start_var 
                                << "," << stop_mean << "," << stop_var;
                    } else {
                        // Four commas for missing aligned positions
                        align_file << ",,,,";  
                    }
                    
                    if (element.misaligned_positions) {
                        const auto& [mis_start_mean, mis_start_var] = element.misaligned_positions->start_stats;
                        const auto& [mis_stop_mean, mis_stop_var] = element.misaligned_positions->stop_stats;
                        align_file << "," << mis_start_mean << "," << mis_start_var 
                                << "," << mis_stop_mean << "," << mis_stop_var;
                    } else {
                        // Four commas for missing misaligned positions
                        align_file << ",,,,";
                    }
                } else {
                    // If no misalignment threshold, need commas for all missing fields
                    align_file << ",,,";  // Three commas for threshold fields
                    align_file << ",,,,";  // Four commas for aligned positions
                    align_file << ",,,,";  // Four commas for misaligned positions
                }
                align_file << "\n";
            }
        }
    }

    // Import layout from CSVs
    void import_from_csv(const std::string& base_path = "", bool verbose = false) {
        // Read main layout file
        csv::CSVFormat format;
        format.delimiter(',').quote('"').header_row(0);
        format.variable_columns(csv::VariableColumnPolicy::THROW);
        csv::CSVReader layout_reader(base_path + "_layout.csv", format);

        //int order_counter = 1;
        for (auto& row : layout_reader) {
            ReadElement elem(
                row["id"].get<std::string>(),
                row["seq"].get<std::string>(),
                row["masked_seq"].get<std::string>(),
                row["expected_length"].is_null() ? std::nullopt : 
                    std::optional<int>(std::stoi(row["expected_length"].get<std::string>())),
                row["type"].get<std::string>(),
                std::stoi(row["order"].get<std::string>()),
                row["direction"].get<std::string>(),
                row["class"].get<std::string>(),
                row["whitelist"].get<std::string>()
            );
            layout.insert(elem);
        }

        // Read position and alignment info
        csv::CSVFormat new_format;
        new_format.delimiter(',').quote('"').header_row(0);
        new_format.variable_columns(csv::VariableColumnPolicy::KEEP);
        csv::CSVReader pos_reader(base_path + "_position_map.csv", new_format);
        for (auto& row : pos_reader) {
            std::string id = row["id"].get<std::string>();
            if(verbose){
                std::cout << "[import_from_both] Processing ID: " << id << "\n";
            }
            auto element_it = by_id().find(id);
            if (element_it == by_id().end()){
                if(verbose){
                    std::cout << "[import_from_both] Skipping ID not found in layout: " << id << "\n";
                }
                continue;
            }
        // Parse position info
            if (!row["primary_start"].get<std::string>().empty()) {
            ReferencePositions pos;
            if(verbose){
                std::cout << "[import_from_both] This did work!\n";
                std::cout << "[import_from_both] Row data for " << id << ":\n";
                std::cout << "[import_from_both] primary_start value: '" << row["primary_start"].get<std::string>() << "'\n";
            }
            try {
                    pos.primary_start = ParsedPosition::from_string(row["primary_start"].get<std::string>(), verbose);
                    pos.primary_stop = ParsedPosition::from_string(row["primary_stop"].get<std::string>(), verbose);
                    pos.secondary_start = ParsedPosition::from_string(row["secondary_start"].get<std::string>(), verbose);
                    pos.secondary_stop = ParsedPosition::from_string(row["secondary_stop"].get<std::string>(), verbose);
            
                    position_map[id] = pos;
                    if(verbose){
                        std::cout << "[import_from_both] Processing element " << id << " with positions:\n"
                        << "Primary: " << pos.primary_start.to_string() << " -> " << pos.primary_stop.to_string() << "\n"
                        << "Secondary: " << pos.secondary_start.to_string() << " -> " << pos.secondary_stop.to_string() << "\n";
                    }
                    
                    by_id().modify(element_it, [pos](ReadElement& elem) {
                        elem.ref_pos = pos;
                    });
                } catch (const std::exception& e) {
                    std::cerr << "Error processing positions for " << id << ": " << e.what() << "\n";
                }
            } else {
                if(verbose){
                    std::cout << "[import_from_both] No primary start value!\n";
                    std::cout << "[import_from_both] Row data for " << id << ":\n";    
                }
        }
        // Parse misalignment and alignment stats
            if (!row["misalign_lower"].is_null()) {
                std::tuple<int, int, int> threshold{
                    std::stoi(row["misalign_lower"].get<std::string>()),
                    std::stoi(row["misalign_mean"].get<std::string>()),
                    std::stoi(row["misalign_upper"].get<std::string>())
                };

                AlignmentPositions align_pos{
                    {std::stod(row["align_start"].get<std::string>()), 
                        std::stod(row["var_start"].get<std::string>())},
                    {std::stod(row["align_stop"].get<std::string>()), 
                        std::stod(row["var_stop"].get<std::string>())}
                };

                AlignmentPositions misalign_pos{
                    {std::stod(row["misalign_start"].get<std::string>()), 
                        std::stod(row["mvar_start"].get<std::string>())},
                    {std::stod(row["misalign_stop"].get<std::string>()), 
                        std::stod(row["mvar_stop"].get<std::string>())}
                };

                by_id().modify(element_it, [&](ReadElement& elem) {
                    elem.misalignment_threshold = threshold;
                    elem.aligned_positions = align_pos;
                    elem.misaligned_positions = misalign_pos;
                });
            }
        }
    }
    
    // import a read layout .csv
    void import_read_layout(const std::string& layout_csv, bool verbose){
        using namespace csv;
        CSVFormat fmt;
        fmt.delimiter(',').quote('"').header_row(0)
        .variable_columns(VariableColumnPolicy::THROW);
        CSVReader reader(layout_csv, fmt);
        for (auto &row : reader) {
            ReadElement elem(
                row["id"].get<std::string>(),
                row["seq"].get<std::string>(),
                row["masked_seq"].get<std::string>(),
                row["expected_length"].is_null() ? std::nullopt : std::optional<int>(std::stoi(row["expected_length"].get<std::string>())),
                row["type"].get<std::string>(),
                std::stoi(row["order"].get<std::string>()),
                row["direction"].get<std::string>(),
                row["class"].get<std::string>(),
                row["whitelist"].get<std::string>()
            );
            layout.insert(elem);
        }
        if (verbose) {
            std::vector<ReadElement> all;
            all.reserve(layout.size());
            for (auto const &e : layout) {
                all.push_back(e);
            }
            std::vector<ReadElement> fwd, rev;
            for (auto const &e : all) {
                if (e.direction == "forward")  fwd.push_back(e);
                if (e.direction == "reverse")  rev.push_back(e);
            }
            // sort by order field
            auto cmp_order = [](auto const &a, auto const &b){
                return a.order < b.order;
            };
            std::sort(fwd.begin(), fwd.end(), cmp_order);
            std::sort(rev.begin(), rev.end(), cmp_order);

            // helper to print a list
            auto print_list = [&](auto const &vec, const char *name){
                std::cout << "[read_layout]:" << name;
                for (auto const &e : vec) {
                    std::string seq = e.seq.empty() ? "" : e.seq;
                    std::cout << "["
                            << e.class_id << ":"
                            << seq
                            << "]";
                }
                std::cout << "\n";
            };
            print_list(fwd, "F:");
            print_list(rev, "R:");
        }
    }
    
    // import a position map .csv
    void import_position_map(const std::string& map_csv,bool verbose){
        using namespace csv;
        // set up CSV reader for the position map
        CSVFormat fmt;
        fmt.delimiter(',').quote('"').header_row(0)
            .variable_columns(VariableColumnPolicy::KEEP);
        CSVReader reader(map_csv, fmt);
        // iterate rows
        for (auto &row : reader) {
            std::string id = row["id"].get<std::string>();
            // find the corresponding element from the layout
            auto element_it = by_id().find(id);
            if (element_it == by_id().end()) {
                if (verbose) {
                    std::cout << "[pos_map] skipping unknown id: " << id << "\n";
                }
                continue;
            }

            // only proceed if we have a non-empty primary_start
            const std::string &ps = row["primary_start"].get<std::string>();
            if (!ps.empty()) {
                // parse reference positions
                ReferencePositions rp;
                rp.primary_start   = ParsedPosition::from_string(ps, verbose);
                rp.primary_stop    = ParsedPosition::from_string(
                                            row["primary_stop"].get<std::string>(), verbose);
                rp.secondary_start = ParsedPosition::from_string(
                                            row["secondary_start"].get<std::string>(), verbose);
                rp.secondary_stop  = ParsedPosition::from_string(
                                            row["secondary_stop"].get<std::string>(), verbose);

                // store in the map and modify the element
                position_map[id] = rp;
                by_id().modify(element_it, [rp](ReadElement &e){
                    e.ref_pos = rp;
                });
            if (verbose) {
                std::cout << "[pos_map] " << id
                            << " -> primary("
                            << rp.primary_start.to_string() << "-"
                            << rp.primary_stop.to_string()   << ")\n";
            }
        }

        // now, if misalignment stats are present, parse and attach
        if (!row["misalign_lower"].is_null()) {
                auto lower = std::stoi(row["misalign_lower"].get<std::string>());
                auto mean  = std::stoi(row["misalign_mean" ].get<std::string>());
                auto upper = std::stoi(row["misalign_upper"].get<std::string>());
                std::tuple<int, int, int> threshold{lower, mean, upper};
                if(verbose){
                    std::cout << "[pos_map] " << id
                                << " -> misalignment threshold: "
                                << lower << "|" << mean << "|" << upper << "\n";
                }

                // alignment positions
                AlignmentPositions align_pos{
                    {
                        std::stod(row["align_start"].get<std::string>()),
                        std::stod(row["var_start"].get<std::string>())
                    },
                    {
                        std::stod(row["align_stop"].get<std::string>()),
                        std::stod(row["var_stop"].get<std::string>())
                    }
                };
                // misalignment positions
                AlignmentPositions misalign_pos{
                    {
                        std::stod(row["misalign_start"].get<std::string>()),
                        std::stod(row["mvar_start"].get<std::string>())
                    },
                    {
                        std::stod(row["misalign_stop" ].get<std::string>()),
                        std::stod(row["mvar_stop"].get<std::string>())
                    }
                };

                // attach them both in one modify call
                by_id().modify(element_it, [&](ReadElement &e){
                    e.misalignment_threshold = threshold;
                    e.aligned_positions      = align_pos;
                    e.misaligned_positions   = misalign_pos;
                });

                if (verbose) {
                    std::cout << "[pos_map] " << id
                                << " -> misalignment threshold: "
                                << lower << "|" << mean << "|" << upper << "\n";
                }
            }
        }
    }

    // Generate position map
    void generate_position_mapping(){
        auto& dir_ordered = by_dir_order();
        auto current_direction = dir_ordered.begin();
        while (current_direction != dir_ordered.end()) {
            auto direction = current_direction->direction;
            auto end_of_direction = dir_ordered.upper_bound(
                boost::make_tuple(direction, std::numeric_limits<int>::max())
            );
            process_reads(current_direction, end_of_direction);
                auto& id_index = by_id();
                for (const auto& pos_pair : position_map) {
                    auto elem_it = id_index.find(pos_pair.first);
                    if (elem_it != id_index.end()) {
                        id_index.modify(elem_it, [positions=pos_pair.second](ReadElement& elem) {
                            elem.ref_pos = positions;
                        });
                    }
                }
            current_direction = end_of_direction;
        }
    }

    //load whitelists from whitelist path in read layout
    /*
    void load_wl(std::optional<std::string> wl_path, std::optional<int> shift, std::optional<int> mut, bool verbose, int nthreads) {
        using namespace std::chrono;
        auto t0 = high_resolution_clock::now();
    
        // collect all static seqs for filter_bcs
        std::vector<std::string> static_seqs;
        for (auto const &elem : layout) {
            if (elem.type == "static" && !elem.seq.empty())
                static_seqs.push_back(elem.seq);
        }
    
        // map each unique path to one wl_entry
        std::unordered_map<std::string, whitelist::wl_entry> path2entry;
    
        for (auto const &elem : layout) {
            if (elem.global_class != "barcode" || elem.type != "variable")
                continue;
    
            std::string key = elem.class_id;
            // normalize the key (strip rc_)
            if (seq_utils::is_rc(key))
                key = seq_utils::remove_rc(key);
    
            // resolve kit or path
            std::string spec;
            if(wl_path.has_value()){
                spec = wl_path.value();
            } else {
                spec = elem.whitelist_path;
            }
            std::string path = whitelist_utils::kit_to_path(spec);
    
            std::cout << "[load_wl] "
                      << "[" << elem.class_id << "]" 
                      << "[" << spec << "] @ " << path << "\n";
    
            // if we've never loaded this file before, do so now
            auto pit = path2entry.find(path);
            if (pit == path2entry.end()) {
                if(verbose) std::cout << "[load_wl] Loading from " << path << " ...\n";
                auto t1 = high_resolution_clock::now();
                // get expected length
                size_t default_length = 16;
                if (auto lit = layout.find(elem.class_id);
                    lit != layout.end() && lit->expected_length)
                {
                    default_length = *lit->expected_length;
                }
    
                // import & possibly generate mismatches
                auto entry = wl_map.import_whitelist(spec, verbose, default_length);
                //this is to generate mismatches for subsets that have been passed in--otherwise we're SOL
                if (entry.true_bcs.size() <= 15000){
                    int shift_amt = shift.value_or(2);
                    int mut_amt   = mut.value_or(2);
                    std::cout << "[load_wl] " << 
                    "[" << elem.class_id << ":" << key << "]"
                              << ", true barcodes = "    << entry.true_bcs.size()
                              << "\n";        
                    entry.generate_mismatch_barcodes(shift_amt, mut_amt, verbose, nthreads);
                }
    
                // build filter_bcs
                entry.filter_bcs.clear();
                size_t def_len = 0;
                if (auto lit = layout.find(elem.class_id);
                    lit != layout.end() && lit->expected_length)
                {
                    def_len = *lit->expected_length;
                }
                for (auto const &s : static_seqs) {
                    for (auto const &kmer : seq_utils::circ_kmerize(s, def_len)) {
                        int64_seq bits(kmer);
                        entry.filter_bcs.insert_bc_entry(bits);
                    }
                }
                auto t2 = high_resolution_clock::now();
                double ms = duration<double, std::milli>(t2 - t1).count();
                std::cout << "[load_wl] loaded in " << (ms/1000.0) << " s\n";
    
                pit = path2entry.emplace(path, std::move(entry)).first;
            }
    
            // map normalized class_id -> the already‐loaded wl_entry
            auto [it, inserted] = wl_map.lists.try_emplace(key, pit->second);
            if(inserted){
                std::cout << "[load_wl] " << "New entry for " << key << "\n";
            } else {
                std::cout << "[load_wl] " << "Reusing entry for " << key << "\n";
            }

            // report
            auto &E = wl_map.lists[key];
            std::cout << "[load_wl] " << 
            "[" << elem.class_id << ":" << key << "]"
                      << ": global = "  << E.global_bcs.size()
                      << ", true = "    << E.true_bcs.size()
                      << ", filter = "  << E.filter_bcs.size()
                      << "\n";
        }

        for (auto & [key, entry] : wl_map.lists) {
            // 1) collect unique barcodes in true_bcs
            std::unordered_set<int64_seq> to_remove;
            to_remove.reserve(entry.true_bcs.size());
            for (auto const & [bc, be] : entry.true_bcs) {
                to_remove.insert(bc);
            }
        
            // 2) drop them all from global_bcs
            for (auto const & bc : to_remove) {
                entry.global_bcs.remove_bc_entry(bc);
            }

            auto &E = wl_map.lists[key];
            std::cout << "[load_wl] " << 
            "[" << key << "]"
                      << ": global = "  << E.global_bcs.size()
                      << ", true = "    << E.true_bcs.size()
                      << ", filter = "  << E.filter_bcs.size()
                      << "\n";

       }
    
        auto t3 = high_resolution_clock::now();
        double total_s = duration<double>(t3 - t0).count();
        std::cout << "[load_wl] finished loading all whitelists in " << total_s << " s\n";
    }
    */

    void load_wl(std::optional<std::string> wl_path,std::optional<int> shift, std::optional<int> mut, bool verbose, int nthreads) {
        using namespace std::chrono;
        auto t0 = high_resolution_clock::now();

        // Clear any old state
        wl_map.lists.clear();
        wl_map.maps.clear();

        // 1) collect static seqs for filter_bcs
        std::vector<std::string> static_seqs;
        for (auto const &elem : layout) {
            if (elem.type=="static" && !elem.seq.empty())
                static_seqs.push_back(elem.seq);
        }

        // 2) import & alias
        for (auto const &elem : layout) {
            if (elem.global_class!="barcode" || elem.type!="variable")
                continue;

            // normalize class_id → key
            std::string key = elem.class_id;
            if (seq_utils::is_rc(key))
                key = seq_utils::remove_rc(key);

            // resolve the on-disk path
            std::string spec = wl_path ? *wl_path : elem.whitelist_path;
            std::string path = whitelist_utils::kit_to_path(spec);

            std::cout << "[load_wl] ["<<elem.class_id<<"] ["<<spec<<"] @ "<<path<<"\n";

            // 2a) import into the pool if not already
            auto pit = wl_map.lists.find(path);
            if (pit == wl_map.lists.end()) {
                if (verbose) std::cout << "[load_wl] Loading from " << path << " ...\n";
                auto t1 = high_resolution_clock::now();

                // get expected_length
                size_t default_length = 16;
                if (auto lit = layout.find(elem.class_id);
                    lit != layout.end() && lit->expected_length)
                {
                    default_length = *lit->expected_length;
                }

                // import & possibly generate mismatches
                auto entry = wl_map.import_whitelist(spec, verbose, default_length);
                if (entry.true_bcs.size() <= 15000) {
                    entry.generate_mismatch_barcodes(shift.value_or(2), mut.value_or(2), verbose, nthreads);
                }

                // build filter_bcs from static_seqs
                entry.filter_bcs.clear();
                size_t def_len = default_length;
                if (auto lit = layout.find(elem.class_id);
                    lit != layout.end() && lit->expected_length)
                {
                    def_len = *lit->expected_length;
                }
                for (auto const &s : static_seqs) {
                    for (auto const &kmer : seq_utils::circ_kmerize(s, def_len)) {
                        int64_seq bits(kmer);
                        entry.filter_bcs.insert_bc_entry(bits);
                    }
                }

                auto t2 = high_resolution_clock::now();
                double ms = duration<double, std::milli>(t2 - t1).count();
                std::cout << "[load_wl] loaded in " << (ms/1000.0) << " s\n";

                pit = wl_map.lists.emplace(path, std::move(entry)).first;
            }

            // 2b) alias into maps by reference
            auto &master = pit->second;
            auto [mit, inserted] =
            wl_map.maps.try_emplace(key, std::ref(master));

            if (inserted) {
                std::cout << "[load_wl] New entry for key ='"<< key << "'\n";
            } else if (verbose) {
                std::cout << "[load_wl] Reusing entry for key ='"<< key << "'\n";
            }

            {
                auto &E = mit->second.get();
                std::cout << "[load_wl]"
                        << " ["<<elem.class_id<<":"<<key<<"]"
                        << ": global ="<<E.global_bcs.size()
                        << ", true ="  <<E.true_bcs.size()
                        << ", filter ="<<E.filter_bcs.size()
                        << "\n";
            }
        }

        // 3) final prune: for each alias, remove true_bcs from global_bcs
        for (auto & [key, entry_ref] : wl_map.maps) {
            auto &E = entry_ref.get();
            std::unordered_set<int64_seq> to_remove;
            to_remove.reserve(E.true_bcs.size());
            for (auto const & [bc,be] : E.true_bcs)
                to_remove.insert(bc);
            for (auto const &bc : to_remove)
                E.global_bcs.remove_bc_entry(bc);

            // final report
            std::cout << "[load_wl]"
                    << " ["<<key<<"]"
                    << ": global ="<<E.global_bcs.size()
                    << ", true ="  <<E.true_bcs.size()
                    << ", filter ="<<E.filter_bcs.size()
                    << "\n";
        }

        auto t3 = high_resolution_clock::now();
        double total_s = duration<double>(t3 - t0).count();
        std::cout << "[load_wl] finished loading all whitelists in "
                << total_s << " s\n";
    }

    // save the whitelist to a CSV 
    void save_wl(std::ostream &out, bool verbose, bool full = false) const {
        //Loop over every whitelist in wl_map.lists
        for (auto const& [class_id, wrap] : wl_map.maps) {
            auto &entry = wrap.get();
            // return the true_bcs
            if(verbose) std::cout << "[save_wl]"
            << " ["<<class_id<<"]"
            << ": global = "<<entry.global_bcs.size()
            << ", true = "  <<entry.true_bcs.size()
            << ", filter = "<<entry.filter_bcs.size()
            << "\n";

            {
            if(!full){
                entry.true_bcs.write_wl_summary(out, class_id, &entry.true_ref);
                continue;
            }
            out << "class_id,source,barcode,"
            "total_count,corrected_count,"
            "forw_count,rev_count,"
            "forw_concat_count,rev_concat_count,"
            "filtered_count\n";
            auto rows = entry.true_bcs.summarize_counts(&entry.true_ref);
                for (auto const &tpl : rows) {
                    out << class_id << ",true,"
                      << std::get<0>(tpl) << ','
                      << std::get<1>(tpl) << ','
                      << std::get<2>(tpl) << ','
                      << std::get<3>(tpl) << ','
                      << std::get<4>(tpl) << ','
                      << std::get<5>(tpl) << ','
                      << std::get<6>(tpl) << ','
                      << std::get<7>(tpl) << '\n';
                }
            }
            // return the global_bcs *only* for keys seen (total_count>0)
            {
                // build the set of “seen” keys
                std::unordered_set<int64_seq> seen;
                seen.reserve(entry.true_bcs.size());
                for (auto const &kv : entry.true_bcs) {
                    auto cnt = kv.second.count.load(barcode_counts::total);
                    if (cnt > 0) seen.insert(kv.first);
                }
    
                auto rows = entry.global_bcs.summarize_counts(&seen);
                for (auto const &tpl : rows) {
                    out << class_id << ",global,"
                      << std::get<0>(tpl) << ','
                      << std::get<1>(tpl) << ','
                      << std::get<2>(tpl) << ','
                      << std::get<3>(tpl) << ','
                      << std::get<4>(tpl) << ','
                      << std::get<5>(tpl) << ','
                      << std::get<6>(tpl) << ','
                      << std::get<7>(tpl) << '\n';
                }
            }
        }
    }

    // save wl but instead of streaming to output, take a path
    void save_wl(const std::string &path, bool verbose, bool full = false) const {
        std::ofstream out(path);
        if (!out.is_open()) {
            std::cerr << "[error] Error opening file for writing: " << path << "\n";
            return;
        }
        save_wl(out, full);
        out.close();
    }

    // print full layout
    void display_read_layout() {
        std::cout << "[read_layout] COMPLETE_READ_LAYOUT:\n";
        for (const auto& elem : layout) {
            // Base fields
            std::cout
                << "ID: "             << elem.class_id
                << ", Seq: "          << elem.seq
                << ", Type: "         << elem.type
                << ", Order: "        << elem.order
                << ", Direction: "    << elem.direction
                << ", Global Class: " << elem.global_class;
    
            // Optional: misalignment_threshold (tuple<int,int,int>)
            if (elem.misalignment_threshold) {
                auto [min_v, max_v, thr] = *elem.misalignment_threshold;
                std::cout << ", Misalign Thresh: (min=" << min_v
                          << ", max=" << max_v
                          << ", thr=" << thr << ")";
            }
    
            // Optional: aligned_positions
            if (elem.aligned_positions) {
                const auto& ap = *elem.aligned_positions;
                std::cout << ", Aligned Pos: (start: mean="   << ap.start_stats.first
                          << ", var="     << ap.start_stats.second
                          << "; stop: mean=" << ap.stop_stats.first
                          << ", var="     << ap.stop_stats.second << ")";
            }
    
            // Optional: misaligned_positions
            if (elem.misaligned_positions) {
                const auto& mp = *elem.misaligned_positions;
                std::cout << ", Misaligned Pos: (start: mean="   << mp.start_stats.first
                          << ", var="     << mp.start_stats.second
                          << "; stop: mean=" << mp.stop_stats.first
                          << ", var="     << mp.stop_stats.second << ")";
            }

            if(elem.ref_pos) {
                const auto& rp = *elem.ref_pos;
                std::cout << ", Ref Pos: (primary_start: " << rp.primary_start.to_string()
                          << ", primary_stop: " << rp.primary_stop.to_string()
                          << ", secondary_start: " << rp.secondary_start.to_string()
                          << ", secondary_stop: " << rp.secondary_stop.to_string() << ")";
            }
    
            std::cout << "\n";
        }
    }

    // Other utility methods for managing layout
    size_t size() const { return layout.size(); }
    bool empty() const { return layout.empty(); }
    void clear() { layout.clear(); }

private:

    void process_left_side(Layout_Struct::index<dir_order_tag>::type::iterator it) {
        auto& ordered = by_order();
        auto next = ordered.find(it->order + 1);
            auto next_next = ordered.find(it->order + 2);

        auto prev = ordered.find(it->order - 1);
        if (prev == ordered.end()) return;

        ReferencePositions positions;
        int offset = it->expected_length.value_or(1);

        // Always set primary positions relative to previous element
        positions.primary_start = {prev->class_id, false, 1, ""};
        positions.primary_stop = {prev->class_id, false, offset, ""};

        // Handle secondary positions based on types
        if (prev->type == "static") {
            // Check for static anchor
            if (next != ordered.end() && next->type == "static" && 
                next->expected_length.value_or(0) >= 13 && 
                next->global_class != "poly_tail") {
                    positions.secondary_start = {next->class_id, true, -offset, ""};
                    positions.secondary_stop = {next->class_id, true, -1, ""};
            } else {
                if(next->type == "variable" && next_next->global_class != "read" && next_next->type == "static"){
                    int length_offset = next->expected_length.value_or(1);
                    positions.secondary_start = {next_next->class_id, true, -length_offset - offset, ""};
                   positions.secondary_stop = {next_next->class_id, true, -length_offset - 1, ""};
                } else {
                    positions.secondary_start = {prev->class_id, false, 1, "left_terminal_linked"};
                    positions.secondary_stop = {prev->class_id, false, offset, "left_terminal_linked"};
                }
            }
        } else if (prev->type == "variable" || prev->global_class == "poly_tail") {
            std::string flag_type = prev->type == "variable" ? "var_chained" : "poly_chained";
            positions.primary_start.add_flags = flag_type + "_start";
            positions.primary_stop.add_flags = flag_type + "_stop";
            // Look for previous static anchor first
            auto prev_static = ordered.find(it->order - 2);  // Look two positions back
            if (prev_static != ordered.end() && prev_static->type == "static" && 
                prev_static->global_class != "poly_tail") {
                int prev_length = prev->expected_length.value_or(1);
                positions.secondary_start = {prev_static->class_id, false, 1 + prev_length, ""};
                positions.secondary_stop = {prev_static->class_id, false, offset + prev_length, ""};
            } else if (next != ordered.end() && next->type == "static") {
                // Then try next static
                positions.secondary_start = {next->class_id, true, -offset, ""};
                positions.secondary_stop = {next->class_id, true, -1, ""};
            } else if (auto prev_pos = position_map.find(prev->class_id); prev_pos != position_map.end()) {
                // Finally fall back to chaining
                int prev_length = prev->expected_length.value_or(1);
                positions.secondary_start = {prev_pos->second.primary_start.ref_id, false, 1 + prev_length, "chained_start"};
                positions.secondary_stop = {prev_pos->second.primary_stop.ref_id, false, 1 + prev_length + offset, "chained_stop"};
            }
        }

        position_map[it->class_id] = positions;
    }

    void process_right_side(Layout_Struct::index<dir_order_tag>::type::iterator it) {
        auto& ordered = by_order();
        auto next = ordered.find(it->order + 1);

        auto prev = ordered.find(it->order - 1);
            auto prev_prev = ordered.find(it->order - 2);

        if (next == ordered.end()) return;

        ReferencePositions positions;
        int offset = it->expected_length.value_or(1);

        // Always set primary positions relative to next element
        positions.primary_start = {next->class_id, true, -offset, ""};
        positions.primary_stop = {next->class_id, true, -1, ""};

        // Handle secondary positions based on types
        if (next->type == "static") {
            // Check for static anchor
            if (prev != ordered.end() && prev->type == "static" && 
                prev->expected_length.value_or(0) >= 13 && 
                prev->global_class != "poly_tail") {
                    positions.secondary_start = {prev->class_id, false, 1, ""};
                    positions.secondary_stop = {prev->class_id, false, offset, ""};
            } else {

                if(prev->type == "variable" && prev_prev->global_class != "read" && prev_prev->type == "static"){
                    int length_offset = prev->expected_length.value_or(1);
                    positions.secondary_start = {prev_prev->class_id, false, length_offset + 1, ""};
                    positions.secondary_stop = {prev_prev->class_id, false, length_offset + offset, ""};
                } else {
                    positions.secondary_start = {next->class_id, true, -offset, "right_terminal_linked"};
                    positions.secondary_stop = {next->class_id, true, -1, "right_terminal_linked"};
                }
            }
        } else if (next->type == "variable" || next->global_class == "poly_tail") {
            std::string flag_type = prev->type == "variable" ? "var_chained" : "poly_chained";
            positions.primary_start.add_flags = flag_type + "_start";
            positions.primary_stop.add_flags = flag_type + "_stop";

            // Look for next static anchor first
            auto next_static = ordered.find(it->order + 2);  // Look two positions ahead
            if (next_static != ordered.end() && next_static->type == "static" && 
                next_static->global_class != "poly_tail") {
                int next_length = next->expected_length.value_or(1);
                positions.secondary_start = {next_static->class_id, true, -next_length - offset, ""};
                positions.secondary_stop = {next_static->class_id, true,  -next_length - 1, ""};
            } else if (prev != ordered.end() && prev->type == "static") {
                // Then try previous static
                positions.secondary_start = {prev->class_id, false, 1, ""};
                positions.secondary_stop = {prev->class_id, false, offset, ""};
            } else if (auto next_pos = position_map.find(next->class_id); next_pos != position_map.end()) {
                // Finally fall back to chaining
                int next_length = next->expected_length.value_or(1);
                positions.secondary_start = {next_pos->second.primary_start.ref_id, true, -offset - next_length, "chained_start"};
                positions.secondary_stop = {next_pos->second.primary_stop.ref_id, true, -next_length, "chained_stop"};
            }
        }

        position_map[it->class_id] = positions;
    }

    void process_reads(
        Layout_Struct::index<dir_order_tag>::type::iterator start,
        Layout_Struct::index<dir_order_tag>::type::iterator end
    ) {
        auto& ordered = by_order();
        // First find the read in this direction
        auto read_it = start;
        for (; read_it != end; ++read_it) {
            if (read_it->type == "variable" && read_it->global_class == "read") {
                break;
            }
        }
        if (read_it == end) return;
        // Process elements before the read with left-side anchoring
        for (auto it = start; it != read_it; ++it) {
            if (it->type == "variable") {
                process_left_side(it);
            }
        }
        // Process read positions
        auto prev = ordered.find(read_it->order - 1);
        if (prev != ordered.end()) {
            find_read_positions(*read_it, *prev);
        }
        // Process elements after the read with right-side anchoring
        auto after_read = read_it;
        ++after_read;
        for (auto it = after_read; it != end; ++it) {
            if (it->type == "variable") {
                process_right_side(it);
            }
        }
    }
    
    void find_read_positions(
        const ReadElement& read, 
        const ReadElement& prev) {
        auto& ordered = by_order();
        ReferencePositions positions;
        positions.primary_start = {prev.class_id, false, 1, ""};

        auto next_it = ordered.find(read.order + 1);
        auto next_next_it = ordered.find(read.order + 2);
        auto pre_prev_it = ordered.find(read.order - 2);

        if (next_it != ordered.end()) {
            positions.primary_stop = {next_it->class_id, true, -1, ""};
            if (next_next_it != ordered.end()) {
                std::string flag = next_it->type == "static" ? "truncated" : "skipped";
                int offset = next_it->expected_length.value_or(1);
                positions.secondary_stop = {
                    next_next_it->class_id,
                    true,
                    flag == "skipped" ? -(1 + offset) : -1,
                    "right_" + flag
                };
            }
        }

        if (pre_prev_it != ordered.end()) {
            std::string flag = prev.type == "static" ? "truncated" : "skipped";
            int offset = prev.expected_length.value_or(1);
            positions.secondary_start = {
                pre_prev_it->class_id,
                false,
                flag == "skipped" ? 1 + offset : 1,
                "left_" + flag
            };
        }

        position_map[read.class_id] = positions;
    }

    const auto& get_pos_map() const { 
        return position_map; 
    }

};

/**
 * @brief This class contains all of the setup for the misalignment threshold, as well as the masked adapters.
 */
class Misalignment_Setup {
public:
    struct perfect_match {
        read_streaming::sequence read;
        std::string direction; // "forward" or "reverse"
    };

    struct misalignment_stats {
    double mean;
    double sum_squares;
    size_t count;
    size_t total_reads_seen;
    bool is_stable;

    struct PositionStats {
        std::vector<double> start_positions;
        std::vector<double> stop_positions;
        
        std::pair<double, double> get_stats(const std::vector<double>& positions) {
            if (positions.empty()) return {0.0, 0.0};
            
            // Calculate mean of normalized positions
            double mean = std::accumulate(positions.begin(), positions.end(), 0.0) / positions.size();
            
            // Calculate standard deviation
            double sum_squares = std::accumulate(positions.begin(), positions.end(), 0.0,
                [mean](double sum, double val) {
                    double diff = val - mean;
                    return sum + (diff * diff);
                });
                
            double sd = std::sqrt(sum_squares / (positions.size() - 1));
            return {mean, sd};
        }
        AlignmentPositions calculate() {
            return {get_stats(start_positions), get_stats(stop_positions)};
        }
    } position_stats;

    misalignment_stats() : mean(0), sum_squares(0), count(0), total_reads_seen(0), is_stable(false) {}

    //contains defaults for the misalignment stats
    void update_stats(int new_edit_distance, size_t total_reads) {
            double old_mean = mean;
            count++;
            total_reads_seen = total_reads;
            mean = mean + (new_edit_distance - mean) / count;
            
            sum_squares += (new_edit_distance - old_mean) * (new_edit_distance - mean);
            
            double sample_ratio = static_cast<double>(count) / total_reads_seen;
            is_stable = ((std::abs(mean - old_mean) < 0.5) && (sample_ratio > 0.01)|| count >= 10000);
        }

    double get_sd() const {
            if (count < 2) return 0.0;
            return std::sqrt(sum_squares / (count - 1));
        }

    double confidence() const {
            return 1.0 - std::exp(-static_cast<double>(count) / 100.0);
        }
    
    };

    struct consensus_matrix {
        std::vector<std::map<char, int>> position_counts;  // For each position, count of each base
        int max_edit_distance;  // Maximum edit distance for sequences in this matrix
        consensus_matrix(size_t length, int max_ed)
            : position_counts(length), max_edit_distance(max_ed) {}

        void ensure_length(size_t needed_length) {
            if (needed_length > position_counts.size()) {
                position_counts.resize(needed_length + std::max(needed_length/2, size_t(10)));
            }
        }
    };

    Misalignment_Setup(const ReadLayout& layout) {
        // Extract static adapters (non poly-tail, non start/stop)
        for (const auto& elem : layout.by_type()) {
            if (elem.type == "static" &&
                elem.global_class != "poly_tail" &&
                elem.global_class != "start" &&
                elem.global_class != "stop") {
                if (elem.direction == "forward") {
                    forward_adapters.push_back({elem.class_id, elem.seq});
                    init_consensus_matrices(elem.class_id, elem.seq.length(), 20);
                } else {
                    reverse_adapters.push_back({elem.class_id, elem.seq});
                    init_consensus_matrices(elem.class_id, elem.seq.length(), 20);
                }
            }
        }
    }

    void generate_misalignment_data(const std::string& fastq_path, ReadLayout& layout, int num_threads = 1, size_t max_reads = 50000) {
        using ChunkFunc = std::function<bool(const std::vector<read_streaming::sequence>&, const std::string&)>;
        chunk_streaming<read_streaming::sequence, ChunkFunc> streamer;
        
        size_t forward_count = 0;
        size_t reverse_count = 0;
        size_t total_reads_processed = 0;
        
        //collecting statistics of both adapter stats and misalignment stats
        std::unordered_map<std::string, misalignment_stats> adapter_stats;
        std::unordered_map<std::string, misalignment_stats> misalignment_stats;
        ChunkFunc process_func = [&](const std::vector<read_streaming::sequence>& chunk, const std::string& file) -> bool {
            #pragma omp atomic
            total_reads_processed += chunk.size();
            process_misalignment_chunk(
                chunk,
                forward_count,
                reverse_count,
                adapter_stats,
                misalignment_stats,
                total_reads_processed
            );
            // Print progress and current misalignment thresholds.
            #pragma omp critical
            {
                std::cout << "\r[misalignment_stats] Processed reads - Forward perfect: " << forward_count  << ", Reverse perfect: " << reverse_count << "\n";
                for (const auto& [adapter_id, stats] : adapter_stats) {
                    std::cout << "[misalignment_stats] " << adapter_id << " misalignment mean: " << stats.mean  << " (n = " << stats.count << ")"
                            << (stats.is_stable ? " [STABLE]" : "") << "\n";
                }
            }
            // Determine if processing should stop.
            bool all_stable;
                #pragma omp critical(adapter_stats_read)
                {
                    all_stable = std::all_of(
                        adapter_stats.begin(),
                        adapter_stats.end(),
                        [](const auto& pair) { 
                            return pair.second.is_stable; 
                        }
                    );
                }
                if (all_stable || forward_count + reverse_count >= max_reads) {
                    return false;  // Signal to stop processing.
                }
                return true;  // Continue processing.
        };

        streamer.process_chunks(fastq_path, process_func, num_threads, max_reads);
        
        // Write final results
        write_perfect_matches();
        update_read_layout(layout, adapter_stats, misalignment_stats);
        layout.generate_position_mapping();
    }

private:
    std::unordered_map<std::string, std::vector<consensus_matrix>> consensus_matrices;
    std::vector<std::pair<std::string, std::string>> forward_adapters;
    std::vector<std::pair<std::string, std::string>> reverse_adapters;
    std::vector<perfect_match> perfect_matches;

    void process_misalignment_chunk(const std::vector<read_streaming::sequence>& chunk, 
        size_t& forward_count, size_t& reverse_count,
        std::unordered_map<std::string, misalignment_stats>& adapter_stats,
        std::unordered_map<std::string, misalignment_stats>& misalignment_stats,
                      size_t total_reads) {
        for (const auto& read : chunk) {

                bool perfect_forward = find_perfect_match(read.seq, forward_adapters, adapter_stats);
                bool perfect_reverse = find_perfect_match(read.seq, reverse_adapters, adapter_stats);
            
            if (perfect_forward) {
                // Update misalignment stats for reverse adapters
                update_misalignment_stats(read.seq, reverse_adapters, adapter_stats, misalignment_stats, total_reads);
                #pragma omp critical
                {
                    perfect_matches.push_back({read, "forward"});
                    #pragma omp atomic
                    forward_count++;
                }
            }
            if (perfect_reverse) {
                // Update misalignment stats for forward adapters
                update_misalignment_stats(read.seq, forward_adapters, adapter_stats, misalignment_stats, total_reads);
                #pragma omp critical
                {
                    perfect_matches.push_back({read, "reverse"});
                    #pragma omp atomic
                    reverse_count++;
                }
            }
        }
    }

    void update_misalignment_stats(const std::string& read_seq, const std::vector<std::pair<std::string, std::string>>& adapters,
                       std::unordered_map<std::string, misalignment_stats>& adapter_stats,
                       std::unordered_map<std::string, misalignment_stats>& misalignment_stats,
                       size_t total_reads) {
        for (const auto& [id, seq] : adapters) {
            EdlibAlignResult result = edlibAlign(
                seq.c_str(), seq.length(),
                read_seq.c_str(), read_seq.length(),
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0)
            );
            std::string aligned_seq = collect_aligned_seq(read_seq, seq, result);
            //something crashes every once in a while during original generation
            #pragma omp critical
            {
                adapter_stats[id].update_stats(result.editDistance, total_reads);
                //this calculates the misalignment positions for the reverse primers in forward locations, and vice versa
                    double norm_misal_start = (static_cast<double>(result.startLocations[0]) / read_seq.length()) * 100.0;
                    double norm_misal_stop = (static_cast<double>(result.endLocations[0]) / read_seq.length()) * 100.0;
                    misalignment_stats[id].position_stats.start_positions.push_back(norm_misal_start);
                    misalignment_stats[id].position_stats.stop_positions.push_back(norm_misal_stop);

                update_consensus_matrices(id, aligned_seq, result.editDistance);
            }
            edlibFreeAlignResult(result);
        }
    }

    void write_perfect_matches(bool write = false) {
        std::string home = std::getenv("HOME");
        std::string desktop = home + "/Desktop/";
        if(write){
        std::ofstream forward_file(desktop + "perfect_forward_reads.fastq");
        std::ofstream reverse_file(desktop + "perfect_reverse_reads.fastq");
        
        if (!forward_file.is_open() || !reverse_file.is_open()) {
            throw std::runtime_error("Failed to open output files");
        }

        for (const auto& match : perfect_matches) {
            std::ofstream& out_file = (match.direction == "forward") ? 
                                    forward_file : reverse_file;
            
            out_file << "@" << match.read.id;
            if (!match.read.comment.empty()) {
                out_file << " " << match.read.comment;
            }
            out_file << "\n"
                    << match.read.seq << "\n"
                    << "+\n"
                    << match.read.qual << "\n";
        }

        std::cout << "Perfect match reads written to:\n"
                  << desktop + "perfect_forward_reads.fastq\n"
                  << desktop + "perfect_reverse_reads.fastq\n";
        }
    }

    bool find_perfect_match(const std::string& read_seq,const std::vector<std::pair<std::string, std::string>>& adapters,
                            std::unordered_map<std::string, misalignment_stats>& adapter_stats) {
        for (const auto& [id, seq] : adapters) {
            EdlibAlignResult result = edlibAlign(
                seq.c_str(), seq.length(),
                read_seq.c_str(), read_seq.length(),
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0)
            );
            bool is_perfect = (result.editDistance == 0 && result.status == EDLIB_STATUS_OK);
            if(is_perfect){
                double normalized_start = (static_cast<double>(result.startLocations[0]) / read_seq.length()) * 100.0;
                double normalized_stop = (static_cast<double>(result.endLocations[0]) / read_seq.length()) * 100.0;
                #pragma omp critical
                {
                    adapter_stats[id].position_stats.start_positions.push_back(normalized_start);
                    adapter_stats[id].position_stats.stop_positions.push_back(normalized_stop);    
                }
            }
            edlibFreeAlignResult(result);
            if(!is_perfect) return false;
        }
        return true;
    }

    void init_consensus_matrices(const std::string& adapter_id, size_t adapter_length, int max_errors) {
        consensus_matrices[adapter_id].clear();
        for(int ed = 0; ed <= max_errors; ed++) {
            consensus_matrices[adapter_id].emplace_back(adapter_length, ed);
        }
    }

    void update_consensus_matrices(const std::string& adapter_id, const std::string& aligned_seq, int edit_distance) {
            // Add sequence to all matrices that accept this edit distance or higher
            for(auto& matrix : consensus_matrices[adapter_id]) {
                if(edit_distance <= matrix.max_edit_distance) {
                    matrix.ensure_length(aligned_seq.length());
                    for(size_t pos = 0; pos < aligned_seq.length(); pos++) {
                        matrix.position_counts[pos][aligned_seq[pos]]++;
                    }
                }
            }
        }

    std::string collect_aligned_seq(const std::string& read_seq, const std::string& adapter_seq, EdlibAlignResult& result) {
        // Get CIGAR string
        char* cigar = edlibAlignmentToCigar(
            result.alignment,
            result.alignmentLength,
            EDLIB_CIGAR_EXTENDED
        );
        
        // Extract aligned portion
        int start_pos = result.startLocations[0];
        int end_pos = result.endLocations[0] + 1;
        std::string aligned = read_seq.substr(start_pos, end_pos - start_pos);
        
        free(cigar);
        return aligned;
    }

    std::string generate_masked_adapter(const std::string& adapter_id, const std::string adapter_seq,  const misalignment_stats& stats) {
    int threshold = static_cast<int>(std::floor(stats.mean - stats.get_sd()));
    // Define max consecutive N's allowed as (threshold - 1), with a minimum of 1.
    int max_consecutive_N = (threshold > 1 ? threshold - 1 : 1);
    int max_total_N = 12;
    const auto& matrices = consensus_matrices[adapter_id];
    
    // Get matrix where max_edit_distance is >= threshold.
    auto it = std::find_if(matrices.begin(), matrices.end(),
        [threshold](const consensus_matrix& matrix) {
            return matrix.max_edit_distance >= threshold;
        });
    
    if(it == matrices.end()) return "";
    
    // Build a consensus matrix for each position.
    std::vector<std::vector<double>> consensus_matrix(it->position_counts.size(), std::vector<double>(4, 0.0));
    std::string masked_seq;
    
    for (size_t pos = 0; pos < it->position_counts.size() && pos < adapter_seq.length(); pos++) {
        int total_count = 0;
        double max_percentage = 0.0;
        char max_base = 'N';
        char original_base = adapter_seq[pos];
        const char bases[] = {'A', 'C', 'G', 'T'};
    
        for (const auto& [base, count] : it->position_counts[pos]) {
            total_count += count;
        }
        // Calculate percentages and find dominant base.
        for (const auto& [base, count] : it->position_counts[pos]) {
            double percentage = (total_count > 0) ? static_cast<double>(count) / total_count * 100.0 : 0.0;
            size_t idx;
            switch(base) {
                case 'A': idx = 0; break;
                case 'C': idx = 1; break;
                case 'G': idx = 2; break;
                case 'T': idx = 3; break;
                default: continue;
            }
            consensus_matrix[pos][idx] = percentage;
            if (percentage > max_percentage) {
                max_percentage = percentage;
                max_base = base;
            }
        }
        // Decide the masked character for this position.
        // masked_seq += (max_percentage > 50.0) ? original_base : 'N';
        // if the character is 'N', enforce a maximum run length.
        char masked_char = (max_percentage > 50.0) ? original_base : 'N';
        if (masked_char == 'N') {
            // Count how many consecutive 'N's are already in masked_seq.
            int n_run = 0;
            for (int i = static_cast<int>(masked_seq.size()) - 1; i >= 0; i--) {
                if (masked_seq[i] == 'N')
                    n_run++;
                else
                    break;
            }
            // Also count the total number of 'N's in the masked sequence.
            int total_N_count = std::count(masked_seq.begin(), masked_seq.end(), 'N');
            // If we already have too many consecutive 'N's, or too many 'N's in total,
            // then use the original base instead.
            if (n_run >= max_consecutive_N || total_N_count >= max_total_N) {
                masked_char = original_base;
                //continue;
            }
        }
        masked_seq.push_back(masked_char);
    }
    return masked_seq;
}

    void update_read_layout(
        ReadLayout& layout,
         std::unordered_map<std::string, misalignment_stats>& adapter_stats,
         std::unordered_map<std::string, misalignment_stats>& misalignment_stats) {
        auto& by_id_index = layout.by_id();
        for(auto& [adapter_id, stats] : adapter_stats) {
            auto it = by_id_index.find(adapter_id);
            if(it != by_id_index.end()) {
                std::string seq = it->seq;
                std::string masked = generate_masked_adapter(adapter_id, seq, stats);
                
                double sd = stats.get_sd();
                std::tuple<int, int, int> threshold{
                    static_cast<int>(std::round(stats.mean - sd)),
                    static_cast<int>(std::floor(stats.mean)),
                    static_cast<int>(std::round(stats.mean + sd))
                };
                
                AlignmentPositions align_pos = stats.position_stats.calculate();
                AlignmentPositions misalign_pos = misalignment_stats[adapter_id].position_stats.calculate();

                by_id_index.modify(it, [&](ReadElement& elem) {
                    elem.masked_seq = masked;
                    elem.misalignment_threshold = threshold;
                    elem.aligned_positions = align_pos;
                    elem.misaligned_positions= misalign_pos;
                });
            }
        }
    }

};
