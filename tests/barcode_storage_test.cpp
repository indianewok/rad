#include "rad/rad_headers.h"
#include <cstdlib>
#include <random>
#include <sstream>

static void require(bool pass, const char* message) {
    if (!pass) { std::cerr << message << '\n'; std::exit(1); }
}
static barcode_entry entry(const std::string& text) {
    barcode_entry result;
    result.barcode.sequence_to_bits(text);
    return result;
}
static std::string random_dna(std::mt19937_64& generator, size_t length) {
    std::string result(length, 'A');
    for (char& base : result) base = "ACTG"[generator() % 4];
    return result;
}
static int infix_distance(const std::string& pattern, const std::string& text) {
    std::vector<int> row(pattern.size() + 1);
    std::iota(row.begin(), row.end(), 0);
    int best = static_cast<int>(pattern.size());
    for (const char base : text) {
        std::vector<int> next(pattern.size() + 1, 0);
        for (size_t i = 1; i <= pattern.size(); ++i)
            next[i] = std::min({row[i] + 1, next[i - 1] + 1,
                                row[i - 1] + (pattern[i - 1] != base)});
        best = std::min(best, next.back());
        row.swap(next);
    }
    return best;
}
static int sequence_distance(const std::string& a, const std::string& b) {
    return a.size() <= b.size() ? infix_distance(a, b) : infix_distance(b, a);
}
int main() {
    static_assert(std::tuple_size<decltype(counter{}.counts)>::value == 3,
                  "Whitelist entries should store only three counters");
    static_assert(sizeof(counter) == 3 * sizeof(std::atomic<int>),
                  "Counter storage should have no removed reporting slots");
    std::mt19937_64 generator(918274);
    const auto short_text = random_dna(generator, 25);
    const auto wide_text = random_dna(generator, 33);
    const auto full_text = random_dna(generator, 64);
    const int64_seq short_key(short_text), wide_key(wide_text), full_key(full_text);
    bc_flat_set catalog;
    catalog.insert_bc_entry(short_key, entry(short_text));
    require(catalog.narrow_size() == 1 && catalog.wide_size() == 0 && !catalog.has_wide_storage(),
            "short catalog allocated a wide table");
    const int64_seq window("ACTG" + short_text + "TGCA");
    require(mutation_tools::int64_lvdist(short_key, window, 0) == 0,
            "25-base catalog key cannot match 33-base temporary window");
    require(!short_key.is_wide() && !catalog.has_wide_storage(),
            "temporary promotion changed stored width");
    catalog.insert_bc_entry(wide_key, entry(wide_text));
    catalog.insert_bc_entry(full_key, entry(full_text));
    require(catalog.narrow_size() == 1 && catalog.wide_size() == 2 && catalog.has_wide_storage(),
            "catalog did not route wide keys separately");
    auto changed_high = full_text;
    changed_high[0] = changed_high[0] == 'A' ? 'C' : 'A';
    require(!catalog.check_wl_for(int64_seq(changed_high)), "wide lookup ignored its upper bits");
    catalog.insert_bc_entry(int64_seq(changed_high), entry(changed_high));
    require(catalog.wide_size() == 3, "wide keys with equal low halves collapsed");
    for (const auto& key : {short_key, wide_key, full_key}) {
        catalog.set_bc_count(key, 3, barcode_counts::raw);
        catalog.update_bc_count(key, barcode_counts::raw);
        catalog.subtract_bc_count(key, barcode_counts::raw);
        require(catalog.get_bc_count(key, barcode_counts::raw) == 3,
                "typed catalog counter operations changed");
    }
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 2000; ++i) {
        catalog.update_bc_count(short_key);
        catalog.update_bc_count(wide_key);
    }
    require(catalog.get_bc_count(short_key) == 2000 && catalog.get_bc_count(wide_key) == 2000,
            "parallel typed catalog counters lost updates");

    bc_map aliases;
    aliases.insert_bc_entry(short_key, entry(wide_text));
    aliases.insert_bc_entry(wide_key, entry(short_text));
    require(aliases.return_putative_correct_bcs(short_key) == std::unordered_set<int64_seq>{wide_key},
            "narrow observed alias did not resolve to wide canonical key");
    require(aliases.return_putative_correct_bcs(wide_key) == std::unordered_set<int64_seq>{short_key},
            "wide observed alias did not resolve to narrow canonical key");
    require(!aliases.has_identity_mapping(short_key) && !aliases.has_identity_mapping(wide_key),
            "mutation alias was treated as an exact identity");
    aliases.insert_bc_entry(short_key, entry(short_text));
    require(aliases.has_identity_mapping(short_key), "exact identity was lost beside an alias");
    require(aliases.unique_val_size() == 2, "cross-width aliases duplicated canonical values");
    auto copied = aliases;
    copied.set_bc_count(wide_key, 7);
    require(aliases.get_bc_count(wide_key) == 0, "copied alias map shares source counters");
    require(copied.return_putative_correct_bcs(short_key).count(wide_key), "copy lost cross-width alias");
    bc_map concurrent;
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 2000; ++i)
        concurrent.insert_bc_entry(i % 2 ? short_key : full_key, entry(wide_text));
    require(concurrent.unique_val_size() == 1 && concurrent.size() == 2,
            "parallel cross-width insertion duplicated canonical entries");
    require(concurrent.return_putative_correct_bcs(full_key).count(wide_key),
            "parallel insertion lost wide alias association");

    whitelist::wl_entry indexes;
    std::vector<std::string> reference;
    for (const int length : {8, 16, 25, 32, 33, 40, 63, 64}) {
        for (int n = 0; n < 3; ++n) {
            auto text = random_dna(generator, length);
            indexes.true_bcs.insert_bc_entry(int64_seq(text), entry(text));
            reference.push_back(std::move(text));
        }
    }
    indexes.build_seed_idx(false);
    indexes.build_part_idx(false);
    for (const auto& text : reference) {
        const int64_seq packed(text);
        const auto candidates = indexes.query_bc_seeds(packed, "true");
        require(std::any_of(candidates.begin(), candidates.end(), [&](auto handle) {
                    return handle.barcode() == packed;
                }), "mixed-length seed index lost exact wide barcode");
        for (int variation = 0; variation < 3; ++variation) {
            auto query = text;
            if (variation == 1) query[text.size() / 2] = query[text.size() / 2] == 'A' ? 'C' : 'A';
            if (variation == 2) query.erase(query.size() / 2, 1);
            if (query.size() + 8 <= 64) query = "ACTG" + query + "TGCA";
            for (int cap = 0; cap <= 2; ++cap) {
                std::vector<std::pair<int64_seq, int>> matches;
                indexes.true_parts.scan(int64_seq(query), cap, matches);
                std::map<std::string, int> actual, expected;
                for (const auto& match : matches) actual.emplace(match.first.bits_to_sequence(), match.second);
                for (const auto& target : reference) {
                    const int distance = sequence_distance(query, target);
                    if (distance <= cap) expected.emplace(target, distance);
                }
                require(actual == expected, "mixed-width part index differs from independent string DP");
            }
        }
    }
    std::ostringstream summary;
    catalog.write_wl_summary(summary, "barcode");
    require(summary.str().find(full_text) != std::string::npos &&
            summary.str().find(short_text) != std::string::npos,
            "summary did not decode both resident widths");

    std::istringstream exported(summary.str());
    std::string line;
    require(bool(std::getline(exported, line)) &&
            line == "class_id,true_barcode,raw_count,total_count,corrected_count",
            "whitelist summary retained removed counter columns");
    while (std::getline(exported, line)) {
        std::istringstream row(line);
        std::vector<std::string> fields;
        for (std::string field; std::getline(row, field, ',');) fields.push_back(field);
        require(fields.size() == 5, "whitelist summary has unexpected counter fields");
        const int64_seq key(fields[1]);
        require(std::stoi(fields[2]) == catalog.get_bc_count(key, barcode_counts::raw) &&
                std::stoi(fields[3]) == catalog.get_bc_count(key, barcode_counts::total) &&
                std::stoi(fields[4]) == catalog.get_bc_count(key, barcode_counts::corrected),
                "remaining counter columns changed meaning");
    }
    whitelist::wl_entry quality;
    quality.global_bcs.insert_bc_entry(short_key, entry(short_text));
    quality.global_bcs.set_bc_count(short_key, 1, barcode_counts::raw);
    quality.global_bcs.set_bc_count(short_key, 10, barcode_counts::total);
    quality.global_bcs.set_bc_count(short_key, 10, barcode_counts::corrected);
    require(!barcode_correction::passes_quality_check(short_key, &quality, "global", "offensive", false) &&
            quality.global_bcs.get_bc_count(short_key, barcode_counts::total) == -1,
            "three-counter QC lost its signed rejection state");
    quality.global_bcs.set_bc_count(short_key, 2, barcode_counts::raw);
    require(barcode_correction::passes_quality_check(short_key, &quality, "global", "offensive", false) &&
            quality.global_bcs.get_bc_count(short_key, barcode_counts::total) == 2 &&
            quality.global_bcs.get_bc_count(short_key, barcode_counts::corrected) == 10,
            "QC recovery reset the wrong counter");
    quality.global_bcs.set_bc_count(short_key, 10, barcode_counts::total);
    quality.global_bcs.set_bc_count(short_key, 9, barcode_counts::corrected);
    require(!barcode_correction::passes_quality_check(short_key, &quality, "global", "offensive", false),
            "three-counter QC ignored a high correction ratio");
    quality.global_bcs.set_bc_count(short_key, 7, barcode_counts::corrected);
    require(barcode_correction::passes_quality_check(short_key, &quality, "global", "offensive", false),
            "three-counter QC rejected the passing correction ratio");

    const std::string numeric_path = "/tmp/rad-wide-numeric-" + std::to_string(getpid()) + ".tsv";
    const std::string mixed_path = "/tmp/rad-wide-mixed-" + std::to_string(getpid()) + ".tsv";
    {
        std::ofstream numeric(numeric_path), mixed(mixed_path);
        numeric << "barcode\n340282366920938463463374607431768211455\n0\n";
        mixed << short_text << '\n' << wide_text << '\n' << full_text << '\n';
    }
    whitelist loader;
    auto imported = loader.import_whitelist(numeric_path, false, 64);
    require(imported.true_bcs.narrow_size() == 0 && imported.true_bcs.wide_size() == 2 &&
            imported.true_bcs.has_identity_mapping(int64_seq(std::string(64, 'G'))) &&
            imported.true_bcs.has_identity_mapping(int64_seq(std::string(64, 'A'))),
            "numeric whitelist loader lost full unsigned128 identities");
    auto scan = load_scan_whitelist(numeric_path, 64, false);
    require(scan.long_barcodes.narrow.empty() && scan.long_barcodes.wide.size() == 2 &&
            scan.long_barcodes.contains(int64_seq(std::string(64, 'G'))),
            "scan whitelist loader lost unsigned128 identities");
    auto mixed = loader.import_whitelist(mixed_path, false, 25);
    require(mixed.true_bcs.narrow_size() == 1 && mixed.true_bcs.wide_size() == 2,
            "mixed-width DNA catalog chose one width for all entries");
    scan = load_scan_whitelist(mixed_path, 25, false);
    require(scan.size() == 3 && scan.packed_barcodes.size() == 1 && scan.long_barcodes.wide.size() == 2,
            "scan whitelist did not preserve the short sentinel path beside wide keys");
    std::remove(numeric_path.c_str());
    std::remove(mixed_path.c_str());

    ReadLayout layout;
    whitelist::wl_entry only;
    only.true_bcs.insert_bc_entry(short_key, entry(short_text));
    layout.wl_map.lists.emplace("test", std::move(only));
    layout.wl_map.maps.emplace("barcode", std::ref(layout.wl_map.lists.at("test")));
    const std::string read_text = "ACGTA" + short_text + "TGC";
    seq_element element("barcode", "barcode", std::nullopt, {5, 29}, "variable", 1, "forward");
    element.seq = read_text.substr(4, 25);
    read_streaming::sequence read{"shift", "", read_text, "", false};
    barcode_correction::bc_position position;
    const auto corrected = barcode_correction::correct_barcode(element, layout, read, false, 2,
                                                               "offensive", &position);
    require(corrected && *corrected == short_key && position.kind == barcode_correction::bc_position::located &&
            position.end == 29 && position.shift == 1,
            "expanded 33-base window did not refine the barcode position");
    require(layout.wl_map["barcode"].true_bcs.wide_size() == 0,
            "position refinement widened resident catalog");
    spat_mask mask;
    mask.rows = mask.cols = 1;
    mask.bits = {true};
    mask.bc1_to_idx.set(short_key, 0);
    mask.bc2_to_idx.set(short_key, 0);
    require(mask.check(short_key, short_key) && !mask.bc1_to_idx.has_wide_storage(),
            "short spatial mask axis initialized wide storage");
    mask.bc1_to_idx.set(wide_key, 0);
    require(mask.check(wide_key, short_key) && mask.bc1_to_idx.has_wide_storage(),
            "spatial mask lost mixed-width coordinates");
    std::cout << "Resident key bytes: " << sizeof(packed_barcode<uint64_t>) << "/"
              << sizeof(packed_barcode<uint128_t>) << "; entry bytes: "
              << sizeof(basic_barcode_entry<uint64_t>) << "/"
              << sizeof(basic_barcode_entry<uint128_t>) << '\n';
    std::cout << "Typed storage, mixed-width indexes, aliases, counters and position refinement passed\n";
}
