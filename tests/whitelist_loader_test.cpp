#include "rad/rad_headers.h"
#include <cstdlib>
#include <iostream>
#include <random>

static void require(bool pass, const char* message) {
    if (!pass) { std::cerr << message << '\n'; std::exit(1); }
}
static std::string dna(uint64_t word, unsigned length = 25) {
    return int64_seq(static_cast<int64_t>(word), length).bits_to_sequence();
}
static void bind(ReadLayout& layout, whitelist::wl_entry&& entry) {
    layout.wl_map.lists.emplace("catalog", std::move(entry));
    layout.wl_map.maps.emplace("barcode", std::ref(layout.wl_map.lists.at("catalog")));
}
static std::optional<int64_seq> call(const ReadLayout& layout, const std::string& raw,
                                   const std::string& expanded, int cap, bool reverse,
                                   const std::string& mode) {
    seq_element elem("barcode", "barcode", std::nullopt, {5, 4 + static_cast<int>(raw.size())},
                     "variable", 1, reverse ? "reverse" : "forward");
    elem.seq = reverse ? seq_utils::revcomp(raw) : raw;
    read_streaming::sequence read{"query", "", expanded, "", false};
    return barcode_correction::correct_barcode_lookup(elem, layout, read, expanded, false, cap, mode);
}
static std::unordered_set<int64_seq> identities(const whitelist::wl_entry& entry) {
    std::unordered_set<int64_seq> result;
    for (const auto e : entry.global_bcs.get_unique_entries()) result.insert(e.barcode());
    for (const auto e : entry.true_bcs.get_unique_entries()) result.insert(e.barcode());
    return result;
}
template <class Table>
static std::unordered_set<int64_seq> resident_identities(const Table& table) {
    std::unordered_set<int64_seq> result;
    for (const auto e : table.get_unique_entries()) result.insert(e.barcode());
    return result;
}
template <class Table>
static bool zero_counts(const Table& table) {
    for (const auto e : table.get_unique_entries()) {
        const auto* counts = e.counts();
        if (!counts || counts->load(raw) || counts->load(total) || counts->load(corrected)) return false;
    }
    return true;
}
static void write_fixture(const std::string& path, const std::vector<std::string>& lines) {
    if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
        gzFile output = gzopen(path.c_str(), "wb");
        require(output != nullptr, "cannot create dual-list gzip fixture");
        for (const auto& line : lines) {
            const std::string row = line + '\n';
            require(gzwrite(output, row.data(), static_cast<unsigned>(row.size())) == static_cast<int>(row.size()),
                    "dual-list gzip fixture write failed");
        }
        require(gzclose(output) == Z_OK, "dual-list gzip fixture close failed");
    } else {
        std::ofstream output(path);
        require(bool(output), "cannot create dual-list fixture");
        for (const auto& line : lines) output << line << '\n';
    }
}
static size_t check_dual(const std::string& label, const std::string& first, const std::string& second,
                         uint16_t length, const std::vector<std::string>& true_sequences,
                         const std::vector<std::string>& global_sequences) {
    whitelist loader;
    std::vector<std::unordered_set<int64_seq>> legacy_sets;
    legacy_sets.push_back(whitelist::load_barcodes(first, length, false));
    legacy_sets.push_back(whitelist::load_barcodes(second, length, false));
    whitelist::wl_entry legacy;
    loader.populate_entry(legacy, std::move(legacy_sets));
    auto imported = loader.import_whitelist(first + ":" + second, false, length);
    std::unordered_set<int64_seq> expected_true, expected_global;
    for (const auto& sequence : true_sequences) expected_true.emplace(sequence);
    for (const auto& sequence : global_sequences) expected_global.emplace(sequence);
    const std::string failure = "dual-list semantics changed: " + label;
    require(resident_identities(imported.true_bcs) == expected_true &&
            resident_identities(imported.global_bcs) == expected_global, failure.c_str());
    require(resident_identities(imported.true_bcs) == resident_identities(legacy.true_bcs) &&
            resident_identities(imported.global_bcs) == resident_identities(legacy.global_bcs), failure.c_str());
    require(imported.true_bcs.size() == expected_true.size() &&
            imported.global_bcs.size() == expected_global.size(), failure.c_str());
    require(zero_counts(imported.true_bcs) && zero_counts(imported.global_bcs) &&
            zero_counts(legacy.true_bcs) && zero_counts(legacy.global_bcs),
            "dual-list import changed initial counters");

    ReadLayout current_layout, legacy_layout;
    bind(current_layout, std::move(imported));
    bind(legacy_layout, std::move(legacy));
    std::vector<std::pair<std::string, std::string>> probes;
    for (const auto& sequence : true_sequences) {
        // The width-limit fixtures verify storage without constructing windows
        // longer than the supported 64-base transient representation.
        if (sequence.size() > 56) continue;
        probes.emplace_back(sequence, "ACGT" + sequence + "TGCA");
        if (sequence.size() == 25) {
            std::string mutated = sequence;
            mutated[10] = sequence[10] == 'A' ? 'C' : 'A';
            probes.emplace_back(mutated, "ACGT" + mutated + "TGCA");
            const auto deleted = sequence.substr(0, 10) + sequence.substr(11);
            probes.emplace_back(deleted, "ACGT" + deleted + "TGCA");
        }
    }
    for (const auto& sequence : global_sequences)
        if (sequence.size() <= 56) probes.emplace_back(sequence, "ACGT" + sequence + "TGCA");
    if (probes.empty()) probes.emplace_back(std::string(25, 'C'), "ACGT" + std::string(25, 'C') + "TGCA");
    size_t checked = 0;
    for (const auto& probe : probes)
        for (int cap : {0, 1, 2})
            for (bool reverse : {false, true})
                for (const std::string mode : {"offensive", "defensive"}) {
                    require(call(current_layout, probe.first, probe.second, cap, reverse, mode) ==
                            call(legacy_layout, probe.first, probe.second, cap, reverse, mode), failure.c_str());
                    ++checked;
                }
    return checked;
}
int main() {
    char directory[] = "/tmp/rad-whitelist-loader-XXXXXX";
    require(mkdtemp(directory) != nullptr, "cannot make test directory");
    const std::string base(directory);
    const std::string plain = base + "/catalog.tsv", gzip = plain + ".gz";
    const std::string numeric = base + "/catalog.bits.tsv";
    std::mt19937_64 random(123456789);
    std::unordered_set<uint64_t> unique;
    while (unique.size() < 100005) unique.insert(random() & ((uint64_t(1) << 50) - 1));
    std::vector<std::string> reference;
    for (uint64_t word : unique) reference.push_back(dna(word));
    std::ofstream out(plain), bitout(numeric);
    gzFile compressed = gzopen(gzip.c_str(), "wb");
    require(bool(out) && bool(bitout) && compressed, "cannot create catalog fixtures");
    for (const auto& sequence : reference) {
        const std::string line = sequence + "\t10\t20\n";
        out << line;
        const int64_seq packed(sequence);
        bitout << static_cast<int64_t>(packed.as_uint64()) << '\n';
        require(gzwrite(compressed, line.data(), static_cast<unsigned>(line.size())) == static_cast<int>(line.size()),
                "gzip fixture write failed");
    }
    // Duplicate rows must not change the source role or retained identities.
    out << reference.front() << '\n';
    out.close(); bitout.close();
    require(gzclose(compressed) == Z_OK, "gzip fixture close failed");

    std::vector<std::pair<std::string, std::string>> queries;
    for (size_t i = 0; i < 12; ++i) {
        const auto& target = reference[i];
        queries.emplace_back(target, "ACGT" + target + "TGCA");
        for (int edits = 1; edits <= 3; ++edits) {
            std::string mutated = target;
            for (int j = 0; j < edits; ++j) mutated[4 + 6 * j] = target[4 + 6 * j] == 'A' ? 'C' : 'A';
            queries.emplace_back(mutated, "ACG" + mutated + "TGC");
            if (edits == 1) queries.emplace_back(mutated, "ACGT" + mutated + "TGCA");
        }
        // An indel/positional error can expose an exact barcode in the expanded
        // window, even when the extracted raw barcode is a shifted sequence.
        queries.emplace_back(target.substr(1) + "T", "ACG" + target + "TGCA");
        queries.emplace_back(seq_utils::revcomp(target), "ACGT" + seq_utils::revcomp(target) + "TGCA");
    }
    whitelist loader;
    ReadLayout full;
    bind(full, loader.import_whitelist(plain, false, 25));
    require(full.wl_map["barcode"].global_bcs.size() == unique.size() &&
            full.wl_map["barcode"].true_bcs.empty(), "full source fixture role or size changed");
    std::string one_edit = reference[0];
    one_edit[4] = one_edit[4] == 'A' ? 'C' : 'A';
    require(call(full, one_edit, "ACG" + one_edit + "TGC", 1, false, "offensive") ==
            std::optional<int64_seq>(int64_seq(reference[0])),
            "fixture must exercise a successful native one-edit correction");
    size_t checked = 0;
    const auto expected_ids = identities(full.wl_map["barcode"]);
    for (const auto& path : {gzip, numeric}) {
        ReadLayout formatted;
        bind(formatted, loader.import_whitelist(path, false, 25));
        require(formatted.wl_map["barcode"].true_bcs.empty() &&
                formatted.wl_map["barcode"].global_bcs.size() == unique.size(),
                "catalog format changed the source role or size");
        require(identities(formatted.wl_map["barcode"]) == expected_ids,
                "catalog format changed imported identities");
        for (int cap : {0, 1, 2, 3}) {
            for (const auto& q : queries) {
                for (bool reverse : {false, true}) for (const std::string mode : {"offensive", "defensive"}) {
                    require(call(full, q.first, q.second, cap, reverse, mode) ==
                            call(formatted, q.first, q.second, cap, reverse, mode),
                            "native correction differs between catalog formats");
                    ++checked;
                }
            }
        }
    }

    // A small catalog remains a true whitelist, including indel targets.
    const std::string small = base + "/small.tsv";
    { std::ofstream f(small); f << reference[0] << '\n' << reference[1] << '\n'; }
    auto small_entry = loader.import_whitelist(small, false, 25);
    require(small_entry.global_bcs.empty() && small_entry.true_bcs.size() == 2,
            "small true whitelist role or size changed");
    ReadLayout small_full;
    bind(small_full, std::move(small_entry));
    for (const auto& raw : {reference[0].substr(0, 10) + reference[0].substr(11),
                            reference[0].substr(0, 10) + "A" + reference[0].substr(10)}) {
        const std::string expanded = "AC" + raw + "GT";
        require(call(small_full, raw, expanded, 2, false, "offensive") ==
                std::optional<int64_seq>(int64_seq(reference[0])),
                "small true-list indel fallback changed");
    }
    const std::string duplicates = base + "/duplicates.tsv";
    { std::ofstream f(duplicates); for (int i = 0; i < 100001; ++i) f << reference[0] << '\n'; }
    auto duplicate_entry = loader.import_whitelist(duplicates, false, 25);
    require(duplicate_entry.global_bcs.empty() && duplicate_entry.true_bcs.size() == 1,
            "row count incorrectly replaced source unique cardinality");

    // The packed dual-input path must retain the legacy parser and role rules,
    // independently of row counts, input order, headers, or sequence width.
    const std::string dual_a = base + "/dual_a.tsv", dual_b = base + "/dual_b.tsv.gz";
    size_t dual_checked = 0;
    write_fixture(dual_a, {reference[0], reference[1]});
    write_fixture(dual_b, {reference[0], reference[2], reference[3], reference[4]});
    dual_checked += check_dual("unequal-overlap", dual_a, dual_b, 25,
                               {reference[0], reference[1]}, {reference[2], reference[3], reference[4]});
    dual_checked += check_dual("unequal-reversed", dual_b, dual_a, 25,
                               {reference[0], reference[1]}, {reference[2], reference[3], reference[4]});
    write_fixture(dual_a, std::vector<std::string>(1001, reference[0]));
    write_fixture(dual_b, {reference[1], reference[2], reference[2]});
    dual_checked += check_dual("duplicate-cardinality", dual_a, dual_b, 25,
                               {reference[0]}, {reference[1], reference[2]});
    write_fixture(dual_a, {reference[0], reference[1], reference[0]});
    write_fixture(dual_b, {reference[1], reference[2]});
    dual_checked += check_dual("equal-union", dual_a, dual_b, 25,
                               {reference[0], reference[1], reference[2]}, {});
    write_fixture(dual_a, {});
    write_fixture(dual_b, {reference[0], reference[1]});
    dual_checked += check_dual("empty-smaller", dual_a, dual_b, 25, {}, {reference[0], reference[1]});
    dual_checked += check_dual("empty-reversed", dual_b, dual_a, 25, {}, {reference[0], reference[1]});
    write_fixture(dual_b, {});
    dual_checked += check_dual("both-empty", dual_a, dual_b, 25, {}, {});
    write_fixture(dual_a, {"whitelist_bcs", "NNNNN", "ACGTN", reference[0] + "\tmetadata"});
    write_fixture(dual_b, {reference[1], reference[2]});
    dual_checked += check_dual("invalid-sentinel-cardinality", dual_a, dual_b, 25,
                               {reference[0], reference[1], reference[2]}, {});
    const std::string thirty_two = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
    const std::string thirty_three = "T" + thirty_two;
    const std::string sixty_four = thirty_two + thirty_two;
    write_fixture(dual_a, {reference[0], thirty_two, thirty_three});
    write_fixture(dual_b, {reference[2], reference[3], thirty_two, sixty_four});
    dual_checked += check_dual("mixed-native-widths", dual_a, dual_b, 25,
                               {reference[0], thirty_two, thirty_three}, {reference[2], reference[3], sixty_four});
    // With default length zero, numeric tokens remain distinct parsed keys but
    // are omitted from resident entries. Their cardinality still selects roles.
    write_fixture(dual_a, {"1", "2", "3", "3"});
    write_fixture(dual_b, {reference[0], reference[1]});
    dual_checked += check_dual("multiple-length-zero-numeric-keys", dual_a, dual_b, 0,
                               {reference[0], reference[1]}, {});
    dual_checked += check_dual("multiple-length-zero-numeric-keys-reversed", dual_b, dual_a, 0,
                               {reference[0], reference[1]}, {});

    // Distinct classes and their RC aliases share one fully imported source.
    ReadLayout shared;
    shared.layout.insert(ReadElement("barcode", "", "", 25, "variable", 1, "forward", "barcode", plain));
    shared.layout.insert(ReadElement("rc_barcode", "", "", 25, "variable", 2, "reverse", "barcode", plain));
    shared.layout.insert(ReadElement("barcode2", "", "", 25, "variable", 3, "forward", "barcode", plain));
    shared.load_wl(std::nullopt, 0, false, 1);
    require(shared.wl_map.lists.size() == 1, "shared source imported more than once");
    require(shared.wl_map["barcode"].global_bcs.size() == unique.size() &&
            identities(shared.wl_map["barcode2"]) == expected_ids,
            "shared source did not retain the full catalog");
    require(&shared.wl_map["barcode"] == &shared.wl_map["barcode2"], "shared catalog aliases diverged");
    boost::filesystem::remove_all(base);
    std::cout << "PASS native whitelist loader: " << checked << " format-equivalent calls, "
              << dual_checked << " legacy-equivalent dual-list calls, role/format/pooling checks\n";
}
