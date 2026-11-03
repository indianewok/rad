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
    std::cout << "PASS native whitelist loader: " << checked << " format-equivalent calls, role/format/pooling checks\n";
}
