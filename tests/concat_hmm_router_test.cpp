#include "rad/sigstring.hpp"
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

concat_hmm::layout_spec curio() {
    concat_hmm::layout_spec s;
    s.name = "router one-sided native capture test";
    auto add = [&](const char* id, const char* seq, bool fixed, const char* klass, int length) {
        concat_hmm::layout_element e;
        e.id = id; e.seq = seq; e.is_static = fixed; e.klass = klass;
        e.order = (int)s.elements.size() + 1;
        if (length) e.length_candidates = {length};
        s.elements.push_back(e);
    };
    add("barcode_1", "", false, "barcode", 8);
    add("forw_primer", "TCTTCAGCGTTCCCGAGA", true, "forw_primer", 18);
    add("barcode_2", "", false, "barcode", 6);
    add("umi", "", false, "umi", 7);
    add("poly_t", "T{8,}+", true, "poly_tail", 8);
    add("read", "", false, "read", 0);
    return s;
}

std::string dna(int count, unsigned seed) {
    std::string s;
    for (int i = 0; i < count; ++i) {
        seed = seed * 1664525u + 1013904223u;
        s += "ACGT"[(seed >> 25) & 3];
    }
    return s;
}

std::string unit(int length, unsigned seed, bool reverse) {
    std::string s = dna(8, seed) + "TCTTCAGCGTTCCCGAGA" + dna(6, seed + 1) +
                    dna(7, seed + 2) + std::string(20, 'T') + dna(length, seed + 3);
    return reverse ? concat_hmm::detail::revcomp(s) : s;
}

read_streaming::sequence read(const std::string& id, const std::string& seq) {
    std::string qual;
    for (size_t i = 0; i < seq.size(); ++i) qual += char('!' + i % 60);
    return {id, "original comment", seq, qual, true};
}

void tests() {
    concat_hmm::model_options opt;
    opt.gate = false;
    opt.foldback_split = false; // deliberately no sequence evidence for the F/R junction
    opt.tau_junction = 1.0f;    // exercises coordinate uncertainty together with posterior abstention
    auto model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto parent = read("unresolved_F_R", unit(355, 12, false) + dna(83, 17) + unit(517, 22, true));
    concat_hmm::Scratch scratch;
    concat_hmm::Result result;
    concat_hmm::segment(*model, parent.seq.data(), (int)parent.seq.size(), scratch, result);
    require(result.k == 2 && result.cut_flags.size() == 1 &&
            (result.cut_flags[0] & concat_hmm::BOUNDARY_UNRESOLVED), "Fixture did not yield an unresolved junction");
    require(result.flags & concat_hmm::RES_ABSTAIN, "Fixture must also exercise legacy-abstain precedence");

    // This is a router test: the callback stands for downstream barcode assignment. No whitelist or
    // variable extraction is needed to prove that blocked children never reach that callback.
    for (bool legacy : {false, true}) {
        ReadLayout layout;
        layout.concat_model = model;
        layout.concat_abstain_legacy = legacy;
        concat_hmm_counters counters;
        concat_layout_info info;
        info.ctr = &counters;
        int assigned = 0, saved = 0;
        const bool handled = SigString::concat_hmm_route(parent, layout, info,
            [&](SigString&, const read_streaming::sequence&) { ++assigned; }, false,
            [&](const read_streaming::sequence& original) {
                ++saved;
                require(original.id == parent.id && original.comment == parent.comment &&
                        original.seq == parent.seq && original.qual == parent.qual && original.is_fastq,
                        "Unresolved parent was altered");
            });
        require(handled && saved == 1 && assigned == 0, "Unresolved junction leaked into barcode assignment/fallback");
        require(counters.boundary_unresolved_parents == 1 && counters.boundary_unresolved_children == 2,
                "Unresolved routing counters disagree");
        require(counters.legacy_abstain == 0, "Legacy route bypassed unresolved-boundary protection");
    }

    // Actual HMM inference on three native reverse capture units with unequal payloads and an unknown
    // scaffold. A located one-sided edge must reach normal downstream processing, not be quarantined.
    opt.tau_junction = 0.9f;
    model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto repeated = read("reverse_repeat", unit(275, 100, true) + dna(253, 101) +
                                               unit(410, 200, true) + dna(37, 201) + unit(620, 300, true));
    concat_hmm::segment(*model, repeated.seq.data(), (int)repeated.seq.size(), scratch, result);
    require(result.k == 3, "Reverse-repeat fixture did not produce three captures");
    for (uint8_t f : result.cut_flags)
        require((f & concat_hmm::BOUNDARY_LEFT_BOUNDED) && !(f & concat_hmm::BOUNDARY_UNRESOLVED),
                "Reverse-repeat safe edge is unresolved");
    ReadLayout layout;
    layout.concat_model = model;
    concat_hmm_counters counters;
    concat_layout_info info;
    info.ctr = &counters;
    int assigned = 0, saved = 0;
    const bool handled = SigString::concat_hmm_route(repeated, layout, info,
        [&](SigString&, const read_streaming::sequence& child) {
            ++assigned;
            require(child.seq.size() == child.qual.size(), "Child sequence/quality lengths differ");
            const size_t pos = repeated.seq.find(child.seq);
            require(pos != std::string::npos && repeated.qual.substr(pos, child.qual.size()) == child.qual,
                    "Child sequence/quality no longer maps to parent");
        }, false, [&](const read_streaming::sequence&) { ++saved; });
    require(handled && assigned == 3 && saved == 0, "One-sided repeat was withheld or sent to legacy");
    require(counters.boundary_unresolved_parents == 0, "Resolved repeat counted as unresolved");
    require(!result.boundary_unresolved_seen, "Unresolved parent status leaked into next read");

    // The only unknown junction is after the retained cap. Dropping the per-cut array must not
    // erase the obligation to preserve the raw parent or permit a legacy re-extraction.
    opt.k_cap = 2;
    model = std::make_shared<concat_hmm::Model>(concat_hmm::Model::build(curio(), nullptr, opt));
    const auto capped = read("unknown_after_cap", unit(275, 410, false) + unit(410, 510, false) +
                                                 unit(355, 610, false) + unit(517, 710, true));
    concat_hmm::segment(*model, capped.seq.data(), (int)capped.seq.size(), scratch, result);
    require(result.k == 2 && (result.flags & concat_hmm::RES_TOO_MANY), "Cap fixture did not exceed cap");
    require(result.boundary_unresolved_seen, "Cap erased unresolved-parent status");
    for (uint8_t f : result.cut_flags)
        require(!(f & concat_hmm::BOUNDARY_UNRESOLVED), "Cap fixture uncertainty must be beyond retained cuts");
    for (bool legacy : {false, true}) {
        ReadLayout capped_layout;
        capped_layout.concat_model = model;
        capped_layout.concat_abstain_legacy = legacy;
        concat_hmm_counters capped_counters;
        concat_layout_info capped_info;
        capped_info.ctr = &capped_counters;
        int written = 0, preserved = 0;
        const bool cap_handled = SigString::concat_hmm_route(capped, capped_layout, capped_info,
            [&](SigString&, const read_streaming::sequence&) { ++written; }, false,
            [&](const read_streaming::sequence& original) {
                ++preserved;
                require(original.seq == capped.seq && original.qual == capped.qual, "Capped parent changed");
            });
        require(cap_handled && preserved == 1 && written == 0, "Cap allowed unresolved parent to escape preservation");
        require(capped_counters.legacy_too_many == 0, "Cap reached legacy despite an unresolved boundary");
    }
}
} // namespace

int main() {
    try { tests(); std::cout << "concat_hmm router regression tests passed\n"; return 0; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
