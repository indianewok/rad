#include "rad/sigstring.hpp"

#include <functional>
#include <iostream>
#include <optional>
#include <string>

namespace {

constexpr const char* kCollisionBarcode = "AAACATCG";

barcode_entry barcode_value(const int64_seq& barcode) {
    barcode_entry entry;
    entry.barcode = barcode;
    entry.filtered = false;
    return entry;
}

struct Fixture {
    ReadLayout layout;
    whitelist::wl_entry* whitelist_entry = nullptr;

    Fixture() {
        auto [it, inserted] =
            layout.wl_map.lists.try_emplace("filter-precedence-fixture");
        (void)inserted;
        whitelist_entry = &it->second;
        layout.wl_map.maps.emplace(
            "barcode_1", std::ref(*whitelist_entry));
    }
};

std::optional<int64_seq> correct(
    Fixture& fixture,
    const std::string& observed,
    const std::string& direction = "forward") {
    seq_element element(
        "barcode_1",
        "barcode",
        0,
        {1, static_cast<int>(observed.size())},
        "variable",
        1,
        direction,
        true,
        true,
        observed,
        std::string(observed.size(), 'I'),
        observed);

    read_streaming::sequence read{
        "filter-precedence-test",
        "",
        observed,
        std::string(observed.size(), 'I'),
        true};

    return barcode_correction::correct_barcode(
        element, fixture.layout, read, false, 1, "offensive");
}

bool expect_barcode(
    const std::string& label,
    const std::optional<int64_seq>& actual,
    const std::string& expected) {
    if (actual.has_value() &&
        actual->bits_to_sequence() == expected) {
        return true;
    }

    std::cerr << "FAIL: " << label << " expected " << expected
              << ", got "
              << (actual.has_value()
                      ? actual->bits_to_sequence()
                      : std::string("<rejected>"))
              << '\n';
    return false;
}

bool expect_rejected(
    const std::string& label,
    const std::optional<int64_seq>& actual) {
    if (!actual.has_value()) {
        return true;
    }

    std::cerr << "FAIL: " << label << " expected rejection, got "
              << actual->bits_to_sequence() << '\n';
    return false;
}

bool exact_true_collision_is_accepted() {
    Fixture fixture;
    const int64_seq collision(kCollisionBarcode);
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(collision);
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        collision, barcode_value(collision));

    return expect_barcode(
        "exact trusted identity overrides static collision",
        correct(fixture, kCollisionBarcode),
        kCollisionBarcode);
}

bool reverse_complement_identity_is_accepted() {
    Fixture fixture;
    const int64_seq canonical(kCollisionBarcode);
    const std::string observed = seq_utils::revcomp(kCollisionBarcode);
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(
        int64_seq(observed));
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        canonical, barcode_value(canonical));

    return expect_barcode(
        "reverse-complement trusted identity overrides static collision",
        correct(fixture, observed),
        kCollisionBarcode);
}

bool identity_wins_over_same_key_alias() {
    Fixture fixture;
    const int64_seq observed(kCollisionBarcode);
    const int64_seq alias_target("AAACATCA");
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(observed);
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        observed, barcode_value(observed));
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        observed, barcode_value(alias_target));

    return expect_barcode(
        "identity wins when the same key also has a mutation alias",
        correct(fixture, kCollisionBarcode),
        kCollisionBarcode);
}

bool reverse_identity_wins_over_direct_alias() {
    Fixture fixture;
    const int64_seq canonical(kCollisionBarcode);
    const std::string observed = seq_utils::revcomp(kCollisionBarcode);
    const int64_seq observed_bits(observed);
    const int64_seq alias_target("CGATGTTA");
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(observed_bits);
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        canonical, barcode_value(canonical));
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        observed_bits, barcode_value(alias_target));

    return expect_barcode(
        "reverse-complement identity wins over a direct mutation alias",
        correct(fixture, observed),
        kCollisionBarcode);
}

bool non_whitelist_static_hit_is_rejected() {
    Fixture fixture;
    const int64_seq collision(kCollisionBarcode);
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(collision);

    return expect_rejected(
        "non-whitelist static k-mer remains filtered",
        correct(fixture, kCollisionBarcode));
}

bool mutation_alias_is_rejected() {
    Fixture fixture;
    const int64_seq observed(kCollisionBarcode);
    const int64_seq canonical("AAACATCA");
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(observed);
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        canonical, barcode_value(canonical));
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        observed, barcode_value(canonical));

    return expect_rejected(
        "mutation association cannot override static collision",
        correct(fixture, kCollisionBarcode));
}

bool global_only_collision_is_rejected() {
    Fixture fixture;
    const int64_seq collision(kCollisionBarcode);
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(collision);
    fixture.whitelist_entry->global_bcs.insert_bc_entry(
        collision, barcode_value(collision));

    return expect_rejected(
        "global-only identity cannot override static collision",
        correct(fixture, kCollisionBarcode));
}

bool low_complexity_identity_is_rejected() {
    Fixture fixture;
    const int64_seq low_complexity("AAAAAAAA");
    fixture.whitelist_entry->filter_bcs.insert_bc_entry(low_complexity);
    fixture.whitelist_entry->true_bcs.insert_bc_entry(
        low_complexity, barcode_value(low_complexity));

    return expect_rejected(
        "low-complexity sequence remains filtered",
        correct(fixture, "AAAAAAAA"));
}

}  // namespace

int main() {
    bool ok = true;
    ok = exact_true_collision_is_accepted() && ok;
    ok = reverse_complement_identity_is_accepted() && ok;
    ok = identity_wins_over_same_key_alias() && ok;
    ok = reverse_identity_wins_over_direct_alias() && ok;
    ok = non_whitelist_static_hit_is_rejected() && ok;
    ok = mutation_alias_is_rejected() && ok;
    ok = global_only_collision_is_rejected() && ok;
    ok = low_complexity_identity_is_rejected() && ok;
    return ok ? 0 : 1;
}
