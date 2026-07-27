#include "rad/rad_headers.h"

#include <iostream>
#include <string>

namespace {
bool require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << message << '\n';
        return false;
    }
    return true;
}
}

int main() {
    ParsedPosition empty;
    if (!require(empty.ref_id.empty(), "default ref_id is not empty") ||
        !require(!empty.is_start, "default is_start is not false") ||
        !require(empty.offset == 0, "default offset is not zero") ||
        !require(empty.add_flags.empty(), "default flags are not empty") ||
        !require(empty.to_string().empty(),
                 "default position did not serialize empty")) {
        return 1;
    }

    const ParsedPosition parsed_empty =
        ParsedPosition::from_string("", false);
    if (!require(parsed_empty.to_string().empty(),
                 "empty position string did not remain empty")) {
        return 1;
    }

    const ParsedPosition legacy_empty =
        ParsedPosition::from_string("|start+1919251551", false);
    if (!require(legacy_empty.to_string().empty(),
                 "legacy empty-reference position was not normalized")) {
        return 1;
    }

    const std::string valid_text = "linker_2|start-18|chained_start";
    const ParsedPosition valid =
        ParsedPosition::from_string(valid_text, false);
    if (!require(valid.to_string() == valid_text,
                 "valid position did not round-trip")) {
        return 1;
    }

    return 0;
}
