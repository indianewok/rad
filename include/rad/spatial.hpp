#pragma once

#include <string>

struct spatial_reformat_options {
    std::string input, output, sample, map_path, bc1_path, bc2_path, metadata_path;
    std::string unit = "bin";
    int bin_size = 2;
    bool tissue_only = false;
    bool strict = false;
};

void usage_reformat_spatial();
int reformat_spatial(const spatial_reformat_options& options);
