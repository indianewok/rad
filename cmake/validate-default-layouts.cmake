if(NOT DEFINED RAD_SOURCE_DIR)
    message(FATAL_ERROR "RAD_SOURCE_DIR is required")
endif()

set(_layout_dir "${RAD_SOURCE_DIR}/resources/read_layout")
set(_mapping_header "${RAD_SOURCE_DIR}/include/rad/misc_utils.hpp")

function(rad_assert_layout_row contract)
    string(REPLACE "|" ";" fields "${contract}")
    list(GET fields 0 layout_file)
    list(GET fields 1 expected_row)
    set(path "${_layout_dir}/${layout_file}")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Missing bundled read layout: ${path}")
    endif()

    file(STRINGS "${path}" rows)
    list(FIND rows "${expected_row}" row_index)
    if(row_index EQUAL -1)
        message(FATAL_ERROR
            "Read-layout contract failed for ${layout_file}: missing row '${expected_row}'")
    endif()
endfunction()

function(rad_assert_layout_lacks layout_file forbidden_text)
    set(path "${_layout_dir}/${layout_file}")
    file(READ "${path}" contents)
    string(FIND "${contents}" "${forbidden_text}" text_index)
    if(NOT text_index EQUAL -1)
        message(FATAL_ERROR
            "Read-layout contract failed for ${layout_file}: unexpected '${forbidden_text}'")
    endif()
endfunction()

function(rad_assert_registered contract)
    string(REPLACE "|" ";" fields "${contract}")
    list(GET fields 0 key)
    list(GET fields 1 relative_path)
    string(REPLACE " " "" expected "{\"${key}\",\"${relative_path}\"}")
    string(FIND "${mapping_text}" "${expected}" mapping_index)
    if(mapping_index EQUAL -1)
        message(FATAL_ERROR
            "Default mapping contract failed: ${key} must resolve to ${relative_path}")
    endif()

    if(NOT EXISTS "${RAD_SOURCE_DIR}/${relative_path}")
        message(FATAL_ERROR
            "Default mapping contract failed: ${relative_path} does not exist")
    endif()
endfunction()

# Chemistry contracts for every bundled non-bulk layout. The scTagger fixture
# intentionally uses a 10 bp UMI even though its barcode catalog is 10x v3.
set(layout_contracts
    "five_prime_read_layout.csv|barcode,,16,variable,barcode,10x_5v1"
    "five_prime_read_layout.csv|umi,,10,variable,umi,"
    "five_prime_read_layout.csv|tso,TTTCTTATATGGG,,static,tso,"
    "five_prime_read_layout.csv|poly_a,,,static,poly_a,"
    "three_prime_read_layout.csv|barcode,,16,variable,barcode,10x_3v3"
    "three_prime_read_layout.csv|umi,,12,variable,umi,"
    "three_prime_read_layout.csv|poly_t,,,static,poly_t,"
    "sctagger_sim_read_layout.csv|barcode,,16,variable,barcode,10x_3v3"
    "sctagger_sim_read_layout.csv|umi,,10,variable,umi,"
    "sctagger_sim_read_layout.csv|poly_a,,,static,poly_a,"
    "splitseq_read_layout.csv|barcode_1,,8,variable,barcode,splitseq_bc2"
    "splitseq_read_layout.csv|barcode_2,,8,variable,barcode,splitseq_bc1"
    "splitseq_read_layout.csv|barcode_3,,8,variable,barcode,splitseq_bc1"
    "splitseq_read_layout.csv|umi,,10,variable,umi,"
    "visium_three_prime_read_layout.csv|barcode,,16,variable,barcode,10x_Vis_V1"
    "visium_three_prime_read_layout.csv|umi,,12,variable,umi,"
    "visium_three_prime_read_layout.csv|poly_t,,,static,poly_t,"
    "visium_hd_read_layout.csv|umi,,9,variable,umi,,"
    "visium_hd_read_layout.csv|barcode_1,,15-16,variable,barcode,visium_hd_bc1,joint_barcode"
    "visium_hd_read_layout.csv|barcode_2,,14-15,variable,barcode,visium_hd_bc2,joint_barcode"
)
foreach(contract IN LISTS layout_contracts)
    rad_assert_layout_row("${contract}")
endforeach()

# Bulk rapid-barcoding is deliberately barcode/UMI-free.
rad_assert_layout_lacks(nanopore_bulk_rapid_bc_read_layout.csv ",variable,barcode,")
rad_assert_layout_lacks(nanopore_bulk_rapid_bc_read_layout.csv ",variable,umi,")

file(READ "${_mapping_header}" mapping_text)
string(REGEX REPLACE "[ \t\r\n]" "" mapping_text "${mapping_text}")

set(mapping_contracts
    "five_prime|resources/read_layout/five_prime_read_layout.csv"
    "sctagger|resources/read_layout/sctagger_sim_read_layout.csv"
    "three_prime|resources/read_layout/three_prime_read_layout.csv"
    "splitseq|resources/read_layout/splitseq_read_layout.csv"
    "visium|resources/read_layout/visium_three_prime_read_layout.csv"
    "visium_hd|resources/read_layout/visium_hd_read_layout.csv"
    "nanopore_rapid_bc|resources/read_layout/nanopore_bulk_rapid_bc_read_layout.csv"
    "10x_3v1|resources/wl/737K-august-2016_bitlist.csv.gz"
    "10x_3v2|resources/wl/737K-august-2016_bitlist.csv.gz"
    "10x_3v3|resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz"
    "10x_3v3.1|resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz"
    "10x_3HTv3.1|resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz"
    "10x_3v4|resources/wl/3M-3pgex-may-2023.txt_bitlist.csv.gz"
    "10x_5v1|resources/wl/737K-august-2016_bitlist.csv.gz"
    "10x_5v2|resources/wl/737K-august-2016_bitlist.csv.gz"
    "10x_5HTv2|resources/wl/737K-august-2016_bitlist.csv.gz"
    "10x_5v3|resources/wl/3M-5pgex-jan-2023.txt_bitlist.csv.gz"
    "10x_Vis_V1|resources/wl/visium-v1_v2_bitlist.csv.gz"
    "10x_Vis_V2|resources/wl/visium-v1_v2_bitlist.csv.gz"
    "10x_Vis_V3|resources/wl/visium-v3_v4_bitlist.csv.gz"
    "10x_Vis_V4|resources/wl/visium-v3_v4_bitlist.csv.gz"
    "10x_Vis_V5|resources/wl/visium-v5_bitlist.csv.gz"
    "splitseq_bc1|resources/wl/splitseq_bc1_bitlist.csv.gz"
    "splitseq_bc2|resources/wl/splitseq_bc2_bitlist.csv.gz"
    "visium_hd_bc1|resources/wl/visium_hd_bc1.csv.gz"
    "visium_hd_bc2|resources/wl/visium_hd_bc2.csv.gz"
)
foreach(contract IN LISTS mapping_contracts)
    rad_assert_registered("${contract}")
endforeach()

message(STATUS "Bundled read-layout and kit contracts are valid")
