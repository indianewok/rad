if(NOT DEFINED RAD_EXE OR NOT EXISTS "${RAD_EXE}")
    message(FATAL_ERROR "RAD_EXE must point to a built rad executable")
endif()
if(NOT DEFINED RAD_TEST_DIR)
    message(FATAL_ERROR "RAD_TEST_DIR is required")
endif()

find_program(GZIP_EXE NAMES gzip)
if(NOT GZIP_EXE)
    message(FATAL_ERROR "gzip is required to inspect RAD's compressed output")
endif()

file(MAKE_DIRECTORY "${RAD_TEST_DIR}")

set(layout_path "${RAD_TEST_DIR}/splitseq_r1_layout.csv")
set(fastq_path "${RAD_TEST_DIR}/splitseq_r1.fastq")
set(whitelist_path "${RAD_TEST_DIR}/component_barcodes.txt")
set(output_name "terminal")
set(output_base "${RAD_TEST_DIR}/${output_name}")

file(WRITE "${layout_path}" [=[Read Layout,,,,,
id,seq,expected_length,type,class,whitelist
read,,,variable,read,
barcode_1,,8,variable,barcode,
linker_1,CCACAGTCTCAAGCACGTGGAT,,static,linker,
barcode_2,,8,variable,barcode,
linker_2,AGTCGTACGCCGATGCGAAACATCGGCCAC,,static,linker,
barcode_3,,8,variable,barcode,
umi,,10,variable,umi,
]=])

file(WRITE "${whitelist_path}" [=[AACCGGTT
CCTTAAGG
ACGTACGT
AAACATCG
]=])

# Truth, in 1-based inclusive coordinates:
# read 1..8, BC1 9..16, linker1 17..38, BC2 39..46,
# linker2 47..76, BC3 77..84, UMI 85..94.
set(read_sequence
    "GATTACAAAACCGGTTCCACAGTCTCAAGCACGTGGATCCTTAAGGAGTCGTACGCCGATGCGAAACATCGGCCACACGTACGTTTGATTACAG")
string(LENGTH "${read_sequence}" read_length)
string(REPEAT "I" ${read_length} read_quality)
set(read_sequence_bc2_collision
    "GATTACAAAACCGGTTCCACAGTCTCAAGCACGTGGATAAACATCGAGTCGTACGCCGATGCGAAACATCGGCCACACGTACGTTTGATTACAG")
set(read_sequence_bc3_collision
    "GATTACAAAACCGGTTCCACAGTCTCAAGCACGTGGATCCTTAAGGAGTCGTACGCCGATGCGAAACATCGGCCACAAACATCGTTGATTACAG")
file(WRITE "${fastq_path}"
    "@splitseq_1\n${read_sequence}\n+\n${read_quality}\n"
    "@splitseq_collision_bc2\n${read_sequence_bc2_collision}\n+\n${read_quality}\n"
    "@splitseq_collision_bc3\n${read_sequence_bc3_collision}\n+\n${read_quality}\n")

file(REMOVE
    "${output_base}_layout.csv"
    "${output_base}_position_map.csv"
    "${output_base}_dbg.fq.gz"
    "${output_base}_dbg.csv.gz"
    "${output_base}_dbg.sig.gz")

execute_process(
    COMMAND "${RAD_EXE}" prep
        --layout "${layout_path}"
        --read-layout
        --position-map
        --fastq "${fastq_path}"
        --output "${output_base}"
        --max-reads 1
        --threads 1
    RESULT_VARIABLE prep_result
    OUTPUT_VARIABLE prep_stdout
    ERROR_VARIABLE prep_stderr
    TIMEOUT 60
)
if(NOT prep_result EQUAL 0)
    message(FATAL_ERROR
        "SPLIT-seq prep failed (${prep_result}):\n"
        "${prep_stdout}\n${prep_stderr}")
endif()

# Exercise import of genuinely empty optional positions. Generated maps use
# the primary position as a harmless duplicate fallback, but user-authored and
# older maps may leave either secondary endpoint empty.
set(position_map "${output_base}_position_map.csv")
file(READ "${position_map}" position_map_content)
string(REPLACE
    "read,seq_start|stop+1,barcode_1|start-1,seq_start|stop+1,linker_1|start-9|right_skipped,"
    "read,seq_start|stop+1,barcode_1|start-1,,linker_1|start-9|right_skipped,"
    position_map_content "${position_map_content}")
string(REPLACE
    "rc_read,rc_barcode_1|stop+1,rc_seq_stop|start-1,rc_linker_1|stop+9|left_skipped,rc_seq_stop|start-1,"
    "rc_read,rc_barcode_1|stop+1,rc_seq_stop|start-1,rc_linker_1|stop+9|left_skipped,,"
    position_map_content "${position_map_content}")
file(WRITE "${position_map}" "${position_map_content}")

foreach(expected_empty_position IN ITEMS
    "read,seq_start|stop+1,barcode_1|start-1,,linker_1|start-9|right_skipped,"
    "rc_read,rc_barcode_1|stop+1,rc_seq_stop|start-1,rc_linker_1|stop+9|left_skipped,,")
    string(FIND "${position_map_content}" "${expected_empty_position}"
        empty_position_index)
    if(empty_position_index EQUAL -1)
        message(FATAL_ERROR
            "failed to construct empty-position import fixture:\n"
            "${position_map_content}")
    endif()
endforeach()

execute_process(
    COMMAND "${RAD_EXE}" demux
        --layout "${layout_path}"
        --fastq "${fastq_path}"
        --global-whitelist "${whitelist_path}"
        --custom-whitelist "${whitelist_path}"
        --output "${output_name}"
        --dir "${RAD_TEST_DIR}"
        --threads 1
        --chunk-size 1
        --max-reads 3
        --whitelist-mutation 0
        --generated-mutation 0
        --min-read-length 0
        --no-umi-rc
        --write-dbg
        --max-verbose
    RESULT_VARIABLE command_result
    OUTPUT_VARIABLE command_stdout
    ERROR_VARIABLE command_stderr
    TIMEOUT 60
)
if(NOT command_result EQUAL 0)
    message(FATAL_ERROR
        "SPLIT-seq terminal-boundary demux failed (${command_result}):\n"
        "${command_stdout}\n${command_stderr}")
endif()

string(FIND "${command_stdout}" "[pos_map] Importing position map..."
    import_message_position)
if(import_message_position EQUAL -1)
    message(FATAL_ERROR
        "demux did not import the prepared map with empty secondary fields:\n"
        "${command_stdout}\n${command_stderr}")
endif()

set(output_fastq "${RAD_TEST_DIR}/${output_name}_dbg.fq.gz")
set(debug_csv "${RAD_TEST_DIR}/${output_name}_dbg.csv.gz")
set(debug_sig "${RAD_TEST_DIR}/${output_name}_dbg.sig.gz")
foreach(output_path IN ITEMS "${output_fastq}" "${debug_csv}" "${debug_sig}")
    if(NOT EXISTS "${output_path}")
        message(FATAL_ERROR "expected output was not created: ${output_path}")
    endif()
endforeach()

function(read_gzip path result_var)
    execute_process(
        COMMAND "${GZIP_EXE}" -dc "${path}"
        RESULT_VARIABLE gzip_result
        OUTPUT_VARIABLE gzip_output
        ERROR_VARIABLE gzip_error
    )
    if(NOT gzip_result EQUAL 0)
        message(FATAL_ERROR "could not read ${path}: ${gzip_error}")
    endif()
    set(${result_var} "${gzip_output}" PARENT_SCOPE)
endfunction()

read_gzip("${output_fastq}" fastq_output)
read_gzip("${debug_csv}" debug_csv_output)
read_gzip("${debug_sig}" debug_sig_output)

set(required_fastq_fragments
    "CB:Z:AACCGGTT-CCTTAAGG-ACGTACGT"
    "CB:Z:AACCGGTT-AAACATCG-ACGTACGT"
    "CB:Z:AACCGGTT-CCTTAAGG-AAACATCG"
    "UB:Z:TTGATTACAG"
    "\nGATTACAA\n"
    "\nIIIIIIII\n"
)
foreach(fragment IN LISTS required_fastq_fragments)
    string(FIND "${fastq_output}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "FASTQ output is missing the expected terminal extraction "
            "'${fragment}':\n${fastq_output}")
    endif()
endforeach()

set(required_debug_csv_fragments
    "splitseq_1,barcode_3,forward,ACGTACGT"
    "splitseq_collision_bc2,barcode_2,forward,AAACATCG"
    "splitseq_collision_bc3,barcode_3,forward,AAACATCG"
    "splitseq_1,read,forward,GATTACAA"
    "splitseq_1,umi,forward,TTGATTACAG"
)
foreach(fragment IN LISTS required_debug_csv_fragments)
    string(FIND "${debug_csv_output}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "debug CSV is missing '${fragment}':\n${debug_csv_output}")
    endif()
endforeach()

set(required_signature_fragments
    "seq_start:0:0:0"
    "read:0:1:8"
    "barcode_3:0:77:84"
    "umi:0:85:94"
    "seq_stop:0:95:95"
)
foreach(fragment IN LISTS required_signature_fragments)
    string(FIND "${debug_sig_output}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "debug signature is missing '${fragment}':\n${debug_sig_output}")
    endif()
endforeach()
