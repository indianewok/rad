if(NOT DEFINED RAD_EXE OR NOT EXISTS "${RAD_EXE}")
    message(FATAL_ERROR "RAD_EXE must point to a built rad executable")
endif()
if(NOT DEFINED RAD_TEST_DIR)
    message(FATAL_ERROR "RAD_TEST_DIR is required")
endif()

file(MAKE_DIRECTORY "${RAD_TEST_DIR}")

set(short_layout "${RAD_TEST_DIR}/five_prime_tso_12.csv")
set(long_layout "${RAD_TEST_DIR}/five_prime_tso_13.csv")
set(fastq_path "${RAD_TEST_DIR}/five_prime.fastq")

set(layout_prefix [=[Read Layout:R1,,,,,
id,seq,expected_length,type,class,whitelist
forw_primer,CTACACGACGCTCTTCCGATCT,,static,primer,
barcode,,16,variable,barcode,
umi,,10,variable,umi,
]=])
set(layout_suffix [=[read,,,variable,read,
]=])

file(WRITE "${short_layout}"
    "${layout_prefix}tso,TTTCTTATATGG,,static,tso,\n${layout_suffix}")
file(WRITE "${long_layout}"
    "${layout_prefix}tso,TTTCTTATATGGG,,static,tso,\n${layout_suffix}")

set(read_sequence
    "CTACACGACGCTCTTCCGATCTAACCGGTTAACCGGTTTGCATGCATGTTTCTTATATGGGGATTACAGATTACA")
string(LENGTH "${read_sequence}" read_length)
string(REPEAT "I" ${read_length} read_quality)
file(WRITE "${fastq_path}" "")
foreach(read_number RANGE 1 10)
    file(APPEND "${fastq_path}"
        "@five_prime_${read_number}\n${read_sequence}\n+\n${read_quality}\n")
endforeach()

function(run_prep layout_path output_base result_var output_var)
    execute_process(
        COMMAND "${RAD_EXE}" prep
            --layout "${layout_path}"
            --read-layout
            --position-map
            --fastq "${fastq_path}"
            --output "${output_base}"
            --max-reads 10
            --threads 1
            --max-verbose
        RESULT_VARIABLE command_result
        OUTPUT_VARIABLE command_stdout
        ERROR_VARIABLE command_stderr
    )
    set(${result_var} "${command_result}" PARENT_SCOPE)
    set(${output_var} "${command_stdout}\n${command_stderr}" PARENT_SCOPE)
endfunction()

run_prep("${short_layout}" "${RAD_TEST_DIR}/short" short_result short_output)
if(NOT short_result EQUAL 0)
    message(FATAL_ERROR "12-base TSO prep failed (${short_result}):\n${short_output}")
endif()

run_prep("${long_layout}" "${RAD_TEST_DIR}/long" long_result long_output)
if(NOT long_result EQUAL 0)
    message(FATAL_ERROR "13-base TSO prep failed (${long_result}):\n${long_output}")
endif()

string(FIND "${short_output}"
    "Secondary selection: previous static fallback primer"
    short_fallback_position)
if(short_fallback_position EQUAL -1)
    message(FATAL_ERROR
        "12-base TSO was not rejected as a secondary anchor:\n${short_output}")
endif()

string(FIND "${short_output}"
    "Secondary selection: next-next static anchor tso"
    short_anchor_position)
if(NOT short_anchor_position EQUAL -1)
    message(FATAL_ERROR
        "12-base TSO was unexpectedly accepted as a secondary anchor:\n${short_output}")
endif()

string(FIND "${long_output}"
    "Secondary selection: next-next static anchor tso"
    long_anchor_position)
if(long_anchor_position EQUAL -1)
    message(FATAL_ERROR
        "13-base TSO was not accepted as a secondary anchor:\n${long_output}")
endif()

set(long_position_map "${RAD_TEST_DIR}/long_position_map.csv")
file(READ "${long_position_map}" long_position_map_content)
string(FIND "${long_position_map_content}"
    "umi,barcode|stop+1|var_chained_start,barcode|stop+10|var_chained_stop,"
    ordinary_chain_position)
if(ordinary_chain_position EQUAL -1)
    message(FATAL_ERROR
        "ordinary nonterminal UMI chain was unexpectedly promoted:\n"
        "${long_position_map_content}")
endif()
