if(NOT DEFINED RAD_EXE OR NOT EXISTS "${RAD_EXE}")
    message(FATAL_ERROR "RAD_EXE must point to a built rad executable")
endif()
if(NOT DEFINED RAD_TEST_DIR)
    message(FATAL_ERROR "RAD_TEST_DIR is required")
endif()

file(MAKE_DIRECTORY "${RAD_TEST_DIR}")

set(layout_path "${RAD_TEST_DIR}/splitseq_r1_layout.csv")
set(fastq_path "${RAD_TEST_DIR}/splitseq_r1.fastq")
set(max_output_base "${RAD_TEST_DIR}/max")
set(verbose_output_base "${RAD_TEST_DIR}/verbose")
file(REMOVE
    "${layout_path}"
    "${fastq_path}"
    "${max_output_base}_layout.csv"
    "${max_output_base}_position_map.csv"
    "${verbose_output_base}_layout.csv"
    "${verbose_output_base}_position_map.csv"
)

file(WRITE "${layout_path}" [=[Read Layout:R1,,,,,
id,seq,expected_length,type,class,whitelist
read,,,variable,read,
barcode_1,,8,variable,barcode,
linker_1,CCACAGTCTCAAGCACGTGGAT,,static,linker,
barcode_2,,8,variable,barcode,
linker_2,AGTCGTACGCCGATGCGAAACATCGGCCAC,,static,linker,
barcode_3,,8,variable,barcode,
umi,,10,variable,umi,
]=])

set(read_sequence
    "GATTACAAAACCGGTTCCACAGTCTCAAGCACGTGGATCCTTAAGGAGTCGTACGCCGATGCGAAACATCGGCCACACGTACGTTTGATTACAG")
string(LENGTH "${read_sequence}" read_length)
string(REPEAT "I" ${read_length} read_quality)
foreach(read_number RANGE 1 10)
    file(APPEND "${fastq_path}"
        "@splitseq_${read_number}\n${read_sequence}\n+\n${read_quality}\n")
endforeach()

function(run_prep verbosity output_base result_var output_var)
    execute_process(
        COMMAND "${RAD_EXE}" prep
            --layout "${layout_path}"
            --read-layout
            --position-map
            --fastq "${fastq_path}"
            --output "${output_base}"
            --max-reads 10
            --threads 1
            "${verbosity}"
        RESULT_VARIABLE command_result
        OUTPUT_VARIABLE command_stdout
        ERROR_VARIABLE command_stderr
    )
    set(${result_var} "${command_result}" PARENT_SCOPE)
    set(${output_var} "${command_stdout}\n${command_stderr}" PARENT_SCOPE)
endfunction()

run_prep("--max-verbose" "${max_output_base}" max_result max_output)
if(NOT max_result EQUAL 0)
    message(FATAL_ERROR "max-verbose prep failed (${max_result}):\n${max_output}")
endif()

run_prep("--verbose" "${verbose_output_base}" verbose_result verbose_output)
if(NOT verbose_result EQUAL 0)
    message(FATAL_ERROR "verbose prep failed (${verbose_result}):\n${verbose_output}")
endif()

set(required_max_verbose_fragments
    "[read_layout] Complete read layout (max verbose):"
    "[read_layout] Element barcode_3:"
)

foreach(fragment IN LISTS required_max_verbose_fragments)
    string(FIND "${max_output}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "max-verbose output is missing expected fragment:\n${fragment}\n\n"
            "Full output:\n${max_output}")
    endif()
endforeach()

set(bc3_trace_start_marker
    "[position_map] Processing right element: barcode_3")
string(FIND "${max_output}" "${bc3_trace_start_marker}" bc3_trace_start)
if(bc3_trace_start EQUAL -1)
    message(FATAL_ERROR "max-verbose output is missing the BC3 trace block")
endif()

string(SUBSTRING "${max_output}" ${bc3_trace_start} -1 bc3_trace_tail)
string(FIND "${bc3_trace_tail}"
    "[position_map] Processing right element: umi"
    bc3_trace_end)
if(bc3_trace_end EQUAL -1)
    message(FATAL_ERROR "could not locate the end of the BC3 trace block")
endif()
string(SUBSTRING "${bc3_trace_tail}" 0 ${bc3_trace_end} bc3_trace_block)

set(required_bc3_trace_fragments
    "[position_map] Processing right element: barcode_3"
    "Previous: linker_2 (order=6, type=static, class=linker, direction=forward)"
    "Next: umi (order=8, type=variable, class=umi, direction=forward)"
    "Next-next: seq_stop (order=9, type=static, class=stop, direction=forward)"
    "Secondary selection: previous static anchor linker_2"
    "Primary selection: promoting static anchor linker_2 over terminal-linked chain"
    "Primary Start: linker_2|stop+1"
    "Secondary Start: umi|start-8|poly_chained_start"
    "Primary Stop: linker_2|stop+8"
    "Secondary Stop: umi|start-1|poly_chained_stop"
)

foreach(fragment IN LISTS required_bc3_trace_fragments)
    string(FIND "${bc3_trace_block}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "BC3 trace block is missing expected fragment:\n${fragment}\n\n"
            "BC3 trace block:\n${bc3_trace_block}")
    endif()
endforeach()

set(umi_trace_start_marker
    "[position_map] Processing right element: umi")
string(FIND "${max_output}" "${umi_trace_start_marker}" umi_trace_start)
if(umi_trace_start EQUAL -1)
    message(FATAL_ERROR "max-verbose output is missing the UMI trace block")
endif()

string(SUBSTRING "${max_output}" ${umi_trace_start} -1 umi_trace_tail)
string(FIND "${umi_trace_tail}"
    "[position_map] Parse Right positions for umi"
    umi_assignment_start)
if(umi_assignment_start EQUAL -1)
    message(FATAL_ERROR "could not locate the UMI position assignment")
endif()
string(SUBSTRING "${umi_trace_tail}" 0 ${umi_assignment_start} umi_context_block)
string(SUBSTRING "${umi_trace_tail}" ${umi_assignment_start} -1 umi_assignment_tail)
string(FIND "${umi_assignment_tail}"
    "[position_map] Processing direction"
    umi_trace_end)
if(umi_trace_end EQUAL -1)
    string(LENGTH "${umi_assignment_tail}" umi_trace_end)
endif()
string(SUBSTRING "${umi_assignment_tail}" 0 ${umi_trace_end} umi_assignment_block)
set(umi_trace_block "${umi_context_block}${umi_assignment_block}")

set(required_umi_trace_fragments
    "[position_map] Processing right element: umi"
    "Previous: barcode_3 (order=7, type=variable, class=barcode, direction=forward)"
    "Next: seq_stop (order=9, type=static, class=stop, direction=forward)"
    "Previous-previous: linker_2 (order=6, type=static, class=linker, direction=forward)"
    "Secondary selection: previous-previous static anchor linker_2"
    "Primary selection: promoting static anchor linker_2 over terminal-linked chain"
    "Primary Start: linker_2|stop+9"
    "Secondary Start: seq_stop|start-10"
    "Primary Stop: linker_2|stop+18"
    "Secondary Stop: seq_stop|start-1"
)

foreach(fragment IN LISTS required_umi_trace_fragments)
    string(FIND "${umi_trace_block}" "${fragment}" fragment_position)
    if(fragment_position EQUAL -1)
        message(FATAL_ERROR
            "UMI trace block is missing expected fragment:\n${fragment}\n\n"
            "UMI trace block:\n${umi_trace_block}")
    endif()
endforeach()

string(REGEX MATCHALL
    "\\[position_map\\] Processing right element: barcode_3"
    bc3_trace_matches
    "${max_output}")
list(LENGTH bc3_trace_matches bc3_trace_count)
if(NOT bc3_trace_count EQUAL 1)
    message(FATAL_ERROR
        "expected one BC3 decision trace, found ${bc3_trace_count}")
endif()

string(FIND "${verbose_output}"
    "[position_map] Processing right element: barcode_3"
    unexpected_verbose_trace)
if(NOT unexpected_verbose_trace EQUAL -1)
    message(FATAL_ERROR
        "ordinary --verbose unexpectedly emitted max-verbose position tracing")
endif()

foreach(suffix IN ITEMS "_layout.csv" "_position_map.csv")
    set(max_file "${max_output_base}${suffix}")
    set(verbose_file "${verbose_output_base}${suffix}")
    if(NOT EXISTS "${max_file}" OR NOT EXISTS "${verbose_file}")
        message(FATAL_ERROR "prep did not create both ${suffix} outputs")
    endif()

    file(SHA256 "${max_file}" max_hash)
    file(SHA256 "${verbose_file}" verbose_hash)
    if(NOT max_hash STREQUAL verbose_hash)
        message(FATAL_ERROR
            "verbosity changed generated ${suffix} content")
    endif()
endforeach()
