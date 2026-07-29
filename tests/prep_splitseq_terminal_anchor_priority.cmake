if(NOT DEFINED RAD_EXE OR NOT EXISTS "${RAD_EXE}")
    message(FATAL_ERROR "RAD_EXE must point to a built rad executable")
endif()
if(NOT DEFINED RAD_SOURCE_DIR)
    message(FATAL_ERROR "RAD_SOURCE_DIR is required")
endif()
if(NOT DEFINED RAD_TEST_DIR)
    message(FATAL_ERROR "RAD_TEST_DIR is required")
endif()

set(layout_path
    "${RAD_SOURCE_DIR}/resources/read_layout/splitseq_read_layout.csv")
if(NOT EXISTS "${layout_path}")
    message(FATAL_ERROR "bundled SPLIT-seq layout was not found")
endif()

file(MAKE_DIRECTORY "${RAD_TEST_DIR}")
set(empty_fastq "${RAD_TEST_DIR}/empty.fastq")
set(output_base "${RAD_TEST_DIR}/splitseq")
file(WRITE "${empty_fastq}" "")
file(REMOVE
    "${output_base}_layout.csv"
    "${output_base}_position_map.csv")

execute_process(
    COMMAND "${RAD_EXE}" prep
        --layout "${layout_path}"
        --position-map
        --fastq "${empty_fastq}"
        --output "${output_base}"
        --max-reads 1
        --threads 1
    RESULT_VARIABLE command_result
    OUTPUT_VARIABLE command_stdout
    ERROR_VARIABLE command_stderr
)
if(NOT command_result EQUAL 0)
    message(FATAL_ERROR
        "bundled SPLIT-seq prep failed (${command_result}):\n"
        "${command_stdout}\n${command_stderr}")
endif()

set(position_map "${output_base}_position_map.csv")
if(NOT EXISTS "${position_map}")
    message(FATAL_ERROR "prep did not create ${position_map}")
endif()
file(READ "${position_map}" position_map_content)

set(expected_rows
    "read,seq_start|stop+1,barcode_1|start-1,seq_start|stop+1,linker_1|start-9|right_skipped,"
    "barcode_3,linker_2|stop+1,linker_2|stop+8,umi|start-8|poly_chained_start,umi|start-1|poly_chained_stop,"
    "umi,linker_2|stop+9,linker_2|stop+18,seq_stop|start-10,seq_stop|start-1,"
    "rc_umi,rc_linker_2|start-18,rc_linker_2|start-9,rc_seq_start|stop+1,rc_seq_start|stop+10,"
    "rc_barcode_3,rc_linker_2|start-8,rc_linker_2|start-1,rc_umi|stop+1|var_chained_start,rc_umi|stop+8|var_chained_stop,"
    "rc_read,rc_barcode_1|stop+1,rc_seq_stop|start-1,rc_linker_1|stop+9|left_skipped,rc_seq_stop|start-1,"
)

foreach(expected_row IN LISTS expected_rows)
    string(FIND "${position_map_content}" "${expected_row}" row_position)
    if(row_position EQUAL -1)
        message(FATAL_ERROR
            "SPLIT-seq position map is missing expected row:\n"
            "${expected_row}\n\nFull map:\n${position_map_content}")
    endif()
endforeach()

string(REGEX MATCH ",\\|(start|stop)" orphan_position
    "${position_map_content}")
if(NOT orphan_position STREQUAL "")
    message(FATAL_ERROR
        "position map contains an anchor with an empty reference ID: "
        "${orphan_position}")
endif()
