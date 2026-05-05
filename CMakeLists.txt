cmake_minimum_required(VERSION 3.18)
project(rad VERSION 0.1.0 LANGUAGES C CXX)
include(FetchContent)
# --------------------------------------------------------------------------
#language standard
# --------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
function(set_required_build_settings_for_GCC8)
    # Always link with libstdc++fs.a when using GCC 8.
    # Note: This command makes sure that this option comes pretty late on the cmdline.
    link_libraries( "$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:-lstdc++fs>" )
endfunction()
set_required_build_settings_for_GCC8()

# --------------------------------------------------------------------------
# OpenMP (auto-detect libomp on macOS via Homebrew or MacPorts)
# --------------------------------------------------------------------------
if(APPLE)
    # Check common Homebrew locations (Apple Silicon + Intel)
    foreach(_omp_prefix
        "/opt/homebrew/opt/libomp"
        "/usr/local/opt/libomp"
        "/opt/local"
    )
        if(EXISTS "${_omp_prefix}" AND NOT OpenMP_FOUND)
            list(APPEND CMAKE_PREFIX_PATH "${_omp_prefix}")
        endif()
    endforeach()
endif()
find_package(OpenMP REQUIRED)

# --------------------------------------------------------------------------
# other dependencies (the policy part is so that it doesn't complain)
# --------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)     
endif()

find_program(PIGZ_EXE NAMES pigz)
if (PIGZ_EXE)
  message(STATUS "Using pigz at: ${PIGZ_EXE}")
else()
  message(STATUS "pigz not found; will fall back to zlib at runtime")
endif()

find_package(Boost REQUIRED COMPONENTS filesystem iostreams)

# --------------------------------------------------------------------------
# third party libraries
# --------------------------------------------------------------------------

# Striped Smith‑Waterman (stand‑alone C code)
add_library(ssw STATIC
    include/ssw/ssw.c
    include/ssw/ssw_cpp.cpp
)

#header file for edlib
add_library(edlib STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/edlib/src/edlib.cpp
)

target_include_directories(ssw PUBLIC include/ssw)
target_include_directories(edlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/edlib/include)
# --------------------------------------------------------------------------
# header‑only rad library (collects all include paths & link deps)
# --------------------------------------------------------------------------
file(GLOB_RECURSE RAD_HEADERS CONFIGURE_DEPENDS rad/*.hpp)

add_library(rad_lib INTERFACE)
target_sources(rad_lib INTERFACE ${RAD_HEADERS})
target_compile_definitions(rad_lib INTERFACE HAVE_PIGZ=1 PIGZ_PATH="${PIGZ_EXE}" RAD_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

target_include_directories(rad_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include                          
    include/kseq
    include/kdepp
    include/edlib
    include/sciplot
    include/parallel_hashmap
    include/csv-parser/single_include
)

target_link_libraries(rad_lib INTERFACE
    OpenMP::OpenMP_CXX
    Boost::filesystem
    Boost::iostreams
    Boost::boost
    ZLIB::ZLIB
    edlib
    ssw
    dl
)

# --------------------------------------------------------------------------
# Executables
# --------------------------------------------------------------------------
add_executable(rad src/main.cpp)
add_executable(rad_config src/config_tools.cpp)
add_executable(generate_whitelist src/whitelist_generator.cpp)
add_executable(discover_layout src/discover_layout.cpp)

# pull in every dependency through rad_lib
target_link_libraries(rad PRIVATE rad_lib)
target_link_libraries(rad_config PRIVATE rad_lib)
target_link_libraries(generate_whitelist PRIVATE rad_lib)
target_link_libraries(discover_layout PRIVATE rad_lib)

# --------------------------------------------------------------------------
message(STATUS "Configure-time compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
