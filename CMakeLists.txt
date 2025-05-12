cmake_minimum_required(VERSION 3.18)
project(cdev_rad VERSION 0.1.0 LANGUAGES C CXX)

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
# OpenMP (Homebrew libomp on macOS)
# --------------------------------------------------------------------------
if(APPLE)
# hardcoded--change it
    #list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp")
endif()
find_package(OpenMP REQUIRED)

# --------------------------------------------------------------------------
# other dependencies (the policy part is so that it doesn't complain)
# --------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)     
endif()
find_package(Boost REQUIRED COMPONENTS filesystem system)

# --------------------------------------------------------------------------
# third part libraries
# --------------------------------------------------------------------------
add_subdirectory(include/edlib EXCLUDE_FROM_ALL)

# Make edlib examples not interfere with build
option(EDLIB_BUILD_EXAMPLES "Build example apps" OFF)
option(EDLIB_BUILD_TESTS "Build tests" OFF)

# Striped Smith‑Waterman (stand‑alone C code)
add_library(ssw STATIC
    include/ssw/ssw.c
    include/ssw/ssw_cpp.cpp
)
target_include_directories(ssw PUBLIC include/ssw)

# --------------------------------------------------------------------------
# header‑only rad library (collects all include paths & link deps)
# --------------------------------------------------------------------------
file(GLOB_RECURSE RAD_HEADERS CONFIGURE_DEPENDS rad/*.hpp)

add_library(rad_lib INTERFACE)
target_sources(rad_lib INTERFACE ${RAD_HEADERS})

target_include_directories(rad_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include                          
    include/kseq
    include/csv-parser
    include/progressbar
)

target_link_libraries(rad_lib INTERFACE
    Boost::boost
    Boost::filesystem
    Boost::system
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX
    edlib
    ssw
    dl
)

# --------------------------------------------------------------------------
# Executables
# --------------------------------------------------------------------------
add_executable(rad src/main.cpp)
add_executable(rad_config src/config_tools.cpp)

# pull in every dependency through rad_lib
target_link_libraries(rad PRIVATE rad_lib)
target_link_libraries(rad_config PRIVATE rad_lib)

# --------------------------------------------------------------------------
message(STATUS "Configure-time compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
