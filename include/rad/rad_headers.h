#pragma once

// ————————————————————————————————————————————————————————————————
// C++ Standard Library in alphabetical order
// ————————————————————————————————————————————————————————————————
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ————————————————————————————————————————————————————————————————
// Parallelism
// ————————————————————————————————————————————————————————————————
#include <atomic>
#include <omp.h>

// ————————————————————————————————————————————————————————————————
// Boost
// ————————————————————————————————————————————————————————————————
#include <boost/bimap.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/multi_index_container.hpp>

// ————————————————————————————————————————————————————————————————
// External libraries
// ————————————————————————————————————————————————————————————————
#include "csv-parser/single_include/csv.hpp"
#include "parallel_hashmap/phmap.h"
#include "edlib/include/edlib.h"
#include "sciplot/sciplot.hpp"
#include "gzstream/gzstream.h"
#include "ssw/ssw_cpp.h"
#include "kdepp/kde.h"
#include "kseq/kseq.h"
#include "zlib.h"

// ————————————————————————————————————————————————————————————————
//  rad headers
// ————————————————————————————————————————————————————————————————
#include "misc_utils.hpp"
#include "io_streaming.hpp"
#include "barcode_correction.hpp"
#include "read_layout.hpp"
#include "sigstring.hpp"
