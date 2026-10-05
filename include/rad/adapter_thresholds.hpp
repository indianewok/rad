#pragma once

#include <algorithm>
#include <cstddef>

namespace adapter_thresholds {

constexpr std::size_t kFallbackErrorNumerator = 3;
constexpr std::size_t kFallbackErrorDenominator = 10;

/**
 * Return the shared fallback edit-distance limit for an uncalibrated adapter.
 *
 * Keep this integer-only so layout preparation and read processing cannot
 * drift because of different rounding rules or duplicated constants.
 */
constexpr int fallback_max_edit_distance(std::size_t adapter_length) {
    const auto rounded_up =
        (adapter_length * kFallbackErrorNumerator +
         kFallbackErrorDenominator - 1) /
        kFallbackErrorDenominator;
    return std::max(1, static_cast<int>(rounded_up));
}

static_assert(fallback_max_edit_distance(22) == 7,
              "22 nt adapters must allow seven fallback edits");
static_assert(fallback_max_edit_distance(23) == 7,
              "23 nt adapters must allow seven fallback edits");

/**
 * Accepted share of chance hits at misalign_lower.
 *
 * The misalignment pass accrues, per static element, the whole-read best edit
 * distance of that element in reads that hold an exact copy of every static
 * element of the other direction (the null). misalign_lower is the largest
 * edit distance at which the cumulative share of null reads with a chance hit
 * at or below it is <= this value. The value is an accepted chance-hit rate,
 * not a property of the reads, so it is not calibrated. Misalignment_Setup
 * takes it as a constructor parameter.
 */
constexpr double kMisalignLowerChanceShare = 0.05;

/**
 * Minimum number of chance hits per static element for a calibrated
 * misalignment threshold. With fewer, update_read_layout writes the fallback
 * (fallback_max_edit_distance). The calibration sample grows until every
 * element reaches this count (see below).
 */
constexpr std::size_t kMinMisalignmentObservations = 100;

/**
 * Calibration sample of the misalignment pass.
 *
 * The pass reads the first kCalibrationBlockReads reads (rad demux; rad prep
 * -n sets the size of this first block). When an element has fewer than
 * kMinMisalignmentObservations chance hits, it reads further blocks of
 * kCalibrationBlockReads reads from the same stream. It stops when every
 * element has enough chance hits, at the end of the input, at
 * kCalibrationMaxReads reads, or early when the projected number of reads
 * (reads so far x kMinMisalignmentObservations / chance hits of the slowest
 * element) exceeds kCalibrationMaxReads. Whole blocks keep the sample a fixed
 * prefix of the file, independent of the thread count.
 *
 * kCalibrationMaxReads is a safety bound on the time spent in calibration,
 * not a property of the reads.
 */
constexpr std::size_t kCalibrationBlockReads = 50000;
constexpr std::size_t kCalibrationMaxReads = 1000000;

}  // namespace adapter_thresholds
