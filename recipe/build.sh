#!/usr/bin/env bash
set -euo pipefail

# macOS: keep CMake from preferring system frameworks / app bundles over $PREFIX.
if [[ "$(uname)" == "Darwin" ]]; then
  EXTRA_CMAKE_ARGS=(-DCMAKE_FIND_FRAMEWORK=NEVER -DCMAKE_FIND_APPBUNDLE=NEVER)
else
  EXTRA_CMAKE_ARGS=()
fi

# The CMakeLists install() rules put the real binary in libexec/rad/rad, a
# self-locating wrapper at bin/rad, and the bundled layouts/whitelists under
# share/rad/resources -- so `cmake --install` is all that's needed (no manual
# `install` of the bare binary, which is what shipped a resource-less v0.6.0).
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_COMPILER="${CC}" \
  "${EXTRA_CMAKE_ARGS[@]}"

cmake --build build -j "${CPU_COUNT}"
cmake --install build
