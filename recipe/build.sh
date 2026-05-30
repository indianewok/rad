#!/usr/bin/env bash
set -euo pipefail

# macOS: keep CMake from preferring system frameworks / app bundles over $PREFIX,
# and force the deployment target that std::filesystem needs (>= macOS 10.15).
# bioconda pins osx-64 to 10.13; conda_build_config.yaml raises it for osx-64.
if [[ "$(uname)" == "Darwin" ]]; then
  EXTRA_CMAKE_ARGS=(
    -DCMAKE_FIND_FRAMEWORK=NEVER
    -DCMAKE_FIND_APPBUNDLE=NEVER
    -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-10.15}"
  )
else
  EXTRA_CMAKE_ARGS=()
fi

# The CMakeLists install() rules put the real binary in libexec/rad/rad, a
# self-locating wrapper at bin/rad, and the bundled layouts/whitelists under
# share/rad/resources -- so `cmake --install` is all that's needed (no manual
# `install` of the bare binary, which is what shipped a resource-less v0.6.0).
# ${CMAKE_ARGS} carries the conda toolchain settings (sysroot, deployment target,
# etc.); it must be forwarded to the configure step.
cmake -S . -B build ${CMAKE_ARGS:-} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_C_COMPILER="${CC}" \
  "${EXTRA_CMAKE_ARGS[@]}"

cmake --build build -j "${CPU_COUNT}"
cmake --install build
