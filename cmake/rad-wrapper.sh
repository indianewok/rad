#!/usr/bin/env bash
# Thin launcher installed as <prefix>/bin/rad.
#
# RAD locates its bundled read layouts + whitelists via find_resource_root(),
# which checks $RAD_RESOURCES first. We set it from THIS script's own location
# so it keeps working after conda relocates the environment (the prefix is not
# hard-coded). Installed layout:
#   <prefix>/bin/rad                 <- this wrapper
#   <prefix>/libexec/rad/rad         <- the real binary
#   <prefix>/share/rad/resources/    <- read_layout/ + wl/
set -euo pipefail

# Resolve the directory this script lives in (follow symlinks).
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
BIN_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
PREFIX="$(cd -P "$BIN_DIR/.." && pwd)"

# Only set RAD_RESOURCES if the caller hasn't overridden it.
if [ -z "${RAD_RESOURCES:-}" ] && [ -d "$PREFIX/share/rad/resources" ]; then
  export RAD_RESOURCES="$PREFIX/share/rad"
fi

exec "$PREFIX/libexec/rad/rad" "$@"
