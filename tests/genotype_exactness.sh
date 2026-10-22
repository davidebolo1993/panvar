#!/usr/bin/env bash
# Run the exactness contracts from the configured release build.
set -euo pipefail
build="${1:-build}"
shift "$(( $# > 0 ? 1 : 0 ))"
exec ctest --test-dir "$build" --output-on-failure --no-tests=error -L '^exactness$' "$@"
