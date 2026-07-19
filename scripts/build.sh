#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &amp;&amp; pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/build"

BUILD_TYPE=${1:-Release}
JOBS=${JOBS:-$(nproc 2&gt;/dev/null || echo 4)}

echo "=== NR Link Simulator Build Script ==="
echo "Project dir: $PROJECT_DIR"
echo "Build dir:   $BUILD_DIR"
echo "Build type:  $BUILD_TYPE"
echo "Jobs:        $JOBS"
echo ""

if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
fi

cd "$BUILD_DIR"

echo "--- Configuring with CMake ---"
cmake "$PROJECT_DIR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DNR_BUILD_TESTS=ON \
    -DNR_BUILD_EXAMPLES=ON \
    -DNR_USE_OPENMP=ON

echo ""
echo "--- Building ---"
cmake --build . -j "$JOBS"

echo ""
echo "--- Build complete! ---"
echo ""
echo "To run the example simulation:"
echo "  cd $BUILD_DIR/examples"
echo "  ./pdsch_bler_simulation --help"
echo ""
echo "To run tests:"
echo "  cd $BUILD_DIR"
echo "  ctest --output-on-failure"
