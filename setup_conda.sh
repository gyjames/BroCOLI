#!/usr/bin/env bash
# One-shot setup: create the conda environment and build BroCOLI + preBroCOLI.
#
# Deliberately does NOT depend on Makefile/CMakeLists.txt being present or
# intact: everything needed to build and test is inlined below with plain
# g++ + pkg-config calls, so a corrupted or missing Makefile cannot block
# environment setup. The Makefile is still there for day-to-day development.
set -euo pipefail

ENV_NAME="${1:-brocoli}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v conda >/dev/null 2>&1 || { echo "conda not found on PATH"; exit 1; }

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[BroCOLI] environment '$ENV_NAME' already exists, updating"
    conda env update -n "$ENV_NAME" -f "$HERE/environment.yml" --prune
else
    echo "[BroCOLI] creating environment '$ENV_NAME'"
    conda env create -n "$ENV_NAME" -f "$HERE/environment.yml"
fi

# `conda activate` needs the shell hook when running non-interactively.
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

CXX="${CXX:-g++}"

echo "[BroCOLI] toolchain:"
echo "  CXX     = $(command -v "$CXX")"
"$CXX" --version | head -1

if ! pkg-config --exists htslib 2>/dev/null; then
    echo "[BroCOLI] ERROR: pkg-config cannot find htslib in this environment." >&2
    echo "           Check that environment.yml's htslib package installed correctly:" >&2
    echo "           conda list -n $ENV_NAME htslib" >&2
    exit 1
fi
HTS_VERSION="$(pkg-config --modversion htslib)"
HTS_CFLAGS="$(pkg-config --cflags htslib)"
HTS_LIBS="$(pkg-config --libs htslib)"
echo "  htslib  = $HTS_VERSION"
echo "  cflags  = $HTS_CFLAGS"
echo "  libs    = $HTS_LIBS"

echo "[BroCOLI] source files:"
if [ ! -f "$HERE/src/main.cpp" ]; then
    echo "[BroCOLI] ERROR: $HERE/src/main.cpp not found - re-check your upload/checkout." >&2
    exit 1
fi
NHPP=$(find "$HERE/src" -name '*.hpp' | wc -l)
echo "  found $NHPP header(s) under src/"
if [ "$NHPP" -lt 11 ]; then
    echo "[BroCOLI] ERROR: expected 11 headers under src/, found $NHPP." >&2
    echo "           The source tree looks incomplete - re-check your upload/checkout." >&2
    exit 1
fi

echo "[BroCOLI] running unit tests (no htslib needed for this part)..."
"$CXX" -std=c++17 -O2 -g -I"$HERE/src" "$HERE/tests/test_pipeline.cpp" -o /tmp/brocoli_test
/tmp/brocoli_test
"$CXX" -std=c++17 -O2 -g -I"$HERE/pre" "$HERE/tests/test_seqio.cpp" -o /tmp/prebrocoli_seqio_test -lz
/tmp/prebrocoli_seqio_test

echo "[BroCOLI] building ./brocoli ..."
"$CXX" -std=c++17 -O3 -Wall -Wextra -Wno-unused-parameter \
    -I"$HERE/src" $HTS_CFLAGS \
    "$HERE/src/main.cpp" -o "$HERE/brocoli" \
    $HTS_LIBS -lpthread

echo "[BroCOLI] build OK: $HERE/brocoli"
# NB: never pipe the binary straight into `head` here. It flushes the banner
# and the usage text separately, so `head` can close the pipe in between; the
# binary then dies of SIGPIPE and, with `set -o pipefail`, takes this whole
# script down with it right after a successful build.
BROCOLI_HELP="$("$HERE/brocoli" -h 2>&1 || true)"
printf '%s\n' "$BROCOLI_HELP" | head -7

# --------------------------------------------------------------- preBroCOLI --
# Needs zlib (always there), edlib (vendored) and the WFA2-lib C++ bindings.
# The include layout of the conda package has changed between releases, so
# instead of hard-coding one guess we probe candidate flag sets with a real
# two-line compile and use the first that works.
echo
echo "[preBroCOLI] checking dependencies ..."

EDLIB_DIR="$HERE/third_party/edlib"
if [ ! -f "$EDLIB_DIR/edlib.cpp" ] || [ ! -f "$EDLIB_DIR/edlib.h" ]; then
    echo "[preBroCOLI] edlib not found in third_party/edlib/."
    echo "             Copy edlib.h and edlib.cpp there (they are in the BroCOLI 1.x"
    echo "             tree) and re-run this script. BroCOLI itself is already built."
    exit 0
fi
echo "  edlib   = $EDLIB_DIR"

PROBE=$(mktemp /tmp/wfa_probe_XXXXXX.cpp)
cat > "$PROBE" <<'PROBEEOF'
#include "bindings/cpp/WFAligner.hpp"
#include <string>
int main() {
    wfa::WFAlignerGapAffine a(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    a.alignEnd2End(std::string("ACGT"), std::string("ACGT"));
    return a.getAlignmentScore() == 0 ? 0 : 0;
}
PROBEEOF

WFA_CFLAGS=""
WFA_LIBS=""
found=0
P="${CONDA_PREFIX:-}"
# -L/-rpath are added explicitly: plain `g++` inside a conda env does not
# search $CONDA_PREFIX/lib unless LDFLAGS says so.
LDP="${P:+-L$P/lib -Wl,-rpath,$P/lib}"
for cand_c in "-I$P/include/wfa2lib -I$P/include" "-I$P/include/wfa2lib" \
              "$(pkg-config --cflags libwfa2 2>/dev/null || true)" "-I$P/include/wfa2"; do
    for cand_l in "$LDP -lwfa2cpp -lwfa2" "$LDP -lwfa2cpp" "$LDP -lwfacpp -lwfa"; do
        if "$CXX" -std=c++17 -O0 $cand_c "$PROBE" -o /tmp/wfa_probe $cand_l >/dev/null 2>&1; then
            WFA_CFLAGS="$cand_c"; WFA_LIBS="$cand_l"; found=1; break 2
        fi
    done
done
rm -f "$PROBE" /tmp/wfa_probe

if [ "$found" -ne 1 ]; then
    echo "[preBroCOLI] ERROR: could not compile against the WFA2-lib C++ bindings." >&2
    echo "             Looked for bindings/cpp/WFAligner.hpp under \$CONDA_PREFIX/include." >&2
    echo "             conda list -n $ENV_NAME wfa2-lib" >&2
    echo "             ls \$CONDA_PREFIX/lib | grep wfa" >&2
    exit 1
fi
echo "  WFA2    = $WFA_CFLAGS  $WFA_LIBS"

echo "[preBroCOLI] building ./prebrocoli ..."
"$CXX" -std=c++17 -O3 -Wall -Wextra -Wno-unused-parameter \
    -I"$HERE/src" -I"$HERE/pre" -I"$EDLIB_DIR" $WFA_CFLAGS \
    "$HERE/pre/main.cpp" "$EDLIB_DIR/edlib.cpp" -o "$HERE/prebrocoli" \
    $WFA_LIBS -lz -lpthread

echo "[preBroCOLI] build OK: $HERE/prebrocoli"
#"$HERE/prebrocoli" -h | head -5

cat <<MSG

[BroCOLI] done.

  conda activate $ENV_NAME
  $HERE/brocoli -h
  $HERE/prebrocoli -h

MSG
