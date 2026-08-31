# BroCOLI 2.0 + preBroCOLI 2.0
#
# Under an activated conda environment, `conda activate` exports CXX, CXXFLAGS
# and LDFLAGS. We *append* to them rather than replacing them, otherwise the
# conda sysroot/rpath settings are lost (and, with `CXXFLAGS ?=`, so is
# -std=c++17). Nothing here needs editing for a conda build.
#
#   make            build both binaries (prebrocoli is skipped if edlib is
#                   missing - see `make check-env`)
#   make brocoli    quantification only, needs htslib
#   make prebrocoli barcode extraction only, needs zlib + edlib + WFA2-lib
#   make test       unit tests (no htslib, no WFA2 needed)

CXX      ?= g++
OPT      ?= -O3
STD      := -std=c++17
WARN     := -Wall -Wextra -Wno-unused-parameter

ifeq ($(NATIVE),1)
  OPT += -march=native
endif

ALL_CXXFLAGS := $(STD) $(OPT) $(WARN) $(CXXFLAGS)
ALL_LDFLAGS  := $(LDFLAGS)

# ------------------------------------------------------------------ htslib --
HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS   := $(shell pkg-config --libs   htslib 2>/dev/null)
ifeq ($(strip $(HTS_LIBS)),)
  HTS_LIBS := -lhts
endif

# ------------------------------------------------- WFA2-lib (C++ bindings) --
# The conda package installs headers under $CONDA_PREFIX/include/wfa2lib and
# the C++ binding as libwfa2cpp.
comma := ,
WFA_PREFIX ?= $(CONDA_PREFIX)
WFA_CFLAGS ?= $(if $(WFA_PREFIX),-I$(WFA_PREFIX)/include/wfa2lib -I$(WFA_PREFIX)/include,)
WFA_LIBS   ?= $(if $(WFA_PREFIX),-L$(WFA_PREFIX)/lib -Wl$(comma)-rpath$(comma)$(WFA_PREFIX)/lib,) -lwfa2cpp -lwfa2

# ------------------------------------------------------------------- edlib --
# Two vendored files (MIT). Copy edlib.h and edlib.cpp from the BroCOLI 1.x
# tree into third_party/edlib/ - there is no reliable conda package for the
# C++ library itself.
EDLIB_DIR ?= third_party/edlib
EDLIB_SRC := $(EDLIB_DIR)/edlib.cpp

BIN      := brocoli
PRE_BIN  := prebrocoli
SRC      := src/main.cpp
PRE_SRC  := pre/main.cpp
HDRS     := $(wildcard src/*.hpp)
PRE_HDRS := $(wildcard pre/*.hpp)

.PHONY: all clean test debug check-env pre-if-possible

all: $(BIN) pre-if-possible

$(BIN): $(SRC) $(HDRS)
	$(CXX) $(ALL_CXXFLAGS) -Isrc $(HTS_CFLAGS) $(SRC) -o $@ $(ALL_LDFLAGS) $(HTS_LIBS) -lpthread
	@echo "[BroCOLI] built ./$(BIN)"

$(PRE_BIN): $(PRE_SRC) $(PRE_HDRS) $(EDLIB_SRC)
	$(CXX) $(ALL_CXXFLAGS) -Isrc -Ipre -I$(EDLIB_DIR) $(WFA_CFLAGS) \
	    $(PRE_SRC) $(EDLIB_SRC) -o $@ $(ALL_LDFLAGS) $(WFA_LIBS) -lz -lpthread
	@echo "[preBroCOLI] built ./$(PRE_BIN)"

# Build prebrocoli only when its dependencies are present, so a fresh clone
# without edlib still produces a working brocoli.
pre-if-possible:
	@if [ -f "$(EDLIB_SRC)" ]; then \
	    $(MAKE) --no-print-directory $(PRE_BIN); \
	else \
	    echo "[preBroCOLI] skipped: $(EDLIB_SRC) not found"; \
	    echo "             copy edlib.h and edlib.cpp into $(EDLIB_DIR)/ and re-run make"; \
	fi

debug: ALL_CXXFLAGS := $(STD) -O0 -g -fsanitize=address,undefined $(WARN) $(CXXFLAGS)
debug: clean $(BIN)

# Neither test binary needs htslib, edlib or WFA2.
test:
	$(CXX) $(STD) -O2 -g -Isrc tests/test_pipeline.cpp -o /tmp/brocoli_test
	/tmp/brocoli_test
	$(CXX) $(STD) -O2 -g -Ipre tests/test_seqio.cpp -o /tmp/prebrocoli_seqio_test -lz
	/tmp/prebrocoli_seqio_test

check-env:
	@echo "CXX        = $(CXX)"
	@echo "CXXFLAGS   = $(ALL_CXXFLAGS)"
	@echo "htslib cf  = $(HTS_CFLAGS)"
	@echo "htslib lib = $(HTS_LIBS)"
	@pkg-config --modversion htslib 2>/dev/null && echo "htslib found" || echo "htslib NOT found by pkg-config"
	@echo "WFA cflags = $(WFA_CFLAGS)"
	@echo "WFA libs   = $(WFA_LIBS)"
	@test -f "$(WFA_PREFIX)/include/wfa2lib/bindings/cpp/WFAligner.hpp" \
	    && echo "WFAligner.hpp found" || echo "WFAligner.hpp NOT found"
	@test -f "$(EDLIB_SRC)" && echo "edlib found in $(EDLIB_DIR)" \
	    || echo "edlib NOT found (copy edlib.h/edlib.cpp into $(EDLIB_DIR))"

clean:
	rm -f $(BIN) $(PRE_BIN) /tmp/brocoli_test /tmp/prebrocoli_seqio_test
