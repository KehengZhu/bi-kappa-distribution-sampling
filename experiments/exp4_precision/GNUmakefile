# Experiment 4 — finite-precision / low-kappa audit.
#
# The whole point of the experiment is that the answer may depend on the standard
# library, so the same source is built twice: once against libc++ (Apple clang) and
# once against libstdc++ (Homebrew GCC).  Both binaries are kept side by side.

CXXFLAGS  = -Wall -Wextra -std=c++11 -O2

CXX_LIBCXX     = clang++
CXX_LIBSTDCXX  = $(shell command -v g++-15 || command -v g++-14 || command -v g++-13)

SRC       = exp4_probe.cpp
HDR       = ../../cpp/bi_kappa_distribution.H
RAWDIR    = raw

TARGETS   = exp4_probe_libcxx.exe
ifneq ($(CXX_LIBSTDCXX),)
TARGETS  += exp4_probe_libstdcxx.exe
endif

N ?= 1000000

all: $(TARGETS)

exp4_probe_libcxx.exe: $(SRC) $(HDR)
	$(CXX_LIBCXX) $(CXXFLAGS) -stdlib=libc++ -o $@ $(SRC)

exp4_probe_libstdcxx.exe: $(SRC) $(HDR)
	$(CXX_LIBSTDCXX) $(CXXFLAGS) -o $@ $(SRC)

# Regenerate every raw result from scratch.  One command, no hidden state.
run: all
	mkdir -p $(RAWDIR)
	@for b in $(TARGETS); do \
	  tag=$${b%.exe}; tag=$${tag#exp4_probe_}; \
	  echo "  $$tag: variates"; ./$$b variates $(N) > $(RAWDIR)/variates_$$tag.jsonl; \
	  echo "  $$tag: released"; ./$$b released $(N) > $(RAWDIR)/released_$$tag.jsonl; \
	  echo "  $$tag: capped";   ./$$b capped   $(N) > $(RAWDIR)/capped_$$tag.jsonl; \
	done
	@echo "  logdump (distributional validation of the log-domain path)"
	@for k in 0.55 0.75 2.0 5.0; do \
	  for s in 4001 4002 4003; do \
	    ./exp4_probe_libcxx.exe logdump $$k $$s 200000 $(RAWDIR)/logr_k$${k}_s$${s}_double.bin double; \
	  done; \
	done
	cd $(RAWDIR) && shasum -a 256 *.jsonl logr_*.bin > checksums.sha256

# Verify a regenerated raw/ against the committed checksums.  The .jsonl counters are
# themselves tracked; the logr_*.bin bulk draws are not, so this is the only way to
# confirm a regenerated set matches the one the committed results were computed from.
verify:
	cd $(RAWDIR) && shasum -a 256 --check checksums.sha256 | grep -v ': OK$$' || \
	  echo "all files match checksums.sha256"

clean:
	rm -f $(TARGETS)

distclean: clean
	rm -rf $(RAWDIR)

.PHONY: all run verify clean distclean
