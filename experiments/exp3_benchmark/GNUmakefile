# Experiment 3 -- performance benchmark.
#
# Correctness is checked BEFORE timing, and `make run` enforces that order:
# no timing number is meaningful until every compared method is verified to
# sample the intended target distribution.

CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2
TARGET   = exp3_bench.exe
SRC      = exp3_bench.cpp
HDR      = ../../cpp/bi_kappa_distribution.H
RAWDIR   = raw

# Validation: the kappa ladder includes kappa <= 1, where the Zenitani (2025)
# recommended envelope index n = kappa/2 is inapplicable.  That is recorded, not
# skipped.
VAL_KAPPA   = 0.75 1.0 1.5 2.0 5.0 10.0
VAL_SEEDS   = 3001 3002 3003
VAL_N       = 200000

# Timing: same kappa ladder plus 50 to expose kappa dependence.
TIME_KAPPA  = 0.75 1.0 1.5 2.0 5.0 10.0 50.0
TIME_N      = 1000000
TIME_REPEAT = 10
TIME_SEED   = 3101

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Regenerate everything from scratch, then record checksums.  The .bin validation
# dumps are not committed; checksums.sha256 and the .jsonl summaries are.  `make run`
# overwrites the committed raw/validate.jsonl and raw/timing.jsonl; `git checkout raw/`
# restores them.
run: $(TARGET)
	mkdir -p $(RAWDIR)
	@echo "--- phase 1: distributional validation (must pass before timing) ---"
	@rm -f $(RAWDIR)/validate.jsonl
	@for k in $(VAL_KAPPA); do \
	  for s in $(VAL_SEEDS); do \
	    ./$(TARGET) validate $(RAWDIR)/val_k$${k}_s$${s}.bin $$k 1.0 $$s $(VAL_N) \
	      >> $(RAWDIR)/validate.jsonl; \
	  done; \
	done
	@echo "--- phase 2: timing ---"
	@rm -f $(RAWDIR)/timing.jsonl
	@for k in $(TIME_KAPPA); do \
	  for m in gamma_ratio_spherical scale_mixture_normals pareto_rejection; do \
	    ./$(TARGET) time $$m iso $$k 1.0 $(TIME_SEED) $(TIME_N) $(TIME_REPEAT) \
	      >> $(RAWDIR)/timing.jsonl; \
	  done; \
	  for var in aniso rotated capped20; do \
	    ./$(TARGET) time gamma_ratio_spherical $$var $$k 1.0 $(TIME_SEED) \
	      $(TIME_N) $(TIME_REPEAT) >> $(RAWDIR)/timing.jsonl; \
	  done; \
	done
	cd $(RAWDIR) && shasum -a 256 *.jsonl val_*.bin > checksums_rerun.txt

# Check raw/ against the committed raw/checksums.sha256, which `make run` never overwrites
# (it writes the checksums of the new files to raw/checksums_rerun.txt, not committed).
# Fails on any mismatch and on any missing file.  A rerun is not expected to pass: the
# timings differ from run to run, and see README.md for the validation dumps.
verify:
	@cd $(RAWDIR) && if shasum -a 256 --check --quiet checksums.sha256; then \
	  echo "all files match raw/checksums.sha256"; \
	else \
	  echo "raw/ does not match raw/checksums.sha256 (lines above name each file)" >&2; \
	  exit 1; \
	fi

clean:
	rm -f $(TARGET)

# Also remove the validation dumps.  The committed raw/validate.jsonl, raw/timing.jsonl and
# raw/checksums.sha256 stay.
distclean: clean
	rm -f $(RAWDIR)/val_*.bin $(RAWDIR)/checksums_rerun.txt

.PHONY: all run verify clean distclean
