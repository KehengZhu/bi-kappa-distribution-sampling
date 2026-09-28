CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2
TARGET   = exp2_sample.exe
SRC      = exp2_sample.cpp
RAWDIR   = raw

all: $(TARGET)

$(TARGET): $(SRC) ../../cpp/bi_kappa_distribution.H
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Regenerate every sample file from scratch.  The .bin files are not committed;
# raw/manifest.csv and raw/checksums.sha256 are, and raw/checksums.sha256 is what makes a
# regenerated raw/ verifiable.  `make run` therefore never overwrites it: the checksums of
# the new files go to raw/checksums_rerun.txt, which is not committed, for inspection.
run: $(TARGET)
	mkdir -p $(RAWDIR)
	./$(TARGET) $(RAWDIR)
	cd $(RAWDIR) && shasum -a 256 manifest.csv run_*.bin > checksums_rerun.txt

# Check raw/ against the committed raw/checksums.sha256.  Fails on any mismatch and on any
# missing file.
verify:
	@cd $(RAWDIR) && if shasum -a 256 --check --quiet checksums.sha256; then \
	  echo "all files match raw/checksums.sha256"; \
	else \
	  echo "raw/ does not match raw/checksums.sha256 (lines above name each file)" >&2; \
	  exit 1; \
	fi

clean:
	rm -f $(TARGET)

# Also remove the regenerated sample files.  The committed raw/manifest.csv and
# raw/checksums.sha256 stay.
distclean: clean
	rm -f $(RAWDIR)/run_*.bin $(RAWDIR)/checksums_rerun.txt

.PHONY: all run verify clean distclean
