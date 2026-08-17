CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2
TARGET   = exp2_sample.exe
SRC      = exp2_sample.cpp
RAWDIR   = raw

all: $(TARGET)

$(TARGET): $(SRC) ../../cpp/bi_kappa_distribution.H
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Regenerate every sample file from scratch and record their checksums.  The .bin
# files are gitignored on purpose; checksums.sha256 is whitelisted and is what makes
# a regenerated raw/ verifiable against the tracked results.
run: $(TARGET)
	mkdir -p $(RAWDIR)
	./$(TARGET) $(RAWDIR)
	cd $(RAWDIR) && shasum -a 256 manifest.csv run_*.bin > checksums.sha256

# Verify a regenerated raw/ against the committed checksums.
verify:
	cd $(RAWDIR) && shasum -a 256 --check checksums.sha256 | grep -v ': OK$$' || \
	  echo "all files match checksums.sha256"

clean:
	rm -f $(TARGET)

distclean: clean
	rm -rf $(RAWDIR)

.PHONY: all run verify clean distclean
