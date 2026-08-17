CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2
TARGET   = exp1_sample.exe
SRC      = exp1_sample.cpp
RAWDIR   = raw

all: $(TARGET)

$(TARGET): $(SRC) ../../cpp/bi_kappa_distribution.H
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Regenerate every sample file from scratch, then re-run the analysis.
run: $(TARGET)
	mkdir -p $(RAWDIR)
	./$(TARGET) $(RAWDIR)

clean:
	rm -f $(TARGET)

distclean: clean
	rm -rf $(RAWDIR)

.PHONY: all run clean distclean
