CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = main.exe
SRC = main.cpp
TEST_TARGET = test.exe
TEST_SRC = test_main.cpp

HEADERS = test_suite.H bi_kappa_distribution.H bi_maxwellian_distribution.H \
          general_position_generator.H general_velocity_generator.H \
          field_aligned_velocity_generator.H

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Build and run the regression suite on its own.  main.exe also runs it, but then
# writes ~25 MB of sample files, so this is the target continuous integration uses.
# The exit status is the suite's: 0 if every check passed.
$(TEST_TARGET): $(TEST_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_SRC)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET)

.PHONY: all clean test
