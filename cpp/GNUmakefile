CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = main.exe
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC) test_suite.H bi_kappa_distribution.H bi_maxwellian_distribution.H \
          general_position_generator.H general_velocity_generator.H \
          field_aligned_velocity_generator.H
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: all clean
