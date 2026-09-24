CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2
PY       = uv run --project ../../python python

all: exp1_draw.exe

# One seed of the uncapped sampler at one kappa, written raw; driven by the scripts below.
exp1_draw.exe: exp1_draw.cpp ../../cpp/bi_kappa_distribution.H
	$(CXX) $(CXXFLAGS) -o exp1_draw.exe exp1_draw.cpp

# Cell-count test: nine kappa, 5 x 10^7 draws each (Table III, 400-cell column).
cells: exp1_draw.exe
	$(PY) exp1_cells.py

# Radius test on 14 shells, on the same draws (Table III radius column, Fig. 4).
tail: exp1_draw.exe
	$(PY) exp1_tail.py

# Variances at kappa = 2, 5, 10, on the same draws (Table IV).
moments: exp1_draw.exe
	$(PY) exp1_moments.py

# Invariance under rescaling and rotation, on the same seeds (Sec. VI C).
frame: exp1_draw.exe
	$(PY) exp1_frame.py

# Histograms of the three components (Fig. 5).
marginals: exp1_draw.exe
	$(PY) exp1_marginals.py

clean:
	rm -f exp1_draw.exe

.PHONY: all cells tail moments frame marginals clean
