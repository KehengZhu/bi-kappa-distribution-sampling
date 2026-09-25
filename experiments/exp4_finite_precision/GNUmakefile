# Experiment 4 - finite-precision losses of the direct (LEGACY) and the stabilized
# (CANDIDATE) radius calculation of the bi-Kappa loader.
#
# The same source is built once per standard library.  Within the 2.x line the stabilized
# loader's random stream is a function of the engine alone, so the two builds must agree
# byte for byte.
#
#   make cxx11-check    the sampling sources compile as strict C++11 under both compilers
#   make selftest       accounting identity, replica-vs-class equality, error paths, digests
#   make check-legacy   the vendored comparator is still the released 1.0.0 header
#   make protocol-check config/protocol.json is what config/make_protocol.py produces
#   make preflight      environment + provenance; refuses a dirty production run
#   make p1             scalar radial law, both methods, full ladder
#   make p2             mechanism: paired diagnostic layer + native layer
#   make p3             conditioning: per-attempt intended state vs success
#   make p4             complete anisotropic, rotated, capped 3-D loader, C0-C6
#   make p5             portability: CANDIDATE, full ladder, one tagged environment
#   make p6             performance: timed blocks, B0-B4, seeds 10006-10010
#   make oracle         recompute every recorded attempt at 100 decimal digits
#   make analyze        results/*.csv, results/analysis_report.md, the acceptance-criteria table
#   make figures        figures/*.{pdf,svg,png}
#   make checksums      write both checksum manifests
#   make verify         re-check every checksum manifest
#   make reverify       regenerate the derived artifacts and compare them with the committed ones
#                       (STRICT=1 for byte-for-byte; also accepted by protocol-check)
#
# A production run is:
#   make preflight && make selftest && make p1 p2 p3 p4 p5 p6 && make oracle && \
#   make analyze figures && make checksums && make verify
#
# Every simulation target depends on `preflight`, so a run with uncommitted changes to a
# source file stops before it writes any data.
#
# Every simulation target accepts SMOKE=1, which shrinks the run and diverts *everything* it
# writes -- data, manifest and all -- under raw/smoke/.  SMOKE is forwarded to analyze and
# figures.  A smoke run touches nothing outside raw/smoke/ and results/smoke/, and the
# production targets never look inside either.

CXX_LIBCXX    := clang++
CXX_LIBSTDCXX := $(shell command -v g++-15 || command -v g++-14 || command -v g++-13)
BOOST_INC     ?= /opt/homebrew/include

# -ffp-contract=off: fused multiply-add changes the rounding of the very products whose
# overflow is under study, and it is applied inconsistently across inlining contexts.  With
# contraction on, the instrumented replica stops reproducing the released class bit for bit.
BASE_FLAGS    := -Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src
ORACLE_FLAGS  := -Wall -Wextra -std=c++14 -O2 -ffp-contract=off -I src -I $(BOOST_INC)
MPFR_PREFIX   ?= /opt/homebrew
ORACLE_LIBS   := -L $(MPFR_PREFIX)/lib -lmpfr -lgmp

# The flags are compiled into the binary, so the manifest records what the build actually
# used rather than what the analysis later guesses it used.
CXXFLAGS      := $(BASE_FLAGS) -DEXP7_CXXFLAGS='"$(BASE_FLAGS)"'

# Note on build settings.  `-ffp-contract=off` is on purpose: a fused multiply-add rounds
# once instead of twice, and the sampler has several expressions a compiler may contract, so
# a contracting build differs from a non-contracting one by an ulp or two on a fraction of
# the components.  Pinning it keeps the stream a property of the source text.
#
# The released header draws its direction by a rejection method that uses no trigonometry.
# A `sin(phi)`/`cos(phi)` pair would make the two builds differ: clang fuses it into
# Darwin's `__sincos_stret`, whose sine differs from the standalone `sin` by an ulp on about
# one argument in a thousand, while GCC calls the two functions separately.

PY            := uv run --project ../../python python

PROBE_SRC     := src/exp7_probe.cpp
GENERATED     := src/exp7_protocol.H
PROBE_HDR     := src/exp7_common.H src/exp7_loaders.H src/exp7_digest.H $(GENERATED) \
                 src/legacy/bi_kappa_distribution_v1.H ../../cpp/bi_kappa_distribution.H

TARGETS       := exp7_probe_libcxx.exe
TAGS          := libcxx
ifneq ($(CXX_LIBSTDCXX),)
TARGETS       += exp7_probe_libstdcxx.exe
TAGS          += libstdcxx
endif

# Frozen production sizes (PROTOCOL.md Sec. 4).  Overriding them makes a run exploratory,
# and the probe records the size it was given in every row.
N_P1          ?= 1000000
N_P2          ?= 10000000
N_P3          ?= 1000000
N_P4          ?= 100000
N_P5          ?= 1000000
SEEDS_PROD    ?= 5

# The environment this build is executing on.  `translated` marks a Rosetta or emulated run,
# which PROTOCOL.md Sec. 5.3 records but excludes from the F7 decision.
EXECUTION     ?= native

ifdef SMOKE
SMOKEFLAG := --smoke
SMOKEARG  := --smoke
STAMP     := (smoke)
RAWROOT   := raw/smoke
else
SMOKEFLAG :=
SMOKEARG  :=
STAMP     :=
RAWROOT   := raw
endif

# STRICT=1 makes protocol-check and reverify require byte-for-byte equality.
STRICTARG := $(if $(STRICT),--strict)

.PHONY: all cxx11-check selftest preflight protocol-check check-legacy \
        p1 p2 p3 p4 p5 p6 oracle analyze figures checksums verify \
        reverify clean distclean

all: $(TARGETS) exp7_oracle.exe

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
$(GENERATED): src/gen_protocol_header.py config/protocol.json
	$(PY) src/gen_protocol_header.py

exp7_probe_libcxx.exe: $(PROBE_SRC) $(PROBE_HDR)
	$(CXX_LIBCXX) $(CXXFLAGS) -stdlib=libc++ -o $@ $(PROBE_SRC)

exp7_probe_libstdcxx.exe: $(PROBE_SRC) $(PROBE_HDR)
	$(CXX_LIBSTDCXX) $(CXXFLAGS) -o $@ $(PROBE_SRC)

exp7_oracle.exe: src/exp7_oracle.cpp src/exp7_common.H src/exp7_loaders.H src/exp7_digest.H \
                 $(GENERATED)
	$(CXX_LIBCXX) $(ORACLE_FLAGS) -I $(MPFR_PREFIX)/include -stdlib=libc++ -o $@ $< $(ORACLE_LIBS)

# The sampling sources must keep compiling as strict C++11 with no extensions, because the
# released loader is a header-only C++11 reference implementation.  The oracle is exempt by
# construction: it is a separate program.
cxx11-check: $(GENERATED)
	$(CXX_LIBCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@test -z "$(CXX_LIBSTDCXX)" || \
	  $(CXX_LIBSTDCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@echo "  strict C++11: ok"

# ---------------------------------------------------------------------------
# Checks that run before any data exists
# ---------------------------------------------------------------------------
selftest: $(TARGETS) cxx11-check
	@for b in $(TARGETS); do echo "== $$b"; ./$$b selftest || exit 1; done

# The comparator must still be the released 1.0.0 header, modulo the guard rename and the
# namespace wrap that letting both versions coexist in one translation unit requires.
check-legacy:
	@$(PY) src/check_legacy.py

# config/protocol.json must be what config/make_protocol.py produces.  Regenerating it and
# comparing is what makes "the analysis reads the protocol" mean something: a hand-edited
# JSON would otherwise pass every hash check in the run.  The comparison, and why its only
# random part is compared statistically, is described in compare_regenerated.py.  STRICT=1
# asks for byte-for-byte equality, which holds under the package versions recorded in
# raw/environment.json.  The committed file is never overwritten.
protocol-check:
	@$(PY) compare_regenerated.py protocol $(STRICTARG)

preflight: $(TARGETS) exp7_oracle.exe protocol-check check-legacy
	$(PY) preflight.py $(if $(ALLOW_DIRTY_DEV),--allow-dirty,)

# ---------------------------------------------------------------------------
# Simulation phases.  Each depends on preflight, so a dirty dependency set stops the run
# before it produces data rather than after.
# ---------------------------------------------------------------------------
define RUN_PHASE
	@mkdir -p $(RAWROOT)/$(1)
	@# Drop this phase's manifest rows before rewriting them.  The probe appends, so a
	@# re-run would otherwise keep the old rows beside the new ones and every total taken
	@# over the manifest would double-count.  Production rows are written once, from a
	@# clean tree; a smoke phase is usually re-run many times.  `phase` is column 3.
	@if [ -f $(RAWROOT)/manifest.csv ]; then \
	  head -1 $(RAWROOT)/manifest.csv > $(RAWROOT)/manifest.keep; \
	  awk -F, 'NR>1 && $$3 != "$(1)"' $(RAWROOT)/manifest.csv >> $(RAWROOT)/manifest.keep; \
	  mv $(RAWROOT)/manifest.keep $(RAWROOT)/manifest.csv; \
	fi
	@for b in $(TARGETS); do \
	  tag=$${b%.exe}; tag=$${tag#exp7_probe_}; \
	  echo "  $(1) [$$tag] $(STAMP)"; \
	  ./$$b $(1) --n $(2) --seeds $(3) --tag $$tag --execution $(EXECUTION) \
	     --out raw $(SMOKEFLAG) || exit 1; \
	done
endef

p1: $(TARGETS) preflight
	$(call RUN_PHASE,p1,$(N_P1),$(SEEDS_PROD))

p2: $(TARGETS) preflight
	$(call RUN_PHASE,p2,$(N_P2),$(SEEDS_PROD))

p3: $(TARGETS) preflight
	$(call RUN_PHASE,p3,$(N_P3),$(SEEDS_PROD))

p4: $(TARGETS) preflight
	$(call RUN_PHASE,p4,$(N_P4),$(SEEDS_PROD))

p5: $(TARGETS) preflight
	$(call RUN_PHASE,p5,$(N_P5),$(SEEDS_PROD))

p6: $(TARGETS) preflight
	$(call RUN_PHASE,p6,0,$(SEEDS_PROD))

# ---------------------------------------------------------------------------
# Recomputation at 100 digits
#
# Every attempt the probe recorded is recomputed at 100 decimal digits by a program that
# shares no arithmetic with the probe.  A single disagreement or conversion failure makes
# the target fail.  Smoke streams are recomputed only under SMOKE=1 and never alongside
# production ones.
# ---------------------------------------------------------------------------
# ORACLE_JOBS shards each audit stream across processes.  Sharding is a scheduling device
# only: shard i adjudicates the records congruent to i modulo N, so the union is exactly the
# unsharded run and every record is adjudicated exactly once.  Verified on the smoke stream --
# 2294 records unsharded, 574+574+573+573 across four shards.  At about 25 us a record the
# frozen matrix is about half an hour on one core, and there are only a handful of audit
# files, so parallelising over files alone would leave most of the machine idle.
ORACLE_JOBS ?= $(shell sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)

oracle: exp7_oracle.exe
	@mkdir -p results $(RAWROOT)
	@rm -f $(RAWROOT)/oracle_audit.jsonl $(RAWROOT)/oracle_disagreements.jsonl
	@found=0; rc=0; J=$(ORACLE_JOBS); \
	for f in $(RAWROOT)/p2/audit_*.bin $(RAWROOT)/p3/audit_*.bin; do \
	  test -f "$$f" || continue; found=1; \
	  echo "  oracle $$f  ($$J shards)"; \
	  pids=""; i=0; \
	  while [ $$i -lt $$J ]; do \
	    ( ./exp7_oracle.exe "$$f" --shard $$i --shards $$J \
	        --disagreements $(RAWROOT)/oracle_dis_$$i.tmp \
	        > $(RAWROOT)/oracle_aud_$$i.tmp ; echo $$? > $(RAWROOT)/oracle_rc_$$i.tmp ) & \
	    pids="$$pids $$!"; i=$$((i+1)); \
	  done; \
	  wait $$pids; \
	  i=0; \
	  while [ $$i -lt $$J ]; do \
	    cat $(RAWROOT)/oracle_aud_$$i.tmp >> $(RAWROOT)/oracle_audit.jsonl; \
	    test -f $(RAWROOT)/oracle_dis_$$i.tmp && \
	      cat $(RAWROOT)/oracle_dis_$$i.tmp >> $(RAWROOT)/oracle_disagreements.jsonl; \
	    if [ "$$(cat $(RAWROOT)/oracle_rc_$$i.tmp 2>/dev/null)" != "0" ]; then \
	      rc=1; echo "    DISAGREEMENTS OR CONVERSION FAILURES in $$f (shard $$i)"; fi; \
	    rm -f $(RAWROOT)/oracle_aud_$$i.tmp $(RAWROOT)/oracle_dis_$$i.tmp \
	          $(RAWROOT)/oracle_rc_$$i.tmp; \
	    i=$$((i+1)); \
	  done; \
	done; \
	if [ $$found -ne 1 ]; then \
	  echo "  no audit streams under $(RAWROOT) - run p2 and p3 first"; exit 1; fi; \
	echo "  oracle: $$(grep -c '"kind":"disagreement"' $(RAWROOT)/oracle_disagreements.jsonl || true) disagreement record(s)"; \
	exit $$rc

# ---------------------------------------------------------------------------
# Analysis.  analyze.py reads the raw records in the formats results/schema.md defines; the
# targets below invoke it and forward SMOKE so that a smoke run can be taken end to end.
# ---------------------------------------------------------------------------
analyze:
	$(PY) analyze.py $(SMOKEARG)

figures:
	$(PY) make_figures.py $(SMOKEARG)

# ---------------------------------------------------------------------------
# Provenance.  The production checksum bundle never descends into raw/smoke/.
# ---------------------------------------------------------------------------
checksums:
	@cd raw && find . -path ./smoke -prune -o -type f \
	   \( -name '*.bin' -o -name '*.jsonl' -o -name 'manifest.csv' \) -print | sort | \
	   xargs shasum -a 256 > raw_checksums.sha256 && \
	   echo "  wrote raw/raw_checksums.sha256 ($$(wc -l < raw_checksums.sha256) files)"
	@# `results/smoke/` is pruned for the same reason as `raw/smoke/`: smoke output is not
	@# part of the production record, and hashing it would make `make verify` fail after
	@# every smoke re-run.
	@find config results figures src PROTOCOL.md GNUmakefile README.md preflight.py \
	   raw/environment.json raw/manifest.csv raw/raw_checksums.sha256 \
	   \( -path './results/smoke' -o -path 'results/smoke' \) -prune -o \
	   -type f ! -name 'checksums.sha256' ! -name '*.pyc' -print 2>/dev/null | sort | \
	   xargs shasum -a 256 > checksums.sha256
	@echo "  wrote checksums.sha256 ($$(wc -l < checksums.sha256) files)"

verify:
	@fail=0; \
	if [ -f raw/raw_checksums.sha256 ]; then \
	  (cd raw && shasum -a 256 --check --quiet raw_checksums.sha256) || fail=1; \
	  echo "  raw/raw_checksums.sha256: $$( [ $$fail -eq 0 ] && echo ok || echo MISMATCH )"; \
	else echo "  raw/raw_checksums.sha256 missing"; fail=1; fi; \
	if [ -f checksums.sha256 ]; then \
	  shasum -a 256 --check --quiet checksums.sha256 || fail=1; \
	  echo "  checksums.sha256: $$( [ $$fail -eq 0 ] && echo ok || echo MISMATCH )"; \
	else echo "  checksums.sha256 missing"; fail=1; fi; \
	exit $$fail

# reverify checks that every derived artifact regenerates from the same raw data.  It copies
# results/ and figures/ to .reverify/, re-runs the analysis and the figures, moves the new
# files to .reverify/*_after, and puts the committed files back.  It never touches raw/.
#
# `results/` -- the CSVs, the JSON and the report -- and `figure_manifest.json`, which records
# each figure's source data, must agree with the committed files: every word, integer and
# row exactly, every floating-point number to a relative tolerance of 1e-9.  A different
# numpy or scipy moves the last digits of some p-values, and that is not a finding; a changed
# count or verdict is.  STRICT=1 requires byte-for-byte equality instead, which holds under
# the package versions recorded in raw/environment.json.  The rendered bytes of the PDF, SVG
# and PNG files are compared and reported but do not fail the target: no result rests on a
# rendered byte, and requiring one would tie the release to a particular matplotlib and font
# stack.  Comparing them at all means that a figure that stops reproducing -- an embedded
# timestamp, a hash-ordered element id, a font path -- is reported, rather than silently
# invalidating checksums.sha256, which hashes it.  On a failure .reverify/ is kept for
# inspection.
reverify:
	@rm -rf .reverify && mkdir -p .reverify
	@cp -R results .reverify/results_before
	@cp -R figures .reverify/figures_before
	@$(PY) analyze.py $(SMOKEARG) && $(PY) make_figures.py $(SMOKEARG); rc=$$?; \
	 mv results .reverify/results_after && mv figures .reverify/figures_after && \
	 cp -R .reverify/results_before results && cp -R .reverify/figures_before figures; \
	 if [ $$rc -ne 0 ]; then echo "  reverify: the analysis or the figures failed to run"; exit 1; fi
	@$(PY) compare_regenerated.py outputs $(STRICTARG) && rm -rf .reverify

clean:
	rm -f $(TARGETS) exp7_oracle.exe $(GENERATED)
	rm -rf .reverify

distclean: clean
	rm -rf raw/smoke raw/p1 raw/p2 raw/p3 raw/p4 raw/p5 raw/p6 \
	       raw/manifest.csv raw/oracle_audit.jsonl raw/oracle_disagreements.jsonl \
	       raw/raw_checksums.sha256 checksums.sha256
