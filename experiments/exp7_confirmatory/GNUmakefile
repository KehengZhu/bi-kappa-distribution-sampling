# Experiment 7 - pre-registered confirmatory test of the log-domain bi-Kappa loader.
#
# The same source is built once per standard library, because one of the claims under test
# is that it no longer matters which one ran it: within version 2.0.0 the stream is a
# function of the engine alone, so the two builds must agree byte for byte.
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
#   make p6             performance: timed blocks, B0-B4, seeds 7006-7010
#   make  why two builds of the candidate can still differ (see below)
#   make oracle         adjudicate every audit stream at 100 decimal digits
#   make analyze        results/*.csv, results/analysis_report.md, the gate table
#   make figures        figures/*.{pdf,svg,png}
#   make checksums      write both checksum manifests
#   make verify         re-check every checksum manifest
#   make reverify       regenerate the derived artifacts and check they are byte-identical
#
# A production run is:
#   make preflight && make selftest && make p1 p2 p3 p4 p5 p6 && make oracle && \
#   make analyze figures && make checksums && make verify
#
# Every simulation target depends on `preflight`.  In Experiment 6 preflight was a
# prerequisite of nothing, so a dirty tree ran to completion and only the manifest said so.
#
# Every simulation target accepts SMOKE=1, which shrinks the run and diverts *everything* it
# writes -- data, manifest and all -- under raw/smoke/.  SMOKE is forwarded to analyze and
# figures.  In Experiment 6 smoke output leaked into `make oracle`, `make checksums` and the
# production manifest; here nothing outside raw/smoke/ is touched by a smoke run, and the
# production targets never look inside it.

CXX_LIBCXX    := clang++
CXX_LIBSTDCXX := $(shell command -v g++-15 || command -v g++-14 || command -v g++-13)
BOOST_INC     ?= /opt/homebrew/include

# -ffp-contract=off: fused multiply-add changes the rounding of the very products whose
# overflow is under study, and it is applied inconsistently across inlining contexts.  With
# contraction on, the instrumented replica stops reproducing the released class bit for bit.
BASE_FLAGS    := -Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src
ORACLE_FLAGS  := -Wall -Wextra -std=c++14 -O2 -ffp-contract=off -I src -I $(BOOST_INC)

# The flags are compiled into the binary, so the manifest records what the build actually
# used rather than what the analysis later guesses it used.
CXXFLAGS      := $(BASE_FLAGS) -DEXP7_CXXFLAGS='"$(BASE_FLAGS)"'

# Note on build settings.  `-ffp-contract=off` is on purpose: a fused multiply-add rounds
# once instead of twice, and the sampler has several expressions a compiler may contract, so
# a contracting build differs from a non-contracting one by an ulp or two on a fraction of
# the components.  Pinning it keeps the stream a property of the source text.
#
# An earlier revision of this experiment carried a third, diagnostic build to locate a
# cross-compiler difference: clang fused the `sin(phi)`/`cos(phi)` pair in the released
# sampler into Darwin's `__sincos_stret`, whose sine differs from the standalone `sin` by an
# ulp on about one argument in a thousand, while GCC called the two separately.  The released
# header now draws its direction by a rejection method that uses no trigonometry at all, so
# the difference no longer exists and the diagnostic has been removed rather than kept as a
# permanently passing control.

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
N_P2          ?= 1000000
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
	$(CXX_LIBCXX) $(ORACLE_FLAGS) -stdlib=libc++ -o $@ $<

# The sampling sources must keep compiling as strict C++11 with no extensions, because the
# released loader is a header-only C++11 reference implementation.  The oracle is exempt by
# construction: it is a separate program.
cxx11-check: $(GENERATED)
	$(CXX_LIBCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@test -z "$(CXX_LIBSTDCXX)" || \
	  $(CXX_LIBSTDCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@echo "  strict C++11: ok"

# ---------------------------------------------------------------------------
# Gates that run before any data exists
# ---------------------------------------------------------------------------
selftest: $(TARGETS) cxx11-check
	@for b in $(TARGETS); do echo "== $$b"; ./$$b selftest || exit 1; done

# The comparator must still be the released 1.0.0 header, modulo the guard rename and the
# namespace wrap that letting both versions coexist in one translation unit requires.
check-legacy:
	@$(PY) src/check_legacy.py

# config/protocol.json must be exactly what config/make_protocol.py produces.  Regenerating
# into a scratch copy and diffing is what makes "the analysis reads the protocol" mean
# something: a hand-edited JSON would otherwise pass every hash check in the run.
protocol-check:
	@cp config/protocol.json config/.protocol.json.committed
	@if $(PY) config/make_protocol.py > /dev/null; then :; else \
	   cp config/.protocol.json.committed config/protocol.json; \
	   rm -f config/.protocol.json.committed; \
	   echo "  protocol-check: config/make_protocol.py failed to run"; \
	   exit 1; \
	 fi
	@if cmp -s config/protocol.json config/.protocol.json.committed; then \
	   rm -f config/.protocol.json.committed; \
	   echo "  protocol-check: config/protocol.json matches config/make_protocol.py"; \
	 else \
	   echo "  protocol-check: config/protocol.json is NOT what config/make_protocol.py produces"; \
	   diff config/.protocol.json.committed config/protocol.json | head -40; \
	   cp config/.protocol.json.committed config/protocol.json; \
	   rm -f config/.protocol.json.committed; \
	   exit 1; \
	 fi

preflight: $(TARGETS) exp7_oracle.exe protocol-check check-legacy
	$(PY) preflight.py $(if $(ALLOW_DIRTY_DEV),--allow-dirty,)

# ---------------------------------------------------------------------------
# Simulation phases.  Each depends on preflight, so a dirty dependency set stops the run
# before it produces data rather than after.
# ---------------------------------------------------------------------------
define RUN_PHASE
	@mkdir -p $(RAWROOT)/$(1)
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
# Adjudication
#
# Every audit stream the probe produced, re-derived at 100 decimal digits by a program that
# shares no arithmetic with it.  A single disagreement fails the target -- Experiment 6's
# equivalent ended in `|| echo`, so a disagreement printed a line and the target still
# succeeded.  Smoke streams are adjudicated only under SMOKE=1 and never alongside
# production ones.
# ---------------------------------------------------------------------------
# ORACLE_JOBS shards each audit stream across processes.  Sharding is a scheduling device
# only: shard i adjudicates the records congruent to i modulo N, so the union is exactly the
# unsharded run and every record is adjudicated exactly once.  Verified on the smoke stream --
# 2294 records unsharded, 574+574+573+573 across four shards.  At about 1.1 ms a record the
# frozen matrix is hours on one core, and there are only a handful of audit files, so
# parallelising over files alone would leave most of the machine idle.
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
# Analysis.  Written by a separate agent against results/schema.md; the targets below only
# invoke it, and forward SMOKE so that a smoke run can be taken end to end.
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
	@find config results figures src PROTOCOL.md GNUmakefile README.md preflight.py \
	   raw/environment.json raw/manifest.csv raw/raw_checksums.sha256 \
	   -type f ! -name 'checksums.sha256' ! -name '*.pyc' 2>/dev/null | sort | \
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

# Gate G6 requires that every derived artifact regenerate byte for byte from the same raw
# data.  This target does exactly that: it re-runs the analysis and the figures into a
# scratch tree and diffs.  It never touches raw/.
reverify:
	@rm -rf .reverify && mkdir -p .reverify
	@cp -R results .reverify/results_before
	@$(PY) analyze.py $(SMOKEARG)
	@$(PY) make_figures.py $(SMOKEARG)
	@if diff -r --brief .reverify/results_before results; then \
	   echo "  reverify: every derived artifact regenerated byte for byte"; \
	   rm -rf .reverify; \
	 else \
	   echo "  reverify: a derived artifact did not reproduce"; exit 1; \
	 fi

clean:
	rm -f $(TARGETS) exp7_oracle.exe $(GENERATED)
	rm -rf .reverify

distclean: clean
	rm -rf raw/smoke raw/p1 raw/p2 raw/p3 raw/p4 raw/p5 raw/p6 \
	       raw/manifest.csv raw/oracle_audit.jsonl raw/oracle_disagreements.jsonl \
	       raw/raw_checksums.sha256 checksums.sha256
