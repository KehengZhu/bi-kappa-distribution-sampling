# Experiment 6 - low-kappa finite-precision validity and log-domain stabilization.
#
# The same source is built once per standard library, because the whole question is
# whether the answer depends on the library's Gamma generator.  Every simulation target
# accepts SMOKE=1, which shrinks the run to 1000 attempts x 1 seed and diverts its output
# to raw/smoke/ so that it can never be mistaken for production data.
#
#   make selftest     accounting identities, replica-vs-class equality, error-path fixtures
#   make preflight    environment + provenance manifest; refuses a dirty production run
#   make pilot        E1  select the LOG primitive
#   make baseline     E0  clean rerun of the QF/SPLIT envelope
#   make mechanism    E2  paired diagnostic layer + native layer
#   make oracle       audit the E2/E3 audit streams at 100 decimal digits
#   make conditioning E3  per-attempt intended state vs success
#   make loader       E4  complete anisotropic, rotated, capped 3-D loader
#   make portability  E5  native scalar envelope, one tagged environment
#   make performance  E6  timed blocks
#   make analyze      results/*.csv, results/exp6_results.json, validation matrix
#   make figures      figures/*.{pdf,svg,png}
#   make verify       re-check every checksum manifest
#
# A production run is:  make preflight && make selftest && make pilot && ... && make verify

CXX_LIBCXX    := clang++
CXX_LIBSTDCXX := $(shell command -v g++-15 || command -v g++-14 || command -v g++-13)
BOOST_INC     ?= /opt/homebrew/include

# -ffp-contract=off: fused multiply-add changes the rounding of the very products whose
# overflow is under study, and it is applied inconsistently across inlining contexts.  With
# contraction on, the instrumented replica stops reproducing the released class bit for bit.
CXXFLAGS      := -Wall -Wextra -std=c++11 -O2 -ffp-contract=off -I src
ORACLE_FLAGS  := -Wall -Wextra -std=c++14 -O2 -ffp-contract=off -I src -I $(BOOST_INC)

PROBE_SRC := src/exp6_probe.cpp
PROBE_HDR := src/exp6_common.H src/exp6_loaders.H src/log_gamma_identity.H \
             src/log_gamma_published.H ../../cpp/bi_kappa_distribution.H

# analyze.py and make_figures.py import these; preflight tracks them as dependencies.
PY_MODULES := exp6_stats.py exp6_io.py

PY        := uv run --project ../../python python

TARGETS   := exp6_probe_libcxx.exe
TAGS      := libcxx
ifneq ($(CXX_LIBSTDCXX),)
TARGETS   += exp6_probe_libstdcxx.exe
TAGS      += libstdcxx
endif

# Frozen production sizes (plan Sec. 7).  Overriding these makes a run exploratory.
N_BASELINE     ?= 1000000
N_PILOT        ?= 200000
N_MECHANISM    ?= 1000000
N_CONDITIONING ?= 1000000
N_LOADER       ?= 100000
N_PORTABILITY  ?= 1000000

SEEDS_PROD  ?= 5
SEEDS_PILOT ?= 3

ifdef SMOKE
SMOKEFLAG := --smoke
STAMP     := (smoke)
else
SMOKEFLAG :=
STAMP     :=
endif

LOG_PRIMITIVE ?= $(shell test -f config/frozen.json && \
                   sed -n 's/.*"log_primitive"[ ]*:[ ]*"\([^"]*\)".*/\1/p' config/frozen.json)
ifeq ($(strip $(LOG_PRIMITIVE)),)
LOG_PRIMITIVE := LOG-ID
endif

.PHONY: all selftest preflight pilot baseline mechanism oracle conditioning loader \
        portability performance analyze figures verify checksums clean distclean cxx11-check \
        check-frozen baseline-contract-check

all: $(TARGETS) exp6_oracle.exe

exp6_probe_libcxx.exe: $(PROBE_SRC) $(PROBE_HDR)
	$(CXX_LIBCXX) $(CXXFLAGS) -stdlib=libc++ -o $@ $(PROBE_SRC)

exp6_probe_libstdcxx.exe: $(PROBE_SRC) $(PROBE_HDR)
	$(CXX_LIBSTDCXX) $(CXXFLAGS) -o $@ $(PROBE_SRC)

exp6_oracle.exe: src/exp6_oracle.cpp src/exp6_common.H
	$(CXX_LIBCXX) $(ORACLE_FLAGS) -stdlib=libc++ -o $@ $<

# The sampling sources must keep compiling as strict C++11 with no extensions, because the
# released loader is a header-only C++11 reference implementation.  The oracle is exempt by
# construction: it is a separate program.
cxx11-check:
	$(CXX_LIBCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@test -z "$(CXX_LIBSTDCXX)" || \
	  $(CXX_LIBSTDCXX) -std=c++11 -pedantic-errors -Wall -Wextra -fsyntax-only -I src $(PROBE_SRC)
	@echo "  strict C++11: ok"

selftest: $(TARGETS) cxx11-check
	@for b in $(TARGETS); do echo "== $$b"; ./$$b selftest || exit 1; done

preflight: $(TARGETS) exp6_oracle.exe
	$(PY) preflight.py $(if $(ALLOW_DIRTY_DEV),--allow-dirty,)

# Once the pilot has frozen a primitive, every later phase refuses to run against a
# modified source for it.  The hash is what is checked, not the name: renaming the file or
# editing it after the selection would otherwise pass silently.
check-frozen:
	@python3 -c "import json,hashlib,sys; \
	  c=json.load(open('config/frozen.json')); \
	  sys.exit(0) if not c.get('frozen') else None; \
	  h=hashlib.sha256(open(c['source_file'],'rb').read()).hexdigest(); \
	  sys.exit(0) if h==c['source_sha256'] else (sys.stderr.write( \
	    'exp6: %s has changed since the pilot froze it\n  frozen  %s\n  current %s\n' \
	    % (c['source_file'], c['source_sha256'], h)) or sys.exit(1))"

define RUN_PHASE
	@mkdir -p raw/$(1) raw/smoke
	@for b in $(TARGETS); do \
	  tag=$${b%.exe}; tag=$${tag#exp6_probe_}; \
	  echo "  $(1) [$$tag] $(STAMP)"; \
	  ./$$b $(1) --n $(2) --seeds $(3) --tag $$tag --log-primitive $(LOG_PRIMITIVE) \
	     --out raw $(SMOKEFLAG) || exit 1; \
	done
endef

pilot: $(TARGETS)
	$(call RUN_PHASE,pilot,$(N_PILOT),$(SEEDS_PILOT))

baseline: $(TARGETS)
	$(call RUN_PHASE,baseline,$(N_BASELINE),$(SEEDS_PROD))

# Experiment 4 was built without an explicit -ffp-contract setting, so it inherited the
# compiler default: fused multiply-add enabled.  This target rebuilds the libstdc++ probe
# the same way and reruns the single-precision baseline, so that any residual difference
# between the clean rerun and the historical numbers can be attributed to that setting
# rather than left unexplained.  It is a provenance control, not a production row.
exp6_probe_contractfast.exe: $(PROBE_SRC) $(PROBE_HDR)
	$(CXX_LIBSTDCXX) -Wall -Wextra -std=c++11 -O2 -I src -o $@ $(PROBE_SRC)

baseline-contract-check: exp6_probe_contractfast.exe
	@mkdir -p raw/baseline
	./exp6_probe_contractfast.exe baseline --n $(N_BASELINE) --seeds $(SEEDS_PROD) \
	   --precision float --tag libstdcxx-contractfast --out raw

mechanism: $(TARGETS) check-frozen
	$(call RUN_PHASE,mechanism,$(N_MECHANISM),$(SEEDS_PROD))

conditioning: $(TARGETS) check-frozen
	$(call RUN_PHASE,conditioning,$(N_CONDITIONING),$(SEEDS_PROD))

loader: $(TARGETS) check-frozen
	$(call RUN_PHASE,loader,$(N_LOADER),$(SEEDS_PROD))

portability: $(TARGETS) check-frozen
	$(call RUN_PHASE,portability,$(N_PORTABILITY),$(SEEDS_PROD))

performance: $(TARGETS) check-frozen
	$(call RUN_PHASE,performance,0,$(SEEDS_PROD))

# Every audit stream the probe produced, re-derived at 100 decimal digits by a program
# that shares no arithmetic with it.  A single disagreement fails the target.
oracle: exp6_oracle.exe
	@mkdir -p results
	@rm -f results/oracle_audit.jsonl results/oracle_disagreements.jsonl
	@found=0; \
	for f in raw/mechanism/audit_*.bin raw/conditioning/audit_*.bin raw/smoke/audit_*.bin; do \
	  test -f "$$f" || continue; found=1; \
	  echo "  oracle $$f"; \
	  ./exp6_oracle.exe "$$f" --disagreements results/oracle_disagreements_tmp.jsonl \
	      >> results/oracle_audit.jsonl || echo "    DISAGREEMENTS FOUND in $$f"; \
	  test -f results/oracle_disagreements_tmp.jsonl && \
	      cat results/oracle_disagreements_tmp.jsonl >> results/oracle_disagreements.jsonl; \
	  rm -f results/oracle_disagreements_tmp.jsonl; \
	done; \
	rm -f results/oracle_disagreements_tmp.jsonl; \
	test $$found -eq 1 || { echo "  no audit streams found - run mechanism/conditioning first"; exit 1; }
	@touch results/oracle_disagreements.jsonl
	@echo "  oracle: $$(wc -l < results/oracle_disagreements.jsonl) disagreement(s)"

analyze:
	$(PY) analyze.py

figures:
	$(PY) make_figures.py

checksums:
	@cd raw && find . -type f \( -name '*.bin' -o -name '*.jsonl' \) | sort | \
	   xargs shasum -a 256 > raw_checksums.sha256 && \
	   echo "  wrote raw/raw_checksums.sha256 ($$(wc -l < raw_checksums.sha256) files)"
	@find config results figures src GNUmakefile analyze.py make_figures.py preflight.py \
	   exp6_stats.py exp6_io.py README.md raw/environment.json raw/manifest.csv \
	   raw/raw_checksums.sha256 -type f ! -name 'checksums.sha256' 2>/dev/null | sort | \
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

clean:
	rm -f $(TARGETS) exp6_oracle.exe

distclean: clean
	rm -rf raw/smoke raw/pilot raw/baseline raw/mechanism raw/conditioning raw/loader \
	       raw/portability raw/performance raw/raw_checksums.sha256 checksums.sha256
