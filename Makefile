INCLUDE = -I. \
	-I./src/core -I./src/rng \
	-I./src/tskit -I./extern/tskit -I./extern/tskit/kastore \
	-I./extern/demes-c -I./extern/libyaml -I./extern/libcyaml

EXTERN_LIB = \
	-Lextern/demes-c -Lextern/libyaml -Lextern/libcyaml \
	-ldemes -lcyaml -lyaml 

CC = gcc
CFLAGS = -O3 -march=native $(INCLUDE)
TEST_CFLAGS = -O2 $(INCLUDE)

# Source directories
SRC_CORE = src/core
SRC_RNG = src/rng
SRC_TSKIT = src/tskit

# Segment pool is now standard
POOL_SOURCES = $(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c

# Tskit source files
TSKIT_SOURCES = extern/tskit/tskit/core.c \
                extern/tskit/tskit/tables.c \
                extern/tskit/tskit/trees.c \
                extern/tskit/tskit/genotypes.c \
                extern/tskit/kastore/kastore.c

all: libyaml demes-c libcyaml discoal symlinks

# Build libyaml
libyaml:
	@echo "Building libyaml..."
	@cd extern/libyaml && $(MAKE) -f Makefile.discoal

# Build demes-c library
demes-c:
	@echo "Building demes-c library..."
	@cd extern/demes-c && $(MAKE) -f Makefile.discoal

# Build libcyaml
libcyaml:
	@echo "Building libcyaml..."
	@cd extern/libcyaml && $(MAKE) -f Makefile.discoal

# Create symlinks in root for backward compatibility with test scripts
symlinks: discoal
	@ln -sf build/discoal discoal
	@ln -sf build/discoal_edited discoal_edited
	@ln -sf build/discoal_legacy_backup discoal_legacy_backup
#
# executable 
#

discoal: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(CFLAGS) -DUSE_XOSHIRO256PP -o build/discoal $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/shapes.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

# The release build with debug symbols and frame pointers, for profilers
# (perf call stacks, valgrind, heaptrack). The same optimization as discoal,
# so that the profile is of the code that users run.
discoal_profile: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(CFLAGS) -g -fno-omit-frame-pointer -DUSE_XOSHIRO256PP -o build/discoal_profile $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/shapes.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

# Build with legacy L'Ecuyer RNG (for regression testing against master)
# This version uses the old RNG so outputs match master exactly for same seeds
discoal_legacy_rng: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o build/discoal_legacy_rng $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/ranlibComplete.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/shapes.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

# Build edited version for testing (same as main but explicit name)
discoal_edited: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(CFLAGS) -DUSE_XOSHIRO256PP -o build/discoal_edited $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/shapes.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

# Phase 3 regression reference: discoal binary at the SHA where Phase 3 began.
# This recipe checks out that SHA in a temp worktree, builds discoal there,
# copies the binary to build/discoal_pre_phase3, and removes the worktree.
# Phase 3 regression reference: discoal binary at the SHA where Phase 3 began.
# build/discoal_pre_phase3 is a real-file target so make skips the rebuild
# if it already exists. Remove the file to force a rebuild.
PHASE3_REF_SHA = 9b7c4e0
.PHONY: discoal_pre_phase3
discoal_pre_phase3: build/discoal_pre_phase3

build/discoal_pre_phase3:
	@mkdir -p build
	@WT=$$(mktemp -d) && \
	  git worktree add --detach "$$WT" $(PHASE3_REF_SHA) && \
	  $(MAKE) -C "$$WT" discoal && \
	  cp "$$WT/build/discoal" build/discoal_pre_phase3 && \
	  git worktree remove "$$WT"
	@echo "Built pre-Phase-3 reference: build/discoal_pre_phase3"

# Phase 5 regression reference: discoal binary at the SHA where Phase 5 began.
# build/discoal_pre_phase5 is a real-file target so make skips the rebuild
# if it already exists. Remove the file to force a rebuild.
PHASE5_REF_SHA = 12415d3
.PHONY: discoal_pre_phase5
discoal_pre_phase5: build/discoal_pre_phase5

build/discoal_pre_phase5:
	@mkdir -p build
	@WT=$$(mktemp -d) && \
	  git worktree add --detach "$$WT" $(PHASE5_REF_SHA) && \
	  $(MAKE) -C "$$WT" discoal && \
	  cp "$$WT/build/discoal" build/discoal_pre_phase5 && \
	  git worktree remove "$$WT"
	@echo "Built pre-Phase-5 reference: build/discoal_pre_phase5"

# Build debug version with ancestry verification
discoal_debug: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c $(TSKIT_SOURCES)
	@mkdir -p build
	$(CC) -O2 $(INCLUDE) -DDEBUG_ANCESTRY -o build/discoal_debug $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/ranlibComplete.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

# Build legacy version from master-backup branch for comparison testing
discoal_legacy_backup:
	@echo "Building legacy version from master-backup branch..."
	@mkdir -p /tmp/discoal_legacy_build
	@if git show master-backup:discoal_multipop.c > /tmp/discoal_legacy_build/discoal_multipop.c 2>/dev/null && \
	   git show master-backup:discoalFunctions.c > /tmp/discoal_legacy_build/discoalFunctions.c 2>/dev/null && \
	   git show master-backup:discoal.h > /tmp/discoal_legacy_build/discoal.h 2>/dev/null && \
	   git show master-backup:discoalFunctions.h > /tmp/discoal_legacy_build/discoalFunctions.h 2>/dev/null && \
	   git show master-backup:ranlibComplete.c > /tmp/discoal_legacy_build/ranlibComplete.c 2>/dev/null && \
	   git show master-backup:alleleTraj.c > /tmp/discoal_legacy_build/alleleTraj.c 2>/dev/null && \
	   git show master-backup:alleleTraj.h > /tmp/discoal_legacy_build/alleleTraj.h 2>/dev/null && \
	   git show master-backup:ranlib.h > /tmp/discoal_legacy_build/ranlib.h 2>/dev/null; then \
		mkdir -p $(CURDIR)/build && cd /tmp/discoal_legacy_build && $(CC) $(CFLAGS) -o discoal_legacy_backup discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm -fcommon && mv discoal_legacy_backup $(CURDIR)/build/ && cd ../..; \
		rm -rf /tmp/discoal_legacy_build; \
		echo "Legacy version built successfully from master-backup branch"; \
	else \
		echo "ERROR: Could not build legacy version from master-backup branch"; \
		echo "This is required for testing. Please ensure master-backup branch exists and contains the legacy code."; \
		rm -rf /tmp/discoal_legacy_build; \
		exit 1; \
	fi

# Build version from mem branch for comparison testing
discoal_mem_branch:
	@echo "Building version from mem branch..."
	@mkdir -p /tmp/discoal_mem_build
	@if git show mem:discoal_multipop.c > /tmp/discoal_mem_build/discoal_multipop.c 2>/dev/null && \
	   git show mem:discoalFunctions.c > /tmp/discoal_mem_build/discoalFunctions.c 2>/dev/null && \
	   git show mem:discoal.h > /tmp/discoal_mem_build/discoal.h 2>/dev/null && \
	   git show mem:discoalFunctions.h > /tmp/discoal_mem_build/discoalFunctions.h 2>/dev/null && \
	   git show mem:ranlibComplete.c > /tmp/discoal_mem_build/ranlibComplete.c 2>/dev/null && \
	   git show mem:alleleTraj.c > /tmp/discoal_mem_build/alleleTraj.c 2>/dev/null && \
	   git show mem:alleleTraj.h > /tmp/discoal_mem_build/alleleTraj.h 2>/dev/null && \
	   git show mem:ranlib.h > /tmp/discoal_mem_build/ranlib.h 2>/dev/null && \
	   git show mem:ancestrySegment.c > /tmp/discoal_mem_build/ancestrySegment.c 2>/dev/null && \
	   git show mem:ancestrySegment.h > /tmp/discoal_mem_build/ancestrySegment.h 2>/dev/null && \
	   git show mem:ancestrySegmentAVL.c > /tmp/discoal_mem_build/ancestrySegmentAVL.c 2>/dev/null && \
	   git show mem:ancestrySegmentAVL.h > /tmp/discoal_mem_build/ancestrySegmentAVL.h 2>/dev/null && \
	   git show mem:activeSegment.c > /tmp/discoal_mem_build/activeSegment.c 2>/dev/null && \
	   git show mem:activeSegment.h > /tmp/discoal_mem_build/activeSegment.h 2>/dev/null && \
	   git show mem:ancestryWrapper.h > /tmp/discoal_mem_build/ancestryWrapper.h 2>/dev/null && \
	   git show mem:ancestryVerify.h > /tmp/discoal_mem_build/ancestryVerify.h 2>/dev/null; then \
		mkdir -p $(CURDIR)/build && cd /tmp/discoal_mem_build && $(CC) $(CFLAGS) -o discoal_mem_branch discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c activeSegment.c -lm -fcommon && mv discoal_mem_branch $(CURDIR)/build/ && cd ../..; \
		rm -rf /tmp/discoal_mem_build; \
		echo "mem branch version built successfully"; \
	else \
		echo "ERROR: Could not build version from mem branch"; \
		echo "This is required for testing. Please ensure mem branch exists."; \
		rm -rf /tmp/discoal_mem_build; \
		exit 1; \
	fi

# make check: the checks that CI runs for every pull request, in one
# command (see testing/README.md). Needs Python with msprime, tskit, numpy,
# scipy, demes and stdpopsim. CHECK_WORKERS sets the parallel workers of the
# statistical checks.
CHECK_WORKERS ?= 4

check: discoal run_tests
	bash testing/yaml_parser_suite.sh
	bash testing/yaml_validation_suite.sh
	bash testing/preflight_validation_suite.sh
	bash testing/demes_validation_suite.sh
	bash testing/demes_output_comparison_suite.sh
	python3 testing/units_convention_check.py
	python3 testing/sweep_snp_check.py
	python3 testing/trajectory_file_check.py
	python3 testing/recurrent_sweep_check.py
	python3 testing/sweep_msprime_check.py --workers $(CHECK_WORKERS)
	python3 testing/gc_crossover_check.py --workers $(CHECK_WORKERS)
	python3 testing/linkage_msprime_check.py --workers $(CHECK_WORKERS)
	python3 testing/fuzz_msprime_check.py --models 150 --reps 500 --workers $(CHECK_WORKERS)

test: check

.PHONY: check test

# Unit test-related targets
TEST_DIR = test/unit
UNITY_DIR = extern/Unity/src

# Unity source files
UNITY_SOURCES = $(UNITY_DIR)/unity.c


# Individual test executables
test_node: $(TEST_DIR)/test_node.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_node $(TEST_DIR)/test_node.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_event: $(TEST_DIR)/test_event.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.c $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_event $(TEST_DIR)/test_event.c $(TEST_DIR)/test_globals.c \
		$(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c \
		$(SRC_TSKIT)/tskitInterface.c $(TSKIT_SOURCES) $(UNITY_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_shapes: $(TEST_DIR)/test_shapes.c $(SRC_CORE)/shapes.c $(SRC_CORE)/shapes.h $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -DUNITY_INCLUDE_DOUBLE -o build/test_shapes $(TEST_DIR)/test_shapes.c $(SRC_CORE)/shapes.c $(TEST_DIR)/test_globals.c $(SRC_RNG)/xoshiro256pp_compat.c $(UNITY_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_alleleTraj: $(TEST_DIR)/test_alleleTraj.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/alleleTraj.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUNITY_INCLUDE_DOUBLE -o build/test_alleleTraj $(TEST_DIR)/test_alleleTraj.c $(SRC_CORE)/alleleTraj.c $(SRC_RNG)/xoshiro256pp_compat.c $(UNITY_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_node_operations: $(TEST_DIR)/test_node_operations.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_node_operations $(TEST_DIR)/test_node_operations.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_mutations: $(TEST_DIR)/test_mutations.c $(SRC_CORE)/discoal.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_mutations $(TEST_DIR)/test_mutations.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_ancestry_segment: $(TEST_DIR)/test_ancestry_segment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/segmentPool.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_ancestry_segment $(TEST_DIR)/test_ancestry_segment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_fenwick: $(TEST_DIR)/test_fenwick.c $(SRC_CORE)/fenwick.c $(SRC_CORE)/fenwick.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_fenwick $(TEST_DIR)/test_fenwick.c $(SRC_CORE)/fenwick.c $(UNITY_SOURCES) -I$(UNITY_DIR)

test_active_segment: $(TEST_DIR)/test_active_segment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_active_segment $(TEST_DIR)/test_active_segment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_trajectory: $(TEST_DIR)/test_trajectory.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_trajectory $(TEST_DIR)/test_trajectory.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_config_interface: demes-c $(TEST_DIR)/test_config_interface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(TEST_DIR)/test_globals.c $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -DUNITY_INCLUDE_DOUBLE -o build/test_config_interface $(TEST_DIR)/test_config_interface.c $(SRC_CORE)/configInterface.c \
		$(SRC_CORE)/demesInterface.c $(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c \
		$(SRC_CORE)/segmentPool.c $(SRC_CORE)/fenwick.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c \
		$(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c $(TSKIT_SOURCES) \
		$(UNITY_SOURCES) -I$(UNITY_DIR) $(EXTERN_LIB) -lm -fcommon




# Run individual unit tests
run_tests: test_node test_event test_shapes test_alleleTraj test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory test_config_interface test_fenwick
	@echo "=== Running Unit Tests ==="
	./build/test_node
	./build/test_event
	./build/test_shapes
	./build/test_alleleTraj
	./build/test_node_operations
	./build/test_mutations
	./build/test_ancestry_segment
	./build/test_active_segment
	./build/test_trajectory
	./build/test_config_interface
	./build/test_fenwick
	@echo "=== All Unit Tests Passed ==="


#
# msUtils targets
#

niceStats: extern/msUtils/niceStats.c extern/msUtils/msGeneralStats.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o build/niceStats extern/msUtils/niceStats.c extern/msUtils/msGeneralStats.c -lm
	@ln -sf build/niceStats niceStats

# Build ms from reference_code for growth comparison testing
extern/ms: reference_code/msdir/ms.c reference_code/msdir/streec.c reference_code/msdir/rand1.c
	$(CC) -O3 -o extern/ms reference_code/msdir/ms.c reference_code/msdir/streec.c reference_code/msdir/rand1.c -lm

# Build alleleTrajTest utility (simple trajectory test program)
alleleTrajTest: test/alleleTrajTest.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/alleleTraj.h $(SRC_RNG)/ranlibComplete.c
	@mkdir -p build
	$(CC) -O3 -I. -I./src/core -I./src/rng -DSTANDALONE_ALLELETRAJ -o build/alleleTrajTest test/alleleTrajTest.c $(SRC_CORE)/alleleTraj.c $(SRC_RNG)/ranlibComplete.c -lm -fcommon

#
# clean
#

clean:
	rm -rf build/*
	rm -f discoaldoc.aux discoaldoc.bbl discoaldoc.blg discoaldoc.log discoaldoc.out
	# Remove any executables/symlinks in root (for backward compatibility)
	rm -f discoal discoal_edited discoal_legacy_backup discoal_mem_branch discoal_debug discoal_legacy_rng
	rm -f test_node test_event test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory
	rm -f alleleTrajTest niceStats
	# Clean libyaml
	@cd extern/libyaml && $(MAKE) -f Makefile.discoal clean
	# Clean demes-c
	@cd extern/demes-c && $(MAKE) -f Makefile.discoal clean
	# Clean libcyaml
	@cd extern/libcyaml && $(MAKE) -f Makefile.discoal clean
