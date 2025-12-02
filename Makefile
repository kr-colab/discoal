CC = gcc
CFLAGS = -O3 -march=native -I. -I./src/core -I./src/rng -I./src/tskit -I./extern/tskit -I./extern/tskit/kastore -I./extern/demes-c
TEST_CFLAGS = -O2 -I. -I./src/core -I./src/rng -I./src/tskit -I./test/unit -I./extern/tskit -I./extern/tskit/kastore -I./extern/demes-c

# Source directories
SRC_CORE = src/core
SRC_RNG = src/rng
SRC_TSKIT = src/tskit

# Segment pool is now standard
POOL_SOURCES = $(SRC_CORE)/segmentPool.c

# Tskit source files
TSKIT_SOURCES = extern/tskit/tskit/core.c \
                extern/tskit/tskit/tables.c \
                extern/tskit/tskit/trees.c \
                extern/tskit/tskit/genotypes.c \
                extern/tskit/kastore/kastore.c

all: demes-c discoal symlinks

# Build demes-c library
demes-c:
	@echo "Building demes-c library..."
	@cd extern/demes-c && $(MAKE) -f Makefile.discoal

# Create symlinks in root for backward compatibility with test scripts
symlinks: discoal
	@ln -sf build/discoal discoal
	@ln -sf build/discoal_edited discoal_edited
	@ln -sf build/discoal_legacy_backup discoal_legacy_backup
#
# executable 
#

discoal: demes-c $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/ancestrySegmentAVL.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h
	@mkdir -p build
	$(CC) $(CFLAGS) -DUSE_XOSHIRO256PP -o build/discoal $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(POOL_SOURCES) $(TSKIT_SOURCES) -Lextern/demes-c -ldemes -lyaml -lm -fcommon

# Build with legacy L'Ecuyer RNG (for regression testing against master)
# This version uses the old RNG so outputs match master exactly for same seeds
discoal_legacy_rng: demes-c $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/ancestrySegmentAVL.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h
	@mkdir -p build
	$(CC) $(CFLAGS) -o build/discoal_legacy_rng $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/ranlibComplete.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(POOL_SOURCES) $(TSKIT_SOURCES) -Lextern/demes-c -ldemes -lyaml -lm -fcommon

# Build edited version for testing (same as main but explicit name)
discoal_edited: demes-c $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/ancestrySegmentAVL.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h
	@mkdir -p build
	$(CC) $(CFLAGS) -DUSE_XOSHIRO256PP -o build/discoal_edited $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(POOL_SOURCES) $(TSKIT_SOURCES) -Lextern/demes-c -ldemes -lyaml -lm -fcommon

# Build debug version with ancestry verification
discoal_debug: $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/ancestrySegmentAVL.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(TSKIT_SOURCES)
	@mkdir -p build
	$(CC) -O2 -I. -I./src/core -I./src/rng -I./src/tskit -I./extern/tskit -I./extern/tskit/kastore -DDEBUG_ANCESTRY -o build/discoal_debug $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/ranlibComplete.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(TSKIT_SOURCES) -lm -fcommon

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

#
# tests - for comprehensive population genetics testing
#

test: discoal_legacy_backup discoal_edited
	@echo "Running comprehensive test suite..."
	@ln -sf build/discoal_legacy_backup discoal_legacy_backup
	@ln -sf build/discoal_edited discoal_edited
	cd testing && ./comprehensive_validation_suite.sh

test_comprehensive: discoal_legacy_backup discoal_edited
	@echo "Running comprehensive test suite..."
	@ln -sf build/discoal_legacy_backup discoal_legacy_backup
	@ln -sf build/discoal_edited discoal_edited
	cd testing && ./comprehensive_validation_suite.sh

test_quick: discoal_legacy_backup discoal_edited
	@echo "Running quick focused validation..."
	@ln -sf build/discoal_legacy_backup discoal_legacy_backup
	@ln -sf build/discoal_edited discoal_edited
	cd testing && ./focused_validation_suite.sh

test_reps: discoal_legacy_backup discoal_edited
	@echo "Running statistical validation with custom replicate count..."
	@read -p "Enter number of replicates (default 100): " reps; \
	reps=$${reps:-100}; \
	echo "Running with $$reps replicates..."; \
	cd testing && ./statistical_validation_suite.sh $$reps

# Build version from HEAD of current branch as legacy_backup for comparison
discoal_legacy_backup_head:
	@echo "Building version from HEAD of current branch as legacy_backup..."
	@mkdir -p /tmp/discoal_head_build
	@git archive HEAD | tar -x -C /tmp/discoal_head_build
	@cd /tmp/discoal_head_build && $(CC) $(CFLAGS) -o discoal_legacy_backup discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c tskitInterface.c $(TSKIT_SOURCES) -lm -fcommon && mv discoal_legacy_backup $(CURDIR)/
	@rm -rf /tmp/discoal_head_build
	@echo "HEAD version built successfully as discoal_legacy_backup"

# Build both versions for HEAD comparison testing
test_binaries_head: discoal_edited discoal_legacy_backup_head
	@echo "Built both optimized and HEAD versions for testing"

# Run comprehensive test comparing optimized vs HEAD of branch
test_comprehensive_head: test_binaries_head
	@echo "Running comprehensive validation suite (optimized vs HEAD)..."
	@echo "This will compare the performance improvements of our optimizations"
	@echo "Optimized = current working directory with trajectory optimizations"
	@echo "Legacy = HEAD of current branch (before optimizations)"
	cd testing && ./comprehensive_validation_suite.sh

test_binaries_mem: discoal_edited discoal_mem_branch
	@echo "Built both optimized and mem branch versions for testing"
	cp build/discoal_mem_branch build/discoal_legacy_backup 
	@echo "copying discoal_mem_branch to discoal_legacy_backup for testing"

test_comprehensive_mem: test_binaries_mem
	@echo "Running comprehensive validation suite (optimized vs HEAD)..."
	@echo "This will compare the performance improvements of our optimizations"
	@echo "Optimized = current working directory with trajectory optimizations"
	@echo "Legacy = mem branch (before optimizations)"
	cd testing && ./comprehensive_validation_suite.sh

test_msprime: discoal_edited
	@echo "Running msprime comparison suite..."
	cd testing && ./msprime_comparison_suite.sh

test_nicestats: discoal_legacy_backup discoal_edited niceStats
	@echo "Running niceStats comparison suite (1000 replicates in 10 chunks)..."
	cd testing && ./nicestats_comparison_suite.sh 1000 50

test_nicestats_quick: discoal_legacy_backup discoal_edited niceStats
	@echo "Running quick niceStats comparison (100 replicates in 5 chunks)..."
	cd testing && ./nicestats_comparison_suite.sh 100 5

test_nicestats_large: discoal_legacy_backup discoal_edited niceStats
	@echo "Running large niceStats comparison (10000 replicates in 20 chunks)..."
	cd testing && ./nicestats_comparison_suite.sh 10000 20

test_nicestats_mem: discoal_mem_branch discoal_edited niceStats
	@echo "Running niceStats comparison suite (current vs mem branch)..."
	cd testing && ./nicestats_comparison_suite.sh 1000 10 ../discoal_mem_branch

test_nicestats_xoshiro: discoal discoal_xoshiro niceStats
	@echo "Running niceStats comparison suite (legacy RNG vs xoshiro256++)..."
	@# Use discoal as legacy and discoal_xoshiro as edited for comparison
	@cp discoal discoal_legacy_backup
	@cp discoal_xoshiro discoal_edited
	cd testing && ./nicestats_comparison_suite.sh 1000 50
	@# Clean up temporary copies

test_growth: discoal extern/ms niceStats
	@echo "Running exponential growth comparison against ms baseline..."
	cd testing && ./growth_comparison_suite.sh
	@rm -f discoal_legacy_backup discoal_edited

# Compare current tskit-integration branch with mem branch
test_tskit_vs_mem: discoal_edited discoal_mem_branch
	@echo "Comparing tskit-integration branch (current) vs mem branch..."
	@echo "This will show the differences between the memory optimizations and tskit integration"
	@# Temporarily rename binaries for the test suite
	@mv discoal_mem_branch discoal_legacy_backup
	cd testing && ./comprehensive_validation_suite.sh
	@# Restore original name
	@mv discoal_legacy_backup discoal_mem_branch

# Unit test-related targets
TEST_DIR = test/unit
UNITY_DIR = extern/Unity/src

# Unity source files
UNITY_SOURCES = $(UNITY_DIR)/unity.c


# Individual test executables
test_node: $(TEST_DIR)/test_node.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_node $(TEST_DIR)/test_node.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_event: $(TEST_DIR)/test_event.c $(SRC_CORE)/discoal.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_event $(TEST_DIR)/test_event.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_node_operations: $(TEST_DIR)/test_node_operations.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_node_operations $(TEST_DIR)/test_node_operations.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_mutations: $(TEST_DIR)/test_mutations.c $(SRC_CORE)/discoal.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_mutations $(TEST_DIR)/test_mutations.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_ancestry_segment: $(TEST_DIR)/test_ancestry_segment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/ancestrySegmentAVL.h $(SRC_CORE)/segmentPool.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_ancestry_segment $(TEST_DIR)/test_ancestry_segment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_active_segment: $(TEST_DIR)/test_active_segment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -o build/test_active_segment $(TEST_DIR)/test_active_segment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_trajectory: $(TEST_DIR)/test_trajectory.c $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_trajectory $(TEST_DIR)/test_trajectory.c $(UNITY_SOURCES) \
		$(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c \
		$(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c \
		$(TSKIT_SOURCES) -I$(UNITY_DIR) -lm -fcommon

test_config_interface: demes-c $(TEST_DIR)/test_config_interface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -o build/test_config_interface $(TEST_DIR)/test_config_interface.c $(SRC_CORE)/configInterface.c \
		$(SRC_CORE)/demesInterface.c $(TEST_DIR)/test_globals.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/ancestrySegment.c \
		$(SRC_CORE)/ancestrySegmentAVL.c $(SRC_CORE)/segmentPool.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c \
		$(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(TSKIT_SOURCES) \
		$(UNITY_SOURCES) -I$(UNITY_DIR) -Lextern/demes-c -ldemes -lyaml -lm -fcommon




# Run individual unit tests
run_tests: test_node test_event test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory test_config_interface
	@echo "=== Running Unit Tests ==="
	./build/test_node
	./build/test_event
	./build/test_node_operations
	./build/test_mutations
	./build/test_ancestry_segment
	./build/test_active_segment
	./build/test_trajectory
	./build/test_config_interface
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
	# Clean demes-c
	@cd extern/demes-c && $(MAKE) -f Makefile.discoal clean