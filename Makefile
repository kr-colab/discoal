CC = gcc
CFLAGS = -O3 -march=native -I. -I./extern/tskit -I./extern/tskit/kastore
TEST_CFLAGS = -O2 -I. -I./test/unit -I./extern/tskit -I./extern/tskit/kastore

# Tskit source files
TSKIT_SOURCES = extern/tskit/tskit/core.c \
                extern/tskit/tskit/tables.c \
                extern/tskit/tskit/trees.c \
                extern/tskit/tskit/genotypes.c \
                extern/tskit/kastore/kastore.c

all: discoal
#
# executable 
#

discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h tskitInterface.c tskitInterface.h $(TSKIT_SOURCES)
	$(CC) $(CFLAGS) -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c tskitInterface.c $(TSKIT_SOURCES) -lm -fcommon

# Build edited version for testing (same as main but explicit name)
discoal_edited: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h tskitInterface.c tskitInterface.h $(TSKIT_SOURCES)
	$(CC) $(CFLAGS) -o discoal_edited discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c tskitInterface.c $(TSKIT_SOURCES) -lm -fcommon

# Build debug version with ancestry verification
discoal_debug: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h $(TSKIT_SOURCES)
	$(CC) -O2 -I. -I./extern/tskit -I./extern/tskit/kastore -DDEBUG_ANCESTRY -o discoal_debug discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c tskitInterface.c $(TSKIT_SOURCES) -lm -fcommon

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
		cd /tmp/discoal_legacy_build && $(CC) $(CFLAGS) -o discoal_legacy_backup discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm -fcommon && mv discoal_legacy_backup $(CURDIR)/ && cd ../..; \
		rm -rf /tmp/discoal_legacy_build; \
		echo "Legacy version built successfully from master-backup branch"; \
	else \
		echo "ERROR: Could not build legacy version from master-backup branch"; \
		echo "This is required for testing. Please ensure master-backup branch exists and contains the legacy code."; \
		rm -rf /tmp/discoal_legacy_build; \
		exit 1; \
	fi

#
# tests - for comprehensive population genetics testing
#

test: discoal_legacy_backup discoal_edited
	@echo "Running comprehensive test suite..."
	cd testing && ./comprehensive_validation_suite.sh

test_quick: discoal_legacy_backup discoal_edited
	@echo "Running quick focused validation..."
	cd testing && ./focused_validation_suite.sh

test_reps: discoal_legacy_backup discoal_edited
	@echo "Running statistical validation with custom replicate count..."
	@read -p "Enter number of replicates (default 100): " reps; \
	reps=$${reps:-100}; \
	echo "Running with $$reps replicates..."; \
	cd testing && ./statistical_validation_suite.sh $$reps

test_msprime: discoal_edited
	@echo "Running msprime comparison suite..."
	cd testing && ./msprime_comparison_suite.sh

# Unit test-related targets
TEST_DIR = test/unit
UNITY_DIR = extern/Unity/src

# Unity source files
UNITY_SOURCES = $(UNITY_DIR)/unity.c

# Test source files
TEST_SOURCES = $(TEST_DIR)/test_runner.c

# Individual test executables
test_node: $(TEST_DIR)/test_node.c discoal.h ancestrySegment.o
	$(CC) $(TEST_CFLAGS) -o test_node $(TEST_DIR)/test_node.c $(UNITY_SOURCES) ancestrySegment.o -I$(UNITY_DIR) -fcommon

test_event: $(TEST_DIR)/test_event.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_event $(TEST_DIR)/test_event.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_node_operations: $(TEST_DIR)/test_node_operations.c discoal.h ancestrySegment.o
	$(CC) $(TEST_CFLAGS) -o test_node_operations $(TEST_DIR)/test_node_operations.c $(UNITY_SOURCES) ancestrySegment.o -I$(UNITY_DIR) -fcommon

test_mutations: $(TEST_DIR)/test_mutations.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_mutations $(TEST_DIR)/test_mutations.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_ancestry_segment: $(TEST_DIR)/test_ancestry_segment.c ancestrySegment.c ancestrySegmentAVL.c ancestrySegment.h ancestrySegmentAVL.h
	$(CC) $(TEST_CFLAGS) -o test_ancestry_segment $(TEST_DIR)/test_ancestry_segment.c ancestrySegment.c ancestrySegmentAVL.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_active_segment: $(TEST_DIR)/test_active_segment.c activeSegment.c activeSegment.h
	$(CC) $(TEST_CFLAGS) -o test_active_segment $(TEST_DIR)/test_active_segment.c activeSegment.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_trajectory: $(TEST_DIR)/test_trajectory.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_trajectory $(TEST_DIR)/test_trajectory.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

test_coalescence_recombination: $(TEST_DIR)/test_coalescence_recombination.c discoal.h ancestrySegment.o
	$(CC) $(TEST_CFLAGS) -o test_coalescence_recombination $(TEST_DIR)/test_coalescence_recombination.c $(UNITY_SOURCES) ancestrySegment.o -I$(UNITY_DIR) -fcommon

test_memory_management: $(TEST_DIR)/test_memory_management.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_memory_management $(TEST_DIR)/test_memory_management.c $(UNITY_SOURCES) -I$(UNITY_DIR) -fcommon

# Build the unified test runner executable
test_runner: $(TEST_SOURCES) test_node test_event test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory test_coalescence_recombination test_memory_management
	$(CC) $(TEST_CFLAGS) -o test_runner $(TEST_SOURCES) -I$(UNITY_DIR) -fcommon

# Object files needed by tests
ancestrySegment.o: ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h
	$(CC) $(TEST_CFLAGS) -c ancestrySegment.c ancestrySegmentAVL.c

# Run individual unit tests
run_tests: test_node test_event test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory test_coalescence_recombination test_memory_management
	@echo "=== Running Unit Tests ==="
	./test_node
	./test_event
	./test_node_operations
	./test_mutations
	./test_ancestry_segment
	./test_active_segment
	./test_trajectory
	./test_coalescence_recombination
	./test_memory_management
	@echo "=== All Unit Tests Passed ==="

# Run all tests using the unified runner
run_all_tests: test_runner
	./test_runner

#
# clean
#

clean:
	rm -f discoal discoal_edited discoal_legacy_backup *.o test_node test_event test_node_operations test_mutations test_ancestry_segment test_active_segment test_trajectory test_coalescence_recombination test_memory_management test_runner alleleTrajTest
	rm -f discoaldoc.aux discoaldoc.bbl discoaldoc.blg discoaldoc.log discoaldoc.out