CC = gcc
CFLAGS = -O2 -I.
TEST_CFLAGS = $(CFLAGS) -I./test/unit

all: discoal
#
# executable 
#



discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestryVerify.c ancestryVerify.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestryVerify.c -lm -fcommon

# Build edited version for testing (same as main but explicit name)
discoal_edited: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestryVerify.c ancestryVerify.h
	$(CC) $(CFLAGS)  -o discoal_edited discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestryVerify.c -lm -fcommon

# Build debug version with ancestry verification
discoal_debug: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestryVerify.c ancestryVerify.h
	$(CC) $(CFLAGS) -DDEBUG_ANCESTRY -o discoal_debug discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestryVerify.c -lm -fcommon

# Build legacy version from master branch for comparison testing
discoal_legacy_backup:
	@echo "Building legacy version from master branch..."
	@mkdir -p /tmp/discoal_legacy_build
	@if git show master:discoal_multipop.c > /tmp/discoal_legacy_build/discoal_multipop.c 2>/dev/null && \
	   git show master:discoalFunctions.c > /tmp/discoal_legacy_build/discoalFunctions.c 2>/dev/null && \
	   git show master:discoal.h > /tmp/discoal_legacy_build/discoal.h 2>/dev/null && \
	   git show master:discoalFunctions.h > /tmp/discoal_legacy_build/discoalFunctions.h 2>/dev/null && \
	   git show master:ranlibComplete.c > /tmp/discoal_legacy_build/ranlibComplete.c 2>/dev/null && \
	   git show master:alleleTraj.c > /tmp/discoal_legacy_build/alleleTraj.c 2>/dev/null && \
	   git show master:alleleTraj.h > /tmp/discoal_legacy_build/alleleTraj.h 2>/dev/null && \
	   git show master:ranlib.h > /tmp/discoal_legacy_build/ranlib.h 2>/dev/null; then \
		cd /tmp/discoal_legacy_build && $(CC) $(CFLAGS) -o discoal_legacy_backup discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm -fcommon && mv discoal_legacy_backup $(CURDIR)/ && cd ../..; \
		rm -rf /tmp/discoal_legacy_build; \
		echo "Legacy version built successfully from master branch"; \
	else \
		echo "ERROR: Could not build legacy version from master branch"; \
		echo "This is required for testing. Please ensure master branch exists and contains the legacy code."; \
		rm -rf /tmp/discoal_legacy_build; \
		exit 1; \
	fi

# Build both versions needed for testing
test_binaries: discoal_edited discoal_legacy_backup
	@echo "Built both edited and legacy versions for testing"

# Run the comprehensive testing suite
test_comprehensive: test_binaries
	@echo "Running comprehensive validation suite..."
	cd testing && ./comprehensive_validation_suite.sh

# Run the focused testing suite  
test_focused: test_binaries
	@echo "Running focused validation suite..."
	cd testing && ./focused_validation_suite.sh
	
test: alleleTrajTest.c alleleTraj.c alleleTraj.h discoalFunctions.c
	$(CC) $(CFLAGS)  -o alleleTrajTest alleleTrajTest.c alleleTraj.c ranlibComplete.c discoalFunctions.c -lm

# unit tests
test_node: test/unit/test_node.c discoalFunctions.c ranlibComplete.c alleleTraj.c discoal.h discoalFunctions.h
	$(CC) $(TEST_CFLAGS) -o test_node test/unit/test_node.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm

test_event: test/unit/test_event.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_event test/unit/test_event.c -lm

test_node_operations: test/unit/test_node_operations.c discoalFunctions.c ranlibComplete.c alleleTraj.c discoal.h discoalFunctions.h
	$(CC) $(TEST_CFLAGS) -o test_node_operations test/unit/test_node_operations.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm

test_mutations: test/unit/test_mutations.c discoalFunctions.c ranlibComplete.c alleleTraj.c discoal.h discoalFunctions.h
	$(CC) $(TEST_CFLAGS) -o test_mutations test/unit/test_mutations.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm

run_tests: test_node test_event test_node_operations test_mutations
	./test_node || exit 1
	./test_event || exit 1
	./test_node_operations || exit 1
	./test_mutations || exit 1

#
# clean
#

clean:
	rm -f discoal discoal_edited discoal_legacy_backup *.o test_node test_event test_node_operations test_mutations alleleTrajTest

