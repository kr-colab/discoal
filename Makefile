CC = gcc
CFLAGS = -O2 -I.
TEST_CFLAGS = $(CFLAGS) -I./test/unit

all: discoal
#
# executable 
#



discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c -lm -fcommon

# Build edited version for testing (same as main but explicit name)
discoal_edited: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h
	$(CC) -O3 -march=native -I.  -o discoal_edited discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c -lm -fcommon

# Build debug version with ancestry verification
discoal_debug: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h ancestrySegment.c ancestrySegment.h ancestrySegmentAVL.c ancestrySegmentAVL.h ancestryVerify.c ancestryVerify.h activeSegment.c activeSegment.h
	$(CC) $(CFLAGS) -DDEBUG_ANCESTRY -o discoal_debug discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c -lm -fcommon

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

# Build version from HEAD of current branch as legacy_backup for comparison
discoal_legacy_backup_head:
	@echo "Building version from HEAD of current branch as legacy_backup..."
	@mkdir -p /tmp/discoal_head_build
	@git archive HEAD | tar -x -C /tmp/discoal_head_build
	@cd /tmp/discoal_head_build && $(CC) $(CFLAGS) -o discoal_legacy_backup discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c ancestrySegment.c ancestrySegmentAVL.c ancestryVerify.c activeSegment.c -lm -fcommon && mv discoal_legacy_backup $(CURDIR)/
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

# Generate PDF documentation from LaTeX source
doc: discoaldoc.pdf

discoaldoc.pdf: discoaldoc.tex texrefs.bib
	@echo "Generating PDF documentation..."
	@if command -v pdflatex >/dev/null 2>&1; then \
		pdflatex discoaldoc.tex && \
		bibtex discoaldoc && \
		pdflatex discoaldoc.tex && \
		pdflatex discoaldoc.tex && \
		rm -f discoaldoc.aux discoaldoc.bbl discoaldoc.blg discoaldoc.log discoaldoc.out && \
		echo "Documentation generated: discoaldoc.pdf"; \
	else \
		echo "ERROR: pdflatex not found. Please install LaTeX (e.g., texlive-latex-base texlive-latex-extra)"; \
		exit 1; \
	fi
	
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
	rm -f discoaldoc.aux discoaldoc.bbl discoaldoc.blg discoaldoc.log discoaldoc.out

