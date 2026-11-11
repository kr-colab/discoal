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

# The sources of discoal (and of discoal_profile), and the files a build
# depends on
DISCOAL_SOURCES = $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/xoshiro256pp_compat.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(SRC_CORE)/shapes.c $(POOL_SOURCES) $(TSKIT_SOURCES)
DISCOAL_DEPS = $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_TSKIT)/tskitInterface.c $(SRC_TSKIT)/tskitInterface.h $(SRC_RNG)/xoshiro256pp_compat.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(SRC_CORE)/demesInterface.c $(SRC_CORE)/demesInterface.h $(SRC_CORE)/configInterface.c $(SRC_CORE)/configInterface.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c

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

# A link to build/discoal in the root, for commands and scripts that run
# ./discoal
symlinks: discoal
	@ln -sf build/discoal discoal
#
# executable 
#

discoal: libyaml demes-c libcyaml $(DISCOAL_DEPS)
	@mkdir -p build
	$(CC) $(CFLAGS) -DUSE_XOSHIRO256PP -o build/discoal $(DISCOAL_SOURCES) $(EXTERN_LIB) -lm -fcommon

# The release build with debug symbols and frame pointers, for profilers
# (perf call stacks, valgrind, heaptrack). The same optimization as discoal,
# so that the profile is of the code that users run.
discoal_profile: libyaml demes-c libcyaml $(DISCOAL_DEPS)
	@mkdir -p build
	$(CC) $(CFLAGS) -g -fno-omit-frame-pointer -DUSE_XOSHIRO256PP -o build/discoal_profile $(DISCOAL_SOURCES) $(EXTERN_LIB) -lm -fcommon

# Build debug version with ancestry verification
discoal_debug: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/ancestrySegment.h $(SRC_CORE)/activeSegment.c $(SRC_CORE)/activeSegment.h $(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c $(TSKIT_SOURCES)
	@mkdir -p build
	$(CC) -O2 $(INCLUDE) -DDEBUG_ANCESTRY -o build/discoal_debug $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_RNG)/ranlibComplete.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/ancestrySegment.c $(SRC_CORE)/activeSegment.c $(SRC_CORE)/shapes.c $(SRC_TSKIT)/tskitInterface.c $(SRC_CORE)/demesInterface.c $(SRC_CORE)/configInterface.c $(POOL_SOURCES) $(TSKIT_SOURCES) $(EXTERN_LIB) -lm -fcommon

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
# clean
#

clean:
	rm -rf build/*
	# Remove the link in the root (make symlinks)
	rm -f discoal
	# Clean libyaml
	@cd extern/libyaml && $(MAKE) -f Makefile.discoal clean
	# Clean demes-c
	@cd extern/demes-c && $(MAKE) -f Makefile.discoal clean
	# Clean libcyaml
	@cd extern/libcyaml && $(MAKE) -f Makefile.discoal clean
