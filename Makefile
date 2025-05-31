CC = gcc
CFLAGS = -O2 -I.
TEST_CFLAGS = $(CFLAGS) -I./test/unit

all: discoal
#
# executable 
#



discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm -fcommon
	
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
	rm -f discoal *.o test_node test_event test_node_operations test_mutations alleleTrajTest

