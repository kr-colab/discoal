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

# new unit tests
test_node: test/unit/test_node.c discoal.h
	$(CC) $(TEST_CFLAGS) -o test_node test/unit/test_node.c -lm

run_tests: test_node
	./test_node

#
# clean
#

clean:
	rm -f discoal *.o test_node alleleTrajTest

