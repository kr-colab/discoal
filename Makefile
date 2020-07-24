CC = gcc
CFLAGS = -O2

all: discoal
#
# executable 
#



discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm -fcommon
	
test: alleleTrajTest.c alleleTraj.c alleleTraj.h discoalFunctions.c
	$(CC) $(CFLAGS)  -o alleleTrajTest alleleTrajTest.c alleleTraj.c ranlibComplete.c discoalFunctions.c -lm
#
# clean
#

clean:
	rm -f discoal *.o

