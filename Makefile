CC = gcc
CFLAGS =  -O2

all: discoal discoal_multipop
#
# executable 
#


discoal: discoal.c discoalFunctions.c discoal.h discoalFunctions.h alleleTraj.c
	$(CC) $(CFLAGS) -lm -o discoal discoal.c discoalFunctions.c ranlibComplete.c alleleTraj.c


discoal_multipop: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h
	$(CC) $(CFLAGS)  -o discoal_multipop discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm
#
# clean
#

clean:
	rm -f discoal discoal_multipop *.o

