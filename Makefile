CC = gcc
CFLAGS =  -O3

all: discoal
#
# executable 
#



discoal: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm
#
# clean
#

clean:
	rm -f discoal *.o

