CC = gcc
CFLAGS =  -O2

all: discoal_multipop
#
# executable 
#



discoal_multipop: discoal_multipop.c discoalFunctions.c discoal.h discoalFunctions.h
	$(CC) $(CFLAGS)  -o discoal discoal_multipop.c discoalFunctions.c ranlibComplete.c alleleTraj.c -lm
#
# clean
#

clean:
	rm -f discoal discoal_multipop *.o

