# Makefile for validation tests

CC = gcc
CFLAGS = -Wall -g -I../../src/core -I../../extern/demes-c
LDFLAGS = -L../../build -lyaml -lm

# Source files from core
CORE_SRCS = ../../src/core/params.c \
            ../../src/core/yaml_loader.c \
            ../../src/core/demes_loader.c

# Object files
CORE_OBJS = $(CORE_SRCS:.c=.o)

# Demes library
DEMES_LIB = ../../extern/demes-c/libdemes.a

# Test programs
TESTS = test_sample_specification

all: $(TESTS)

test_sample_specification: test_sample_specification.o $(CORE_OBJS) $(DEMES_LIB)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(CORE_OBJS) $(TESTS)

.PHONY: all clean