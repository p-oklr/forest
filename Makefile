# Makefile for the project
SHELL := /bin/bash
CC=mpif90 -fdefault-real-8 -fdefault-integer-8
CFLAGS=-O2
EXEC=main.bin

all: $(EXEC)

init:
	mkdir -p data/ figures/

$(EXEC): parameters.o random.o models.o diagnostics.o main.o
	$(CC) $^ -o $(EXEC)
	git rev-parse HEAD > git_hash


%.o: %.f03
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.mod
