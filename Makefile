CC=g++

CFLAGS=-c -Wall

all:
	g++ main.cpp ESTAR_init.cpp ESTAR.h -o ESTAR -lm

clean:
	rm -rf *o ESTAR

