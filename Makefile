CC=g++

CFLAGS= -Ofast -ffloat-store -c -Wall  -std=c++11
#CFLAGS= -c -g -O0 -Wall  -std=c++11

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: cr3bp

cr3bp: main.o body.o orbitalIntegrator.o
	$(CC) main.o body.o orbitalIntegrator.o -o cr3bp
#	$(CC) main.o -o cr3bp

main.o: src\main.cpp
	$(CC) $(CFLAGS) src\main.cpp

body.o: src\body.cpp
	$(CC) $(CFLAGS) src\body.cpp

orbitalIntegrator.o: src\orbitalIntegrator.cpp
	$(CC) $(CFLAGS) src\orbitalIntegrator.cpp


clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid
