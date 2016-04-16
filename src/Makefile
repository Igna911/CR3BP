CC=g++

CFLAGS= -Ofast -ffloat-store -lquadmath -c -Wall  -std=c++11

SRCDIR = /source/ODIS/
BUILDDIR = /source/build/

all: cr3bp

cr3bp: main.o body.o orbitalIntegrator.o
	$(CC) main.o body.o orbitalIntegrator.o -o cr3bp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

body.o: body.cpp
	$(CC) $(CFLAGS) body.cpp

orbitalIntegrator.o: orbitalIntegrator.cpp
		$(CC) $(CFLAGS) orbitalIntegrator.cpp


clean:
	rm -r *o ODIS NorthVelocity EastVelocity Displacement Grid
