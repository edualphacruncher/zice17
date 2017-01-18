CC = g++
OPTC = -O3 -fopenmp -ffast-math -mtune=native
LIBTSG = -L../ -ltasmaniansparsegrid
ITSG = -I../
OMPENABLE = -fopenmp

all: example

example: example.cpp
	$(CC) $(OPTC) $(ITSG) $(OMPENABLE) example.cpp -o example $(LIBTSG)
	
clean:
	rm -fr example
