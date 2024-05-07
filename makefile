all: STExa TD
STExa: STExa.o
	g++ -std=c++14 -fopenmp -O3 STExa.o -o STExa
STExa.o: STExa.cpp
	g++ -std=c++14 -fopenmp -O3 -c  STExa.cpp
TD: TD.o
	g++ -std=c++14 -fopenmp -O3 TD.o -o TD
TD.o: TD.cpp 
	g++ -std=c++14 -fopenmp -O3 -c  TD.cpp
