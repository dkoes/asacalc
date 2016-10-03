

all: asacalc.cpp
	g++ -o asacalc asacalc.cpp -O3 -I/usr/local/include/openbabel-2.0 -lopenbabel -lann -lboost_system
