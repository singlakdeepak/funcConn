all: fconn.cpp
	g++ -L /usr/local/atlas -I . -std=c++11 -Wdeprecated-declarations fconn.cpp -o fconn.o -fopenmp -llapack -llapack -lcblas
