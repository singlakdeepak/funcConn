all: fconn.cpp
	g++ -L /usr/lib/OpenBLAS/lib \
	-I /usr/lib/OpenBLAS/include -I . \
	-std=c++11 -Wno-deprecated-declarations  fconn.cpp \
	-o fconn.o -lopenblas -llapack -fopenmp
