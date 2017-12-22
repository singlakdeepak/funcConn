OpenBLAS: fconn.cpp
	g++ -L /usr/lib/OpenBLAS/lib \
	-I /usr/lib/OpenBLAS/include -I . \
	-std=c++11 -Wno-deprecated-declarations  fconn.cpp \
	-o fconn.o -lopenblas -llapack -fopenmp -D BUILD1
MKL: fconn.cpp
	g++ fconn.cpp -o fconn.o -fopenmp -std=c++11 \
	-L /mnt/project1/home1/cs5130287/intel/mkl/lib/intel64 \
	-I /mnt/project1/home1/cs5130287/intel/mkl/include/ -I . \
	-L /mnt/project1/home1/cs5130287/intel/lib/intel64 -lmkl_blas95_lp64 \
	-lmkl_intel_lp64 -lmkl_lapack95_lp64 -liomp5 -lm -dl -lmkl_intel_thread -lmkl_core -lpthread
CBLAS: fconn.cpp
	 g++ -L /usr/local/atlas -I . -std=c++11 -Wdeprecated-declarations fconn.cpp \
	-o fconn.o -fopenmp -llapack -llapack -lcblas -D BUILD1

