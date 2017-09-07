# Logs Date-Wise

## 10-08-2017: Atlas Problem
* Calculating the info about the atlases, I found that I need to ensure that all the variables in matlab are in **double** format. Otherwise, you'll get the wrong answers. 
* The BrainNetome 2mm atlas gets overlayed on MNI space and they have same dimensions also. It can be checked by fsleyes. 
* If the data is given in 3mm format and while you are having the atlas in 2mm, then get the Transformation matrix by registering the subjects on **MNI152_T1_2mm_brain.nii.gz** with nearest neighbor interpolation ensured and then multiply the ROIs with the transformation matrix in order to get the ROI atlas in 3mm. Otherwise, you can straightaway work with atlas of 3mm.

## 22-08-2017: ATLAS Software problem
I had been given the task to use [Ankita's toolbox](https://github.com/Ank0905/funcConn) for the correlations between the voxels. For compiling it, I required ATLAS software for optimizing the operations and getting the correlations faster. But it didn't work and showed the absence of some .so files. This consumed a lot of time but I later on starting working on the toolkit without the optimized versions. So, I needed OpenBLAS, LAPACK and OpenMP. I also had to made some changes in the MakeFile of the toolkit. The new MakeFile is : 
```
all: fconn.cpp
        g++ -L /usr/lib/OpenBLAS/lib \
        -I /usr/lib/OpenBLAS/include -I . \
        -std=c++11 -Wno-deprecated-declarations  fconn.cpp \
        -o fconn.o -lopenblas -llapack -fopenmp                                                      
```
To test on some functional image, type: ` ./fconn.o -i <INPUT> -o <OUTPUTDIR will be formed> -a 1`.

I also tested its performance by `perf` tool of Linux. For the current input, it showed:
```
     506706.508340      task-clock (msec)         #    1.039 CPUs utilized          
         12,34,059      context-switches          #    0.002 M/sec                  
             4,650      cpu-migrations            #    0.009 K/sec                  
         25,37,734      page-faults               #    0.005 M/sec                  
17,56,70,95,84,390      cycles                    #    3.467 GHz                    
12,06,75,48,96,563      instructions              #    0.69  insn per cycle         
   95,87,57,90,941      branches                  #  189.214 M/sec                  
      29,53,59,245      branch-misses             #    0.31% of all branches        

     487.754769540 seconds time elapsed
```

## 23-08-2017: GFLOPS calculation
The timings that I calculated in the previous iteration also included loading and saving timings, so it was incorrect for calculating the GFLOPS for the program. The right way is :
* Add the clock where the operations are being done in the program i.e. keep on noting down the time for each iteration and total it. 
* In my case, the total number of operations had been calculated this way:
  Let the total number of voxels be **N** and timeseries be **T** while the non-zero voxels are **M**.
  * **Calculating the mean**: (T-1 +1) operations for each voxel. Therefore, *total operations for mean* = T*N
  * **Calculating the standard Deviation**: Formula is *E(X^2) - [E(X)]^2*. Hence, for calculating E(X^2), it takes *T * N +T * N = 2T * N* operations. Also, subtracting and squaring E(X), takes another 2N operations. Therefore, *total operations for standard deviation* = 2T * N + 2N
  * **Normalizing the data**: Now, we need to normalize each and every data point. Hence, *total operations for normalization* = 2T * N
  * **Whole brain to Whole brain matrix multiplication**: Now, since M are the number of voxels with non-zero values. Hence, we shall do matrix multiplication of *M X T with T X M matrix*. Therefore, *total operations for matrix multiplication* = (2T - 1)*M^2 
  * Roughly, we have total number of operations: 5T * N + 2N + (2T - 1) * M^2
  
In my case, the total time it took for whole brain to whole brain correlation was **182.256s** and the number of operations for a brain with N = 271633 , T = 176 and M = 110351 came out to be **4.274 X 10^12**. 

Therefore, 
* **Calculating the mean**: 47807480 ops
* **Calculating the Standard Deviation**: 96158226 ops
* **Normalizing the data**: 95614690 ops
* **Whole brain to Whole brain matrix multiplication**: 4.274 X 10^12 ops
**GFLOPS for this iteration = 23.45 GFLOPs**

while the processor has clock frequency = 3.2 GHz and processor name is Intel (Haswell)Core i7. 

**Rough way for calculating GFLOPS of processor = Number of cores X Flops/Cycle X Frequency.**

For Intel Haswell which I am using, you can read Flops/Cycle [here](https://en.wikipedia.org/wiki/FLOPS) = 16 DP FLOPs/Cycle.

**Rough GFLOPs for my processor = 8 X 16 X 3.2 = 409.6 GFLOPs** or you can download [Linpack](https://serverfault.com/a/822011) from here, a good tool for calculating the GigaFLOPs. From LinPack, the range for GFLOPs for the processor was **145-160 GFLOPs**.
