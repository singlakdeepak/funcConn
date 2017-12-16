# FCONN : A Toolkit on analysis of Functional Connectivity in BOLD fMRI

## ABOUT
FCONN is a multi processing fMRI toolkit for functional connectivity analysis. It aims for finding the all pair correlation 
* between all the voxels of the brain,
* between a given voxel and all the voxels of the brain,and
* between a given ROI and all the voxels of the brain.

Further, it is aimed for within and between the group analysis as well as giving out the results in proper formats. Since the task takes a lot of time, so all the efforts have been made to fasten it by parallelizing the large sets of calculations. The optimization preocedures have been done by parallelizing the computation on the cores and using Intel-MKL (optimized for CBLAS) and Python libraries (NiPype, Multiprocessing). 

It also comes with a GUI interface(made in Qt-Creator) for giving the user freedom and neatness. All the options have been listed in the GUI in the sequential order and the outputs at various stages can also be saved. The speed of the preprocessing analysis has been compared with the FSL and it comes out be approximately 8 times faster.

## PRE-REQUISITES
Before instaling this one must have the following softwares downloaded and configured on their PC :-

* GCC Compiler 5.0
* Either [Automatically Tuned Linear Algebra Software (ATLAS)](http://math-atlas.sourceforge.net/) 
    or Intel - MKL* (For increasing the speed, * optional *)
* QTCreator v5.7
* FSL Software
* [Anaconda](https://www.anaconda.com/download/#linux) complete for Python 3.6
* Nipype (See the instructions [here](http://miykael.github.io/nipype-beginner-s-guide/installation.html))

    *Read [this small blog](https://github.com/singlakdeepak/Wiki/blob/master/InstallingMKL.md) for installation of Intel-MKL.


## COMPILATION
Use the make command to compile the toolkit

## RUNNING THE FCONN TOOLKIT 

Usage : fconn -i <fmri.nii/nii.gz> 
```  
Compulsory arguments(You have to specify the following):-
	 -i	path of the input volume
	 -o	project name(will create a directory with same name with the results in it )
Only one of the arguments must be present
  	 -r 	path of the volume containg the desired ROI
	-s 1	for seed to all voxel mode
		 -x 	x-coordinate for seed (compulosry in -s mode)
		 -y	y-coordinate for seed (compulosry in -s mode)
		 -z	z-coordinate for seed (compulosry in -s mode)
	-a 1	for all voxels to all voxels mode
Optional arguments(You may optionally specify the following)
	-t 	an upper threshold
	-h	display the help message
```
## UNDERSTANDING THE OUTPUT
The toolkit will create a directory named the same as that of the name given as in -o option. Inside the directory there will be:
* ROI directory: Which will contain all the roi to all voxels correlations
* Voxel directory: Which will contain all the voxels to all voxels correlations
* Voxelmap file: Which contains all the voxel to file mappings

Ankita's makefile:
```
 all: fconn.cpp
g++ -L /usr/local/atlas -I . -std=c++11 -Wdeprecated-declarations fconn.cpp -o fconn.o -fopenmp -llapack -llapack -lcblas
```


