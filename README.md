# FCONN : Functional Connectivity Toolkit

## ABOUT
FCONN is a multi processing FMRI toolkit for functional connectivity analysis. It aims for finding the all pair correlation 
* between all the voxels of the brain,
* between a given voxel and all the voxels of the brain,and
* between a given ROI and all the voxels of the brain.

## PRE-REQUISITES
Before instaling this one must have the following softwares downloaded and configured on their PC :-
* GCC Compiler 5.0
* Automatically Tuned Linear Algebra Software (ATLAS) [It must be configured on the pc] *ATLAS can be downloaded from http://math-atlas.sourceforge.net/*
* QTCreator v5.7
* FSL Software

## COMPILATION
Use the make command to compile the toolkit

## RUNNING THE TOOLKIT 

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





