'''
This script has been used for making synthetic data. You 
can use it if you want to make such brain ROI maps. These 
maps don't tell anything and have been made using the 
random function of Python Numpy distribution.
'''

import nibabel as nib
import numpy as np
import os
import sys

np.random.seed(0)

start_name_for_file = int(sys.argv[1])
end_name_for_file = int(sys.argv[2])
for i in range(start_name_for_file,end_name_for_file+1):
	Subject = np.random.randn(91,109,91,246)
	Refer = nib.load('abc.nii.gz')
	header = Refer.header
	header['dim'][4] = 246
	affine = Refer.affine
	array_img = nib.Nifti1Image(Subject, affine, header)
	nib.save(array_img,  'file%d.nii.gz'%i)
