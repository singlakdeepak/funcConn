import nibabel as nib
import numpy as np
import os
from scipy import stats


def div0( a, b ):
	'''
	It is meant for ignoring the places where standard deviation 
	is zero.
	'''
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def calc_mean_and_std(ROICorrMaps, n_subjects):
	'''
	Function for calculating the mean and standard 
	deviation of the data. At a time, only one of the nii
	file is loaded and the elements keep on adding as we 
	enumerate over the subjects.
	'''
	if (n_subjects != 0):
        f = nib.load(ROICorrMaps[0])
        dimensions = f.get_header().get_data_shape()
    else:
        break
    
    Sample_mean_Array = np.zeros(dimensions)
    Sample_std_Array = np.zeros(dimensions)
    for count, subject in enumerate(n_subjects):
        Corr_data = nib.load(subject).get_data()
        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
    Sample_mean_Array /= n_subjects
    Sample_std_Array = np.sqrt(Sample_std_Array/n_subjects - np.square(Sample_mean_Array))
    return Sample_mean_Array,Sample_std_Array

def ttest_1samp_for_all_ROIs(ROICorrMaps, PopMean = 0.0):
    '''
    This is the 1 sample t-test for ROI correlation maps.
    df = no of subjects - 1
    * ROICorrMaps is the list of filepaths of ROI correlation 
    maps for a group.
    * Each ROI correlation map has the 4th dimension equal to 
    the number of ROIs.
    * It calculates both the ttest as well as the p values.
    '''
    n_subjects = len(ROICorrMaps)
    Sample_mean_Array, Sample_std_Array = calc_mean_and_std(ROICorrMaps, n_subjects)
    ttest_1samp_for_all = div0((Sample_mean_Array - PopMean) * np.sqrt(n_subjects), Sample_std_Array)
    pval = stats.t.sf(np.abs(ttest_1samp_for_all), n_subjects-1)*2
    return ttest_1samp_for_all, pval

# def ttest_1samp_for_one_ROI(ROICorrMaps, PopMean = 0.0, ROINo = 0):
    