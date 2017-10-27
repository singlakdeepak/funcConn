import nibabel as nib
import numpy as np
import os
from scipy import stats


def div0( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def calc_mean_and_std(ROICorrMaps, n_subjects):
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
    n_subjects = len(ROICorrMaps)
    Sample_mean_Array, Sample_std_Array = calc_mean_and_std(ROICorrMaps, n_subjects)
    ttest_1samp_for_all = div0((Sample_mean_Array - PopMean) * np.sqrt(n_subjects), Sample_std_Array)
    pval = stats.t.sf(np.abs(ttest_1samp_for_all), n_subjects-1)*2

# def ttest_1samp_for_one_ROI(ROICorrMaps, PopMean = 0.0, ROINo = 0):
    