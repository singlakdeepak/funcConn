import nibabel as nib
import numpy as np
import os
from scipy import stats
from numpy import ma
import scipy.special as special

def div0( a, b ):
    '''
	It is meant for ignoring the places where standard deviation 
	is zero.
	'''
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.divide( a, b )
        # c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def calc_mean_and_std(ROICorrMaps, n_subjects, mask, applyFisher = False):
    '''
	Function for calculating the mean and standard 
	deviation of the data. At a time, only one of the nii
	file is loaded and the elements keep on adding as we 
	enumerate over the subjects.
	'''
    if (n_subjects != 0):
        f = nib.load(ROICorrMaps[0])
        dimensions = f.get_header().get_data_shape()
        print(dimensions)
    else:
        exit
    print(ROICorrMaps)
    Sample_mean_Array = np.zeros(dimensions)
    Sample_std_Array = np.zeros(dimensions)
    Sample_mean_Array = ma.masked_array(Sample_mean_Array, mask = mask)
    Sample_std_Array = ma.masked_array(Sample_std_Array, mask = mask)
    for count, subject in enumerate(ROICorrMaps):

        Corr_data = nib.load(subject).get_data()
        Corr_data = ma.masked_array(Corr_data, mask = mask)
        if applyFisher:
        	Corr_data = np.arctanh(Corr_data)
        
        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
        print('Done subject ', count)                                                                                                                                                                                                       
    Sample_mean_Array /= n_subjects
    Sample_std_Array = np.sqrt(Sample_std_Array/n_subjects - np.square(Sample_mean_Array))
    return Sample_mean_Array,Sample_std_Array

def ttest_1samp_for_all_ROIs(ROICorrMaps, ROIAtlas, PopMean = 0.0, applyFisher = False):
    '''
    This is the 1 sample t-test for ROI correlation maps.
    df = no of subjects - 1
    * ROICorrMaps is the list of filepaths of ROI correlation 
    maps for a group.
    * Each ROI correlation map has the 4th dimension equal to 
    the number of ROIs.
    * It calculates both the ttest as well as the p values.
    
    QUESTIONS???????????????????????????????????????????????
    For application of the Fisher transform, I saw that it is 
    same as the inverse hyperbolic tangent function. 
    Doubt is regarding the standard deviation of the distribution after
    applying Fisher. It was written that the sd is now 1/sqrt(no_of_subjs - 3).
    So, that means for each voxel or variable, the sd now becomes this.

    Ref: https://docs.scipy.org/doc/numpy/reference/generated/numpy.arctanh.html
     https://en.wikipedia.org/wiki/Fisher_transformation
    TO BE ASKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # The ttest will return t value Inf or NaN where the denom is
    # zero. See what to return in these places. Ask tomorrow.
    '''
    mask = nib.load(ROIAtlas).get_data()
    totalROIs = np.amax(mask)
    mask = ma.masked_object(mask,0).mask
    totalROIs = 246
    mask = np.repeat(mask[:, :, :, np.newaxis], totalROIs, axis=3)
    print(mask.shape)
    n_subjects = len(ROICorrMaps)
    Sample_mean_Array, Sample_std_Array = calc_mean_and_std(ROICorrMaps, 
    	                                                    n_subjects,
                                                            mask, 
    	                                                    applyFisher = applyFisher)
    ttest_1samp_for_all = div0((Sample_mean_Array - PopMean) * np.sqrt(n_subjects), Sample_std_Array)
    
    df = n_subjects - 1
    # pval = stats.t.sf(np.abs(ttest_1samp_for_all), df)*2
    pval = special.betainc(0.5*df, 0.5, df/(df + ttest_1samp_for_all*ttest_1samp_for_all)).reshape(ttest_1samp_for_all.shape)
    return ttest_1samp_for_all, pval




def ttest_ind_samples(ROICorrMapsA, ROICorrMapsB, ROIAtlas, equal_var = True, applyFisher = False):
    '''
    Modified from https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html ,
    https://github.com/scipy/scipy/blob/v0.19.1/scipy/stats/stats.py#L3950-L4072

    Since it didn't support if the data is large and everything can't be loaded at once. So,
    such modification has been made.
    '''
    mask = nib.load(ROIAtlas).get_data()
    totalROIs = np.amax(mask)
    mask = ma.masked_object(mask,0).mask
    mask = np.repeat(mask[:, :, :, np.newaxis], totalROIs, axis=3)

    n_subjectsA = len(ROICorrMapsA)
    Sample_mean_ArrayA, Sample_std_ArrayA = calc_mean_and_std(ROICorrMapsA, 
    	                                                      n_subjectsA,
                                                              mask, 
    	                                                      applyFisher = applyFisher)
    Sample_var_ArrayA = np.square(Sample_std_ArrayA)
    del(Sample_std_ArrayA)

    n_subjectsB = len(ROICorrMapsB)
    Sample_mean_ArrayB, Sample_std_ArrayB = calc_mean_and_std(ROICorrMapsB, 
    	                                                      n_subjectsB,
                                                              mask, 
    	                                                      applyFisher = applyFisher)
    Sample_var_ArrayB = np.square(Sample_std_ArrayB)
    del(Sample_std_ArrayB)
    
    if equal_var:
        # force df to be an array for masked division not to throw a warning
        df = ma.asanyarray(n_subjectsA + n_subjectsB - 2.0)
        svar = ((n_subjectsA-1)*Sample_var_ArrayA+(n_subjectsB-1)*Sample_var_ArrayB) / df
        denom = ma.sqrt(svar*(1.0/n_subjectsA + 1.0/n_subjectsB))  # n-D computation here!
    else:
        vn1 = Sample_var_ArrayA/n_subjectsA
        vn2 = Sample_var_ArrayB/n_subjectsB
        with np.errstate(divide='ignore', invalid='ignore'):
            df = (vn1 + vn2)**2 / (vn1**2 / (n_subjectsA - 1) + vn2**2 / (n_subjectsB - 1))

        # If df is undefined, variances are zero.
        # It doesn't matter what df is as long as it is not NaN.
        df = np.where(np.isnan(df), 1, df)
        denom = ma.sqrt(vn1 + vn2)

    with np.errstate(divide='ignore', invalid='ignore'):
        ttest_ind = (Sample_mean_ArrayA - Sample_mean_ArrayB) / denom
    pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return ttest_ind , pvalues




# def ttest_1samp_for_one_ROI(ROICorrMaps, PopMean = 0.0, ROINo = 0):
    