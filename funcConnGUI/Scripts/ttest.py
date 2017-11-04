import nibabel as nib
import numpy as np
import os
from scipy import stats
from numpy import ma
import scipy.special as special
from statsmodels.stats import multitest

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

def calc_mean_and_std(ROICorrMaps, n_subjects, ROIAtlasmask, ddof =1, applyFisher = False):
    '''
	Function for calculating the mean and standard 
	deviation of the data. At a time, only one of the nii
	file is loaded and the elements keep on adding as we 
	enumerate over the subjects.
	'''
    mask = nib.load(ROIAtlasmask).get_data()
    mask = ma.masked_object(mask,0).mask
    if (n_subjects != 0):
        f = nib.load(ROICorrMaps[0])
        dimensions = f.get_header().get_data_shape()
        print(dimensions)
    else:
        exit
    mask = np.repeat(mask[:, :, :, np.newaxis], dimensions[3], axis=3)
    print(ROICorrMaps)

    Sample_mean_Array = np.zeros(dimensions)
    Sample_std_Array = np.zeros(dimensions)
    Sample_mean_Array = ma.masked_array(Sample_mean_Array, 
                                       mask = mask, 
                                        fill_value = 0)
    Sample_std_Array = ma.masked_array(Sample_std_Array, 
                                      mask = mask , 
                                       fill_value = 0)
    for count, subject in enumerate(ROICorrMaps):

        Corr_data = nib.load(subject).get_data()
        Corr_data = ma.masked_array(Corr_data, mask = mask)
        if applyFisher:
            Corr_data = np.arctanh(Corr_data)
        
        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
        print('Done subject ', count)                                                                                                                                                                                                       
    Sample_mean_Array /= n_subjects
    Sample_std_Array = np.sqrt((Sample_std_Array - n_subjects*np.square(Sample_mean_Array))/(n_subjects - ddof))
    return Sample_mean_Array,Sample_std_Array

def calc_mean_and_std_if_npy(ROICorrMaps, n_subjects, ddof =1, applyFisher = False):
    '''
    Function to be used if the file is given in the format 
    No of ROIs versus All brain voxels in the ROI mapped.
    '''
    print(ROICorrMaps)
    initialize = np.load(ROICorrMaps[0])
    initialize = ma.masked_array(initialize)
    if applyFisher:
        initialize = np.arctanh(initialize)
    Sample_mean_Array = ma.masked_array(initialize, 
                                        fill_value = 0)
    Sample_std_Array = ma.masked_array(np.square(initialize), 
                                       fill_value = 0)
    del initialize
    print('Done subject ', 0)
    for count, subject in enumerate(ROICorrMaps[1:]):

        Corr_data = np.load(subject)
        Corr_data = ma.masked_array(Corr_data) 
        if applyFisher:
            Corr_data = np.arctanh(Corr_data)
        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
        print('Done subject ', count)                                                                                                                                                                                               
    Sample_mean_Array /= n_subjects
    Sample_std_Array = np.sqrt((Sample_std_Array - n_subjects*np.square(Sample_mean_Array))/(n_subjects - ddof))
    return Sample_mean_Array,Sample_std_Array

def _ttest_1samp(Sample_mean_Array, Sample_std_Array, n_subjects, PopMean = 0.0):
    ttest_1samp_for_all = div0((Sample_mean_Array - PopMean) \
                            * np.sqrt(n_subjects), Sample_std_Array)
    df = n_subjects - 1
    # pval = stats.t.sf(np.abs(ttest_1samp_for_all), df)*2
    pval = special.betainc(0.5*df, 0.5, df/ \
            (df + ttest_1samp_for_all*ttest_1samp_for_all)).reshape(ttest_1samp_for_all.shape)
    # ttest_1samp_for_all, pval = ma.filled(ttest_1samp_for_all), ma.filled(pval)

    return ttest_1samp_for_all, pval


def ttest_1samp_for_all_ROIs(ROICorrMaps, 
                                ROIAtlasmask, 
                                PopMean = 0.0, 
                                applyFisher = False):
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


    n_subjects = len(ROICorrMaps)
    assert (n_subjects>0)
    Sample_mean_Array, Sample_std_Array = calc_mean_and_std(ROICorrMaps, 
                                                            n_subjects,
                                                            ROIAtlasmask, ddof =1,
                                                            applyFisher = applyFisher)
    ttest_1samp_for_all, pval = _ttest_1samp(Sample_mean_Array, 
                                             Sample_std_Array,
                                             n_subjects,
                                             PopMean = PopMean)
    return ttest_1samp_for_all, pval


def ttest_1samp_ROIs_if_npy(ROICorrMaps,
                            PopMean = 0.0,
                            applyFisher = False):
    n_subjects = len(ROICorrMaps)
    assert (n_subjects>0)
    Sample_mean_Array, Sample_std_Array = \
                                calc_mean_and_std_if_npy( ROICorrMaps,
                                                        n_subjects, ddof =1,
                                                        applyFisher = applyFisher)
    return _ttest_1samp(Sample_mean_Array,
                        Sample_std_Array,
                        n_subjects,
                        PopMean = PopMean)


def _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB,Sample_var_ArrayB, n_subjectsB,
                equal_var = True):
    if equal_var:
        # force df to be an array for masked division not to throw a warning
        df = ma.asanyarray(n_subjectsA + n_subjectsB - 2.0)
        svar = ((n_subjectsA-1)*Sample_var_ArrayA+(n_subjectsB-1)*Sample_var_ArrayB)/ df
        denom = ma.sqrt(svar*(1.0/n_subjectsA + 1.0/n_subjectsB))  # n-D computation here!
    else:
        vn1 = Sample_var_ArrayA/n_subjectsA
        vn2 = Sample_var_ArrayB/n_subjectsB
        df = (vn1 + vn2)**2 / (vn1**2 / (n_subjectsA - 1) + vn2**2 / (n_subjectsB - 1))

        # If df is undefined, variances are zero.
        # It doesn't matter what df is as long as it is not NaN.
        df = np.where(np.isnan(df), 1, df)
        denom = ma.sqrt(vn1 + vn2)

    with np.errstate(divide='ignore', invalid='ignore'):
        ttest_ind = (Sample_mean_ArrayA - Sample_mean_ArrayB) / denom
    pvalues = special.betainc(0.5*df, 0.5, df/(df + ttest_ind*ttest_ind)).reshape(ttest_ind.shape)
    
    # ttest_ind, pvalues = ma.filled(ttest_ind), ma.filled(pvalues)
    return ttest_ind, pvalues

def ttest_ind_samples(ROICorrMapsA, ROICorrMapsB, ROIAtlasmask, 
                    equal_var = True, applyFisher = False):
    '''
    Modified from https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html ,
    https://github.com/scipy/scipy/blob/v0.19.1/scipy/stats/stats.py#L3950-L4072

    Since it didn't support if the data is large and everything can't be loaded at once. So,
    such modification has been made.
    '''

    n_subjectsA = len(ROICorrMapsA)
    n_subjectsB = len(ROICorrMapsB)
    assert (n_subjectsA > 0)
    assert (n_subjectsB > 0)
    Sample_mean_ArrayA, Sample_std_ArrayA = calc_mean_and_std(ROICorrMapsA, 
                                                              n_subjectsA,
                                                              ROIAtlasmask, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayA = np.square(Sample_std_ArrayA)
    del(Sample_std_ArrayA)

    # n_subjectsB = len(ROICorrMapsB)
    Sample_mean_ArrayB, Sample_std_ArrayB = calc_mean_and_std(ROICorrMapsB, 
                                                              n_subjectsB,
                                                              ROIAtlasmask, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayB = np.square(Sample_std_ArrayB)
    del(Sample_std_ArrayB)

    # pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def ttest_ind_samples_if_npy(ROICorrMapsA, ROICorrMapsB, equal_var = True, applyFisher = False):
    '''
    Modified from https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html ,
    https://github.com/scipy/scipy/blob/v0.19.1/scipy/stats/stats.py#L3950-L4072

    Since it didn't support if the data is large and everything can't be loaded at once. So,
    such modification has been made.
    '''

    n_subjectsA = len(ROICorrMapsA)
    n_subjectsB = len(ROICorrMapsB)
    assert (n_subjectsA > 0)
    assert (n_subjectsB > 0)
    Sample_mean_ArrayA, Sample_std_ArrayA = calc_mean_and_std_if_npy(ROICorrMapsA, 
                                                              n_subjectsA, ddof = 1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayA = np.square(Sample_std_ArrayA)
    del(Sample_std_ArrayA)
    Sample_mean_ArrayB, Sample_std_ArrayB = calc_mean_and_std_if_npy(ROICorrMapsB, 
                                                              n_subjectsB, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayB = np.square(Sample_std_ArrayB)
    del(Sample_std_ArrayB)
    # pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def convert_ma_to_np(MaskedArrayObj):
    return ma.filled(MaskedArrayObj)

def fdr_correction(pvalues , type = 'ind_ROIs'):
    '''
    Two types:
    ind_ROIs: When the ROIs are taken independently and the FDR is done considering the 
           the tests only in that ROI. 
    all: When all the tests are treated as one.  
    '''

# def ttest_1samp_for_one_ROI(ROICorrMaps, PopMean = 0.0, ROINo = 0):
    
