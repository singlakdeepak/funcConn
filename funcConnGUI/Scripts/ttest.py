import nibabel as nib
import numpy as np
import os
from scipy import stats
from numpy import ma
import scipy.special as special
from statsmodels.stats import multitest
from multiprocessing import Pool
from functools import partial
import multiprocessing.managers
class MyManager(multiprocessing.managers.BaseManager):
    pass
MyManager.register('ma_empty_like', ma.empty_like, multiprocessing.managers.ArrayProxy)

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

def convert_pvals_to_log_fmt(pvalues,Sample_mean_ArrayA = None, Sample_mean_ArrayB = None):
    if (Sample_mean_ArrayA and Sample_mean_ArrayB) is not None:
        return (-1*np.log10(pvalues)*np.sign(Sample_mean_ArrayA - Sample_mean_ArrayB))
    return (-1*np.log10(pvalues))

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
        print('Done subject ', count +1)                                                                                                                                                                                                       
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
        print('Done subject ', count + 1)                                                                                                                                                                                               
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
    if save_pval_in_log_fmt:
        ttest, pvals = _ttest_1samp(Sample_mean_Array,
                        Sample_std_Array,
                        n_subjects,
                        PopMean = PopMean)
        return ttest, convert_pvals_to_log_fmt(pvals)

    return _ttest_1samp(Sample_mean_Array,
                        Sample_std_Array,
                        n_subjects,
                        PopMean = PopMean)


def ttest_1samp_ROIs_if_npy(ROICorrMaps,
                            PopMean = 0.0,
                            save_pval_in_log_fmt = True,
                            applyFisher = False):
    n_subjects = len(ROICorrMaps)
    assert (n_subjects>0)
    Sample_mean_Array, Sample_std_Array = \
                                calc_mean_and_std_if_npy( ROICorrMaps,
                                                        n_subjects, ddof =1,
                                                        applyFisher = applyFisher)
    if save_pval_in_log_fmt:
        ttest, pvals = _ttest_1samp(Sample_mean_Array,
                        Sample_std_Array,
                        n_subjects,
                        PopMean = PopMean)
        return ttest, convert_pvals_to_log_fmt(pvals)

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
                        save_pval_in_log_fmt = True,
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
    if save_pval_in_log_fmt:
        ttest, pvals = _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)
        return ttest, convert_pvals_to_log_fmt(pvals, 
                                            Sample_mean_ArrayA = Sample_mean_ArrayA,
                                            Sample_mean_ArrayB = Sample_mean_ArrayB)
    # pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def ttest_ind_samples_if_npy(ROICorrMapsA, ROICorrMapsB, 
                            save_pval_in_log_fmt = True,
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
    if save_pval_in_log_fmt:
        ttest, pvals = _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)
        return ttest, convert_pvals_to_log_fmt(pvals,
                                            Sample_mean_ArrayA = Sample_mean_ArrayA,
                                            Sample_mean_ArrayB = Sample_mean_ArrayB)
    return  Sample_mean_ArrayA, \
                Sample_mean_ArrayB, \
                _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def convert_ma_to_np(MaskedArrayObj):
    return ma.filled(MaskedArrayObj)

def fdrcorrect_worker2(A,B,args):
    return fdrcorrect_worker(A,B,*args)

def fdrcorrect_worker(fdrcorrected_brain, 
                    rejected_pvals, 
                    roi_number, pvals, 
                    is_npy, 
                    alpha=0.05, type = 'indep_samps',
                    method='indep'):
    '''
    Taken andd modified from : https://github.com/statsmodels/statsmodels/blob/master/statsmodels/stats/multitest.py
    http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html

    pvalue correction for false discovery rate
    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
    Parameters
    ----------
    pvals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate
    method : {'indep', 'negcorr')
    Returns
    -------
    rejected : array, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : array
        pvalues adjusted for multiple hypothesis testing to limit FDR
    Notes
    -----
    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to alpha * m/m_0 where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)
    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).
    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
    fdr_by.
    '''

    def _ecdf(x):
        '''no frills empirical cdf used in fdrcorrection
        '''
        nobs = len(x)
        return np.arange(1,nobs+1)/float(nobs)
    mask_pvals = pvals.mask
    indices = np.where(mask_pvals ==False)
    Truepvals = ma.compressed(pvals)

    pvals_sortind = np.argsort(Truepvals)
    pvals_sorted = np.take(Truepvals, pvals_sortind)

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
##    elif method in ['n', 'negcorr']:
##        cm = np.sum(np.arange(len(pvals)))
##        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject
    if not is_npy:
        rejected_pvalBrain = ma.masked_all_like(pvals)
        rejected_pvalBrain[indices] = reject_
        if (type == 'indep_samps'):
            rejected_pvals[:,:,:,roi_number] = rejected_pvalBrain
        else:
            rejected_pvals = rejected_pvalBrain

        pvals_correctedBrain = ma.masked_all_like(pvals)
        pvals_correctedBrain[indices] = pvals_corrected_
        if (type =='indep_samps'):
            fdrcorrected_brain[:,:,:,roi_number] = pvals_correctedBrain
        else:
            fdrcorrected_brain = pvals_correctedBrain
        print("Done for ", roi_number)
    else:
        if (type == 'indep_samps'):
            rejected_pvals[roi_number,indices] = reject_
            fdrcorrected_brain[roi_number,indices] = pvals_corrected_
        else:
            rejected_pvals[indices] = reject_
            fdrcorrected_brain[indices] = pvals_corrected_
        print("Done for ", roi_number)

def fdr_correction(pvalues , type = 'indep_samps', is_npy = False):
    '''
    pvalues :: Pvalue maps for all ROIs
    Two types:
    indep_samps: When the ROIs are taken independently and the FDR is done considering the 
           the tests only in that ROI. 
    all: When all the tests are treated as one.  
    '''
    FDR_types = ['indep_samps', 'all']

    if (type == 'indep_samps'):
        procs = 8
        pool = Pool(procs)

        m = MyManager()
        m.start()
        fdrcorrected_brain = m.ma_empty_like(pvalues)
        rejected_pvals = m.ma_empty_like(pvalues)
        func = partial(fdrcorrect_worker2, fdrcorrected_brain, rejected_pvals)
        # no_rois : Total ROIS in the P-value file
        if not is_npy:
            no_rois = pvalues.shape[3]
        else:
            no_rois = pvalues.shape[0]
        print('Total no of ROIs ',no_rois)
        MaxPools = no_rois//procs
        print('MaxPools: ', MaxPools)
        # for roi_number in range(0,MaxPools*procs,procs):
        #     pool_inputs = [] #np.arange(number_of_ROIs)
        #     select_roi = 0
        #     while (select_roi<procs):
        #         print(roi_number +select_roi)
        #         if not is_npy:
        #             pool_inputs.append((roi_number+select_roi, 
        #                             pvalues[:,:,:,roi_number+select_roi], is_npy))
        #         else:
        #             pool_inputs.append((roi_number+select_roi,
        #                             ma.masked_array(pvalues[roi_number + select_roi,:],
        #                             fill_value = 0),
        #                             is_npy))
        #         select_roi+=1

        #     data_outputs = pool.map(func, pool_inputs)

        # if (no_rois%procs!=0):
        #     pool_inputs = [] #np.arange(number_of_ROIs)
        #     for roi_number in range(MaxPools*procs,no_rois):
        #         if not is_npy:
        #             pool_inputs.append((roi_number, 
        #                                 pvalues[:,:,:,roi_number], 
        #                                 is_npy))
        #         else:
        #             pool_inputs.append((roi_number,
        #                                 ma.masked_array(pvalues[roi_number,:],
        #                                     fill_value = 0),
        #                                 is_npy))

        #     data_outputs = pool.map(func, pool_inputs)

        pool_inputs = [] #np.arange(number_of_ROIs)
        for roi_number in range(no_rois):

            print(roi_number)
            if not is_npy:
                pool_inputs.append((roi_number, 
                                pvalues[:,:,:,roi_number], is_npy))
            else:
                pool_inputs.append((roi_number, pvalues[roi_number, :], is_npy))
                                # ma.masked_array(pvalues[roi_number ,:],
                                # fill_value = 0),
                                # is_npy))


        data_outputs = pool.map(func, pool_inputs)

        return rejected_pvals, fdrcorrected_brain
    elif (type == 'all'):
        fdrcorrected_brain = []
        rejected_pvals = []
        fdrcorrect_worker(fdrcorrected_brain, rejected_pvals, 
                            None, pvalues,is_npy, 
                            type = type)
        return rejected_pvals, fdrcorrected_brain
    else:
        return 1


