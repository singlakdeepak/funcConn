import os
import json
import nipype.interfaces.utility as util
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

from nipype.interfaces import fsl
import nipype.workflows.fmri.fsl.preprocess as Preprocessor

def getthreshop(thresh):
    return ['-thr %.10f -Tmin -bin' % (0.1 * val[1]) for val in thresh]

tolist = lambda x: [x]
highpass_lowpass_operand = lambda x: '-bptf %.10f %.10f' % x


def give_Slice_Timer_Node(SliceTimeCorrect,time_repeat):
    if (SliceTimeCorrect ==1):
        slicetimer = MapNode(fsl.SliceTimer(index_dir=True,
                                     interleaved=False,
                                     output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                             iterfield=['in_file'],
                          name="slicetimer")
    elif (SliceTimeCorrect == 2):
        slicetimer = MapNode(fsl.SliceTimer(index_dir=False,
                                     interleaved=False,
                                     output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                              iterfield=['in_file'],
                          name="slicetimer")

    elif (SliceTimeCorrect == 3):
        slicetimer = MapNode(fsl.SliceTimer(index_dir=False,
                                     interleaved=True,
                                     output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                              iterfield=['in_file'],
                          name="slicetimer")
    elif (SliceTimeCorrect == 4):
        slicetimer = MapNode(fsl.SliceTimer(custom_order = '',
                                    output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                              iterfield=['in_file'],
                          name="slicetimer")
    elif (SliceTimeCorrect == 5):
        slicetimer = MapNode(fsl.SliceTimer(custom_timings = '',
                                    output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                              iterfield=['in_file'],
                          name="slicetimer")
    return slicetimer

def create_parallelfeat_preproc(name='featpreproc', highpass= True, 
                                Intensity_Norm = True,
                                BETextract = True,
                                MotionCorrection = 0, 
                                SliceTimeCorrect = 0,
                                time_repeat = 3):
    """
    Preprocess each run with FSL independently of the others
    Parameters
    ----------
    ::
      name : name of workflow (default: featpreproc)
      highpass : boolean (default: True)
    Inputs::
        inputspec.func : functional runs (filename or list of filenames)
        inputspec.fwhm : fwhm for smoothing with SUSAN
        inputspec.highpass : HWHM in TRs (if created with highpass=True)
    Outputs::
        outputspec.reference : volume to which runs are realigned
        outputspec.motion_parameters : motion correction parameters
        outputspec.realigned_files : motion corrected files
        outputspec.motion_plots : plots of motion correction parameters
        outputspec.mask : mask file used to mask the brain
        outputspec.smoothed_files : smoothed functional data
        outputspec.highpassed_files : highpassed functional data (if highpass = True)
        outputspec.normalized_files : Files after Internsity Normalization (if highpass = False)
    Example
    -------
    >>> preproc = create_parallelfeat_preproc()
    >>> preproc.inputs.inputspec.func = ['f3.nii', 'f5.nii']
    >>> preproc.inputs.inputspec.fwhm = 5
    >>> preproc.inputs.inputspec.highpass = 128./(2*2.5)
    >>> preproc.base_dir = '/tmp'
    >>> preproc.run() # doctest: +SKIP
    >>> preproc = create_parallelfeat_preproc(highpass=False)
    >>> preproc.inputs.inputspec.func = 'f3.nii'
    >>> preproc.inputs.inputspec.fwhm = 5
    >>> preproc.base_dir = '/tmp'
    >>> preproc.run() # doctest: +SKIP
    """

    featpreproc = Workflow(name=name)

    """
    Set up a node to define all inputs required for the preprocessing workflow
    """

    if highpass:
        inputnode = Node(interface=util.IdentityInterface(fields=['func',
                                                                     'fwhm',
                                                                     'highpass']),
                                                                name='inputspec')
        outputnode = Node(interface=util.IdentityInterface(fields=['reference',
                                                                  'motion_parameters',
                                                                  'realigned_files',
                                                                  'motion_plots',
                                                                  'mask',
                                                                  'smoothed_files',
                                                                  'highpassed_files',
                                                                  'normalized_files']),
                                                             name='outputspec')
    else:
        inputnode = Node(interface=util.IdentityInterface(fields=['func',
                                                                     'fwhm'
                                                                     ]),
                            name='inputspec')
        outputnode = Node(interface=util.IdentityInterface(fields=['reference',
                                                                  'motion_parameters',
                                                                  'realigned_files',
                                                                  'motion_plots',
                                                                  'mask',
                                                                  'smoothed_files',
                                                                  'normalized_files']),
                                                                 name='outputspec')

    """
    Set up a node to define outputs for the preprocessing workflow
    """

    """
    Convert functional images to float representation. Since there can
    be more than one functional run we use a MapNode to convert each
    run.
    """

    img2float = MapNode(interface=fsl.ImageMaths(out_data_type='float',
                                                 op_string = '',
                                                 suffix='_dtype'),
                           iterfield=['in_file'],
                           name='img2float')
    featpreproc.connect(inputnode, 'func', img2float, 'in_file')

    """
    Extract the first volume of the first run as the reference
    """

    extract_ref = MapNode(interface=fsl.ExtractROI(t_size=1),
                             iterfield=['in_file', 't_min'],
                             name = 'extractref')
    
    
    
    """
    Realign the functional runs to the reference (1st volume of first run)
    """
    
    motion_correct = MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                                  save_plots = True),
                            name='realign',
                            iterfield = ['in_file', 'ref_file'])
    
    """
    Plot the estimated motion parameters
    """

    plot_motion = MapNode(interface=fsl.PlotMotionParams(in_source='fsl'),
                            name='plot_motion',
                            iterfield=['in_file'])
    plot_motion.iterables = ('plot_type', ['rotations', 'translations'])

    """
    Extract the mean volume of the first functional run
    """

    meanfunc = MapNode(interface=fsl.ImageMaths(op_string = '-Tmean',
                                                   suffix='_mean'),
                          iterfield=['in_file'],
                          name='meanfunc')
    """
    Strip the skull from the mean functional to generate a mask
    """

    meanfuncmask = MapNode(interface=fsl.BET(mask = True,
                                             no_output=True,
                                             frac = 0.3),
                              iterfield=['in_file'],
                              name = 'meanfuncmask')        
    """
    Mask the functional runs with the extracted mask
    """

    maskfunc = MapNode(interface=fsl.ImageMaths(suffix='_bet',
                                                   op_string='-mas'),
                          iterfield=['in_file', 'in_file2'],
                          name = 'maskfunc')

    """
    Determine the 2nd and 98th percentile intensities of each functional run
    """

    getthresh = MapNode(interface=fsl.ImageStats(op_string='-p 2 -p 98'),
                           iterfield = ['in_file'],
                           name='getthreshold')        

    """
    Threshold the first run of the functional data at 10% of the 98th percentile
    """

    threshold = MapNode(interface=fsl.ImageMaths(out_data_type='char',
                                                 suffix='_thresh'),
                           iterfield=['in_file', 'op_string'],
                           name='threshold')    

    """
    Determine the median value of the functional runs using the mask
    """

    medianval = MapNode(interface=fsl.ImageStats(op_string='-k %s -p 50'),
                           iterfield = ['in_file', 'mask_file'],
                           name='medianval')

    """
    Dilate the mask
    """
    dilatemask = MapNode(interface=fsl.ImageMaths(suffix='_dil',
                                                  op_string='-dilF'),
                            iterfield=['in_file'],
                            name='dilatemask')

    """
    Mask the motion corrected functional runs with the dilated mask
    """

    maskfunc2 = MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                    op_string='-mas'),
                          iterfield=['in_file', 'in_file2'],
                          name='maskfunc2')
    
    def whether_BETextract_or_not(featpreproc,MotionCorrection = 0, SliceTimeCorrect =0, BETextract = True):

        if ((MotionCorrection ==0)and (SliceTimeCorrect==0)):
            if BETextract:
                featpreproc.connect(img2float, 'out_file', meanfunc, 'in_file')
                featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')
                featpreproc.connect(img2float, 'out_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')
                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                featpreproc.connect(img2float, 'out_file', getthresh, 'in_file')
                featpreproc.connect(img2float, 'out_file', threshold, 'in_file')

        elif ((MotionCorrection ==1)and (SliceTimeCorrect==0)):
            if BETextract:
                featpreproc.connect(motion_correct, 'out_file', meanfunc, 'in_file')
                featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')

                featpreproc.connect(motion_correct, 'out_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')

                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                featpreproc.connect(motion_correct, 'out_file', getthresh, 'in_file')
                featpreproc.connect(motion_correct, 'out_file', threshold, 'in_file')

        else:
            if BETextract:
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', meanfunc, 'in_file')
                featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')
                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', getthresh, 'in_file')
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', threshold, 'in_file')
        return featpreproc 

    """
    Realign the functional runs to the reference (1st volume of first run)
    """
    if (MotionCorrection == 1):
        featpreproc.connect(img2float, 'out_file', extract_ref, 'in_file')
        featpreproc.connect(img2float, ('out_file', Preprocessor.pickmiddle), extract_ref, 't_min')
        featpreproc.connect(extract_ref, 'roi_file', outputnode, 'reference')
        
        featpreproc.connect(img2float, 'out_file', motion_correct, 'in_file')
        featpreproc.connect(extract_ref, 'roi_file', motion_correct, 'ref_file')
        featpreproc.connect(motion_correct, 'par_file', outputnode, 'motion_parameters')
        featpreproc.connect(motion_correct, 'out_file', outputnode, 'realigned_files')

        featpreproc.connect(motion_correct, 'par_file', plot_motion, 'in_file')
        featpreproc.connect(plot_motion, 'out_file', outputnode, 'motion_plots')
        
        
        """
        SliceTimer - correct for slice wise acquisition
        """
        if (SliceTimeCorrect != 0):
            slicetimer = give_Slice_Timer_Node(SliceTimeCorrect,time_repeat)

            featpreproc.connect(motion_correct, 'out_file', slicetimer, 'in_file')
            
            featpreproc = whether_BETextract_or_not(featpreproc,
                                                    MotionCorrection = MotionCorrection, 
                                                    SliceTimeCorrect = SliceTimeCorrect, 
                                                    BETextract = BETextract)
            """
            Define a function to get 10% of the intensity
            """
            featpreproc.connect(getthresh, ('out_stat', getthreshop), threshold, 'op_string')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', medianval, 'in_file')
            featpreproc.connect(threshold, 'out_file', medianval, 'mask_file')
            featpreproc.connect(threshold, 'out_file', dilatemask, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', outputnode, 'mask')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', maskfunc2, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', maskfunc2, 'in_file2')            
            
        else:
            featpreproc = whether_BETextract_or_not(featpreproc,
                                                    MotionCorrection = MotionCorrection, 
                                                    SliceTimeCorrect = SliceTimeCorrect, 
                                                    BETextract = BETextract)         
            """
            Define a function to get 10% of the intensity
            """

            featpreproc.connect(getthresh, ('out_stat', getthreshop), threshold, 'op_string')

            featpreproc.connect(motion_correct, 'out_file', medianval, 'in_file')
            featpreproc.connect(threshold, 'out_file', medianval, 'mask_file')


            featpreproc.connect(threshold, 'out_file', dilatemask, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', outputnode, 'mask')

            featpreproc.connect(motion_correct, 'out_file', maskfunc2, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', maskfunc2, 'in_file2')
    
    else:

        """
        SliceTimer - correct for slice wise acquisition
        """
        if (SliceTimeCorrect != 0):
            slicetimer = give_Slice_Timer_Node(SliceTimeCorrect,time_repeat)
            
            featpreproc.connect(img2float, 'out_file', slicetimer, 'in_file')
            featpreproc = whether_BETextract_or_not(featpreproc,
                                                    MotionCorrection = MotionCorrection, 
                                                    SliceTimeCorrect = SliceTimeCorrect, 
                                                    BETextract = BETextract)
            featpreproc.connect(getthresh, ('out_stat', getthreshop), threshold, 'op_string')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', medianval, 'in_file')
            featpreproc.connect(threshold, 'out_file', medianval, 'mask_file')
            featpreproc.connect(threshold, 'out_file', dilatemask, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', outputnode, 'mask')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', maskfunc2, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', maskfunc2, 'in_file2')
        
        else:

            featpreproc = whether_BETextract_or_not(featpreproc,
                                                    MotionCorrection = MotionCorrection, 
                                                    SliceTimeCorrect = SliceTimeCorrect, 
                                                    BETextract = BETextract)
            """
            Define a function to get 10% of the intensity
            """
            featpreproc.connect(getthresh, ('out_stat', getthreshop), threshold, 'op_string')

            featpreproc.connect(img2float, 'out_file', medianval, 'in_file')
            featpreproc.connect(threshold, 'out_file', medianval, 'mask_file')

            featpreproc.connect(threshold, 'out_file', dilatemask, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', outputnode, 'mask')

            featpreproc.connect(img2float, 'out_file', maskfunc2, 'in_file')
            featpreproc.connect(dilatemask, 'out_file', maskfunc2, 'in_file2')
    


    """
    Smooth each run using SUSAN with the brightness threshold set to 75%
    of the median value for each run and a mask consituting the mean
    functional
    """

    smooth = Preprocessor.create_susan_smooth()
    
    """
    Mask the smoothed data with the dilated mask
    """
    maskfunc3 = MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                    op_string='-mas'),
                          iterfield=['in_file', 'in_file2'],
                          name='maskfunc3')
    concatnode = Node(interface=util.Merge(2),
                         name='concat')
    """
    The following nodes select smooth or unsmoothed data depending on the
    fwhm. This is because SUSAN defaults to smoothing the data with about the
    voxel size of the input data if the fwhm parameter is less than 1/3 of the
    voxel size.
    """
    selectnode = Node(interface=util.Select(),name='select')
    
    """
    Scale the median value of the run is set to 10000
    """

    meanscale = MapNode(interface=fsl.ImageMaths(suffix='_gms'),
                          iterfield=['in_file','op_string'],
                          name='meanscale')
    
    """
    Perform temporal highpass filtering on the data
    """
    highpassfilt = MapNode(interface=fsl.ImageMaths(suffix='_tempfilt'),
                      iterfield=['in_file'],
                      name='highpassfilt')
#     """
#     Generate a mean functional image from the first run
#     """
#     meanfunc3 = MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
#                                                     suffix='_mean'),
#                            iterfield=['in_file'],
#                           name='meanfunc3')    
    

    featpreproc.connect(inputnode, 'fwhm', smooth, 'inputnode.fwhm')
    featpreproc.connect(maskfunc2, 'out_file', smooth, 'inputnode.in_files')
    featpreproc.connect(dilatemask, 'out_file', smooth, 'inputnode.mask_file')


    featpreproc.connect(smooth, 'outputnode.smoothed_files', maskfunc3, 'in_file')

    featpreproc.connect(dilatemask, 'out_file', maskfunc3, 'in_file2')

    featpreproc.connect(maskfunc2,('out_file', tolist), concatnode, 'in1')
    featpreproc.connect(maskfunc3,('out_file', tolist), concatnode, 'in2')

    featpreproc.connect(concatnode, 'out', selectnode, 'inlist')

    featpreproc.connect(inputnode, ('fwhm', Preprocessor.chooseindex), selectnode, 'index')
    featpreproc.connect(selectnode, 'out', outputnode, 'smoothed_files')

    datasink = Node(interface=DataSink(), name="datasink")

    if Intensity_Norm:
        featpreproc.connect(selectnode, 'out', meanscale, 'in_file')
        """
        Define a function to get the scaling factor for intensity normalization
        """
        featpreproc.connect(medianval, ('out_stat', Preprocessor.getmeanscale), meanscale, 'op_string')
        featpreproc.connect(meanscale, 'out_file', outputnode, 'normalized_files')
        
        if highpass:

            featpreproc.connect(inputnode, ('highpass', highpass_lowpass_operand), highpassfilt, 'op_string')
            featpreproc.connect(meanscale, 'out_file', highpassfilt, 'in_file')
            featpreproc.connect(highpassfilt, 'out_file', outputnode, 'highpassed_files')
            featpreproc.connect(outputnode, 'highpassed_files',
                      datasink, 'out_file')
        else:
            featpreproc.connect(outputnode, 'normalized_files', datasink, 'out_file')

    else :
        if highpass:
            featpreproc.connect(inputnode, ('highpass', highpass_lowpass_operand), highpassfilt, 'op_string')
            featpreproc.connect(selectnode, 'out', highpassfilt, 'in_file')
            featpreproc.connect(highpassfilt, 'out_file', outputnode, 'highpassed_files')
            featpreproc.connect(outputnode, 'highpassed_files',
                      datasink, 'out_file')
        else:
          featpreproc.connect(outputnode, 'smoothed_files', datasink, 'out_file')

#     if highpass:
#         featpreproc.connect(highpass, 'out_file', meanfunc3, 'in_file')
#     else:
#         featpreproc.connect(meanscale, 'out_file', meanfunc3, 'in_file')

#    featpreproc.connect(meanfunc3, 'out_file', outputnode, 'mean')
    return featpreproc

def reg_workflow(name = 'registration'):
    """Create a FEAT preprocessing workflow
    Parameters
    ----------
    ::
        name : name of workflow (default: 'registration')
    Inputs::
        inputspec.source_files : files (filename or list of filenames to register)
        inputspec.mean_image : reference image to use
        inputspec.anatomical_image : anatomical image to coregister to
        inputspec.target_image : registration target
    Outputs::
        outputspec.func2anat_transform : FLIRT transform
        outputspec.anat2target_transform : FLIRT+FNIRT transform
        outputspec.transformed_files : transformed files in target space
        outputspec.transformed_mean : mean image in target space
    Example
    -------
    """

    register = Workflow(name=name)

    inputnode = Node(interface=util.IdentityInterface(fields=['source_files',
                                                                 'anatomical_image',
                                                                 'target_image']),
                        name='inputspec')
    outputnode = Node(interface=util.IdentityInterface(fields=['func2anat_transform',
                                                                  'anat2target_transform',
                                                                  'transformed_files',
                                                                  'transformed_mean',
                                                                  ]),
                         name='outputspec')
    meanfunc = Node(fsl.ImageMaths(op_string='-Tmean',suffix='_mean'), name = 'meanfunc')
    register.connect(inputnode, 'source_files', meanfunc, 'in_file')
    """
    Estimate the tissue classes from the anatomical image. But use spm's segment
    as FSL appears to be breaking.
    """

    stripper = Node(fsl.BET(), name='stripper')
    register.connect(inputnode, 'anatomical_image', stripper, 'in_file')
    fast = Node(fsl.FAST(), name='fast')
    register.connect(stripper, 'out_file', fast, 'in_files')

    """
    Binarize the segmentation
    """

    binarize = Node(fsl.ImageMaths(op_string='-nan -thr 0.5 -bin'),
                       name='binarize')
    pickindex = lambda x, i: x[i]
    register.connect(fast, ('partial_volume_files', pickindex, 2),
                     binarize, 'in_file')

    """
    Calculate rigid transform from mean image to anatomical image
    """

    mean2anat = Node(fsl.FLIRT(), name='mean2anat')
    mean2anat.inputs.dof = 6
    register.connect(meanfunc,'out_file', mean2anat, 'in_file')
    register.connect(stripper, 'out_file', mean2anat, 'reference')

    """
    Now use bbr cost function to improve the transform
    """

    mean2anatbbr = Node(fsl.FLIRT(), name='mean2anatbbr')
    mean2anatbbr.inputs.dof = 6
    mean2anatbbr.inputs.cost = 'bbr'
    mean2anatbbr.inputs.schedule = os.path.join(os.getenv('FSLDIR'),
                                                'etc/flirtsch/bbr.sch')
    register.connect(meanfunc,'out_file', mean2anatbbr, 'in_file')
    register.connect(binarize, 'out_file', mean2anatbbr, 'wm_seg')
    register.connect(inputnode, 'anatomical_image', mean2anatbbr, 'reference')
    register.connect(mean2anat, 'out_matrix_file',
                     mean2anatbbr, 'in_matrix_file')

    """
    Calculate affine transform from anatomical to target
    """

    anat2target_affine = Node(fsl.FLIRT(), name='anat2target_linear')
    anat2target_affine.inputs.searchr_x = [-180, 180]
    anat2target_affine.inputs.searchr_y = [-180, 180]
    anat2target_affine.inputs.searchr_z = [-180, 180]
    register.connect(stripper, 'out_file', anat2target_affine, 'in_file')
    register.connect(inputnode, 'target_image',
                     anat2target_affine, 'reference')


#     register.connect(anat2target_affine, 'out_matrix_file',
#                      anat2target_nonlinear, 'affine_file')
#     register.connect(inputnode, 'anatomical_image',
#                      anat2target_nonlinear, 'in_file')
#     register.connect(inputnode, 'config_file',
#                      anat2target_nonlinear, 'config_file')
#     register.connect(inputnode, 'target_image',
#                      anat2target_nonlinear, 'ref_file')

    """
    Transform the mean image. First to anatomical and then to target
    """

    warpmean = Node(fsl.ApplyWarp(interp='spline'), name='warpmean')
    datasink = Node(interface=DataSink(), name="datasink")


    register.connect(meanfunc,'out_file', warpmean, 'in_file')
    register.connect(mean2anatbbr, 'out_matrix_file', warpmean, 'premat')
    register.connect(inputnode, 'target_image', warpmean, 'ref_file')
#     register.connect(anat2target_nonlinear, 'fieldcoeff_file',
#                      warpmean, 'field_file')
    register.connect(anat2target_affine, 'out_matrix_file',
                     warpmean, 'postmat')
    
    """
    Transform the remaining images. First to anatomical and then to target
    """

    warpall = MapNode(fsl.ApplyWarp(interp='spline'),
                         iterfield=['in_file'],
                         nested=True,
                         name='warpall')
    register.connect(inputnode, 'source_files', warpall, 'in_file')
    register.connect(mean2anatbbr, 'out_matrix_file', warpall, 'premat')
    register.connect(inputnode, 'target_image', warpall, 'ref_file')
#     register.connect(anat2target_nonlinear, 'fieldcoeff_file',
#                      warpall, 'field_file')
    register.connect(anat2target_affine, 'out_matrix_file',
                     warpall, 'postmat')

    """
    Assign all the output files
    """

    register.connect(warpmean, 'out_file', outputnode, 'transformed_mean')
    register.connect(warpall, 'out_file', outputnode, 'transformed_files')
    register.connect(mean2anatbbr, 'out_matrix_file',
                     outputnode, 'func2anat_transform')
#     register.connect(anat2target_nonlinear, 'fieldcoeff_file',
#                      outputnode, 'anat2target_transform')
    register.connect(anat2target_affine, 'out_matrix_file',
                     outputnode, 'anat2target_transform')
    register.connect(warpall, 'out_file', datasink,'out_file')

    return register