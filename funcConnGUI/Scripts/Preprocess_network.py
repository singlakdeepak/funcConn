import os
import json
import nipype.interfaces.utility as util
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

from nipype.interfaces import fsl

def getthreshop(thresh):
    return ['-thr %.10f -Tmin -bin' % (0.1 * val[1]) for val in thresh]

tolist = lambda x: [x]
highpass_lowpass_operand = lambda x: '-bptf %.10f %.10f' % x

def pickmiddle(files):
    from nibabel import load
    import numpy as np
    from nipype.utils import NUMPY_MMAP
    middlevol = []
    for f in files:
        middlevol.append(int(np.ceil(load(f, mmap=NUMPY_MMAP).shape[3] / 2)))
    return middlevol

def pickfirst(files):
    if isinstance(files, list):
        return files[0]
    else:
        return files

def getbtthresh(medianvals):
    return [0.75 * val for val in medianvals]

def getusans(x):
    return [[tuple([val[0], 0.75 * val[1]])] for val in x]

def chooseindex(fwhm):
    if fwhm < 1:
        return [0]
    else:
        return [1]

def global_sig_regression(in_file, mask_file):
    '''
    Refer the paper: The Global Signal and Observed Anticorrelated Resting State Brain Networks: 
                      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694109/
    Also see, Towards a consensus regarding global signal regression for 
              resting state functional connectivity MRI: 
              https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5489207/
    Also see, https://fcp-indi.github.io/docs/user/nuisance.html
    '''

    import nibabel as nib
    import numpy as np
    from os.path import join as opj
    import os

    sub_id = in_file.split('/')[-1].split('.')[0]
    glb_reg_file_name = sub_id + '_glb_reg.nii.gz'
    glb_reg_file_name = opj(os.getcwd(),glb_reg_file_name)
    if os.path.exists(glb_reg_file_name):
        print('Saved file in : %s'%glb_reg_file_name)
        
    else:
        brain_data = nib.load(in_file)
        brain = brain_data.get_data()
        brain_affine = brain_data.affine
        x_dim, y_dim, z_dim, num_volumes = brain.shape
        
        mask_Obj = nib.load(mask_file)
        mask_data = mask_Obj.get_data()

        brain_voxels_X,brain_voxels_Y,brain_voxels_Z = np.where(mask_data==1)
        num_brain_voxels = len(brain_voxels_X)

        voxel_matrix = np.zeros((num_volumes, num_brain_voxels))
        # Fill up the voxel_matrix 

        for i in range(num_volumes):
            voxel_matrix[i] = brain[brain_voxels_X,
                                    brain_voxels_Y,
                                    brain_voxels_Z,
                                    i]

        mean_signal = np.expand_dims(np.mean(voxel_matrix, axis=1),axis=1)
        pseudo_inv_mean_sig = np.dot(np.linalg.inv(np.dot(mean_signal.T,mean_signal)),mean_signal.T)
             
        beta_g = np.dot(pseudo_inv_mean_sig, voxel_matrix)
        voxel_matrix -= np.dot(mean_signal, beta_g)
     
        for i in range(num_volumes):
            brain[brain_voxels_X,brain_voxels_Y,brain_voxels_Z,i] = voxel_matrix[i]
        brain_obj = nib.Nifti1Image(brain, affine = brain_affine)
        nib.save(brain_obj,glb_reg_file_name)

        print('Saved GSR file in: ', glb_reg_file_name)
    return glb_reg_file_name        

def getmeanscale(medianvals):
    return ['-mul %.10f' % (10000. / val) for val in medianvals]


def create_susan_smooth(name="susan_smooth", separate_masks=True):
    """Create a SUSAN smoothing workflow
    Parameters
    ----------
    ::
        name : name of workflow (default: susan_smooth)
        separate_masks : separate masks for each run
    Inputs::
        inputnode.in_files : functional runs (filename or list of filenames)
        inputnode.fwhm : fwhm for smoothing with SUSAN
        inputnode.mask_file : mask used for estimating SUSAN thresholds (but not for smoothing)
    Outputs::
        outputnode.smoothed_files : functional runs (filename or list of filenames)
    Example
    -------
    >>> smooth = create_susan_smooth()
    >>> smooth.inputs.inputnode.in_files = 'f3.nii'
    >>> smooth.inputs.inputnode.fwhm = 5
    >>> smooth.inputs.inputnode.mask_file = 'mask.nii'
    >>> smooth.run() # doctest: +SKIP
    """

    susan_smooth = Workflow(name=name)

    """
    Set up a node to define all inputs required for the preprocessing workflow
    """

    inputnode = Node(interface=util.IdentityInterface(fields=['in_files',
                                                                 'fwhm',
                                                                 'mask_file',
                                                                   'median_value']),
                        name='inputnode')

    """
    Smooth each run using SUSAN with the brightness threshold set to 75%
    of the median value for each run and a mask consituting the mean
    functional
    """

    smooth = MapNode(interface=fsl.SUSAN(),
                        iterfield=['in_file', 'brightness_threshold', 'usans'],
                        name='smooth')

    """
    Mask the motion corrected functional runs with the dilated mask
    """

    if separate_masks:
        mask = MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                   op_string='-mas'),
                          iterfield=['in_file', 'in_file2'],
                          name='mask')
    else:
        mask = MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                   op_string='-mas'),
                          iterfield=['in_file'],
                          name='mask')
    susan_smooth.connect(inputnode, 'in_files', mask, 'in_file')
    susan_smooth.connect(inputnode, 'mask_file', mask, 'in_file2')

    """
    Determine the mean image from each functional run
    """

    meanfunc = MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
                                                   suffix='_mean'),
                          iterfield=['in_file'],
                          name='meanfunc2')
    susan_smooth.connect(mask, 'out_file', meanfunc, 'in_file')

    """
    Merge the median values with the mean functional images into a coupled list
    """

    merge = Node(interface=util.Merge(2, axis='hstack'),
                    name='merge')
    susan_smooth.connect(meanfunc, 'out_file', merge, 'in1')
    susan_smooth.connect(inputnode, 'median_value', merge, 'in2')

    """
    Define a function to get the brightness threshold for SUSAN
    """
    susan_smooth.connect(inputnode, 'fwhm', smooth, 'fwhm')
    susan_smooth.connect(inputnode, 'in_files', smooth, 'in_file')
    susan_smooth.connect(inputnode, ('median_value', getbtthresh), smooth, 'brightness_threshold')
    susan_smooth.connect(merge, ('out', getusans), smooth, 'usans')

    outputnode = Node(interface=util.IdentityInterface(fields=['smoothed_files']),
                         name='outputnode')

    susan_smooth.connect(smooth, 'smoothed_file', outputnode, 'smoothed_files')

    return susan_smooth

def give_Slice_Timer_Node(SliceTimeCorrect,time_repeat):
    '''
    It gives the type of slicetimer to be used in the script according to SliceTimeCorrect
    Parameters
    ----------
    Inputs::
    Slicetimcorrect:
    # 1 : From top to down
    # 2 : From bottom to up
    # 3 : Interleaved
    # 4 : Custom order file is given 
    # 5 : Custom timings file is given

    time_repeat : repetition time in seconds
    
    Outputs::
    slicetimer: Slice Timer Node
    '''
    if (SliceTimeCorrect ==1):
        slicetimer = MapNode(fsl.SliceTimer(index_dir=False,
                                     interleaved=False,
                                     output_type='NIFTI_GZ',
                                     time_repetition=time_repeat),
                             iterfield=['in_file'],
                          name="slicetimer")
    elif (SliceTimeCorrect == 2):
        slicetimer = MapNode(fsl.SliceTimer(index_dir=True,
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
                                BETvalue = 0,
                                robustBET = False,
                                GSR = False,
                                MotionCorrection = 0, 
                                SliceTimeCorrect = 0,
                                time_repeat = 3):
    """
    Preprocess each run with FSL independently of the others. It has some modifications. 
    In case we don't want to use some nodes, those are specified to be false.

    Parameters
    ----------
    ::
      name : name of workflow (default: featpreproc)
      highpass : boolean (default: True)
    Inputs::
        inputspec.func : functional runs (filename or list of filenames)
        inputspec.fwhm : fwhm for smoothing with SUSAN
        inputspec.highpass : It is a tuple (Highpass sigma, LowPass sigma). To be passed only when
                            highpass is true.
        highpass : whether Highpass and/or lowpass is to be done or not.
        Intensity_Norm : whether Intensity Normalization is to be done or not.
        BETextract : whether BETextract is to be done or not.
        MotionCorrection: #0 : No motion correction, #1: MCFlirt
        SliceTimeCorrect: 
              #0 : No Slice Time Correction
              # 1 : From top to down
              # 2 : From bottom to up
              # 3 : Interleaved
              # 4 : Custom order file is given 
              # 5 : Custom timings file is given
        time_repeat: Repetion time of the volume or scan.

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
    >>> preproc.inputs.inputspec.highpass = (12,15)
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
    Realign the functional runs to the reference (1st volume of first run).
    The MCFLIRT uses the mean image for aligning all the other volumes to it.
    """
    
    motion_correct = MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                                  save_plots = True,
                                                  interpolation = 'spline'),
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
    Extract the mean volume of the first functional run.
    This extract the mean volume of the functional run. 
    """

    meanfunc = MapNode(interface=fsl.ImageMaths(op_string = '-Tmean',
                                                   suffix='_mean'),
                          iterfield=['in_file'],
                          name='meanfunc')
    """
    Strip the skull from the mean functional to generate a mask
    """
    if robustBET:
        meanfuncmask = MapNode(interface=fsl.BET(args = '-R',
                                              mask = True,
                                             no_output=True,
                                             frac = BETvalue),
                              iterfield=['in_file'],
                              name = 'meanfuncmask')
    else:
        meanfuncmask = MapNode(interface=fsl.BET(mask = True,
                                             no_output=True,
                                             frac = BETvalue),
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
    Threshold the first run of the functional data at 10% of the 98th percentile.

    """

    threshold = MapNode(interface=fsl.ImageMaths(out_data_type='char',
                                                 suffix='_thresh'),
                           iterfield=['in_file', 'op_string'],
                           name='threshold')    

    """
    Determine the median value of the functional runs using the mask.
    The median values means having p 50.
    """

    medianval = MapNode(interface=fsl.ImageStats(op_string='-k %s -p 50'),
                           iterfield = ['in_file', 'mask_file'],
                           name='medianval')

    """
    Dilate the mask. It is used for SUSAN smoothing since dilating the mask 
    makes the smoothing go better.
    """
    dilatemask = MapNode(interface=fsl.ImageMaths(suffix='_dil',
                                                  op_string='-dilF'),
                            iterfield=['in_file'],
                            name='dilatemask')

    """
    Mask the motion corrected functional runs with the dilated mask.
    Dilate mask is also applied on all the volumes and in case the output of 
    SUSAN smoothing is not chosen, it comes out to be the final answer. 
    """

    maskfunc2 = MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                    op_string='-mas'),
                          iterfield=['in_file', 'in_file2'],
                          name='maskfunc2')
    """
    Smooth each run using SUSAN with the brightness threshold set to 75%
    of the median value for each run and a mask consituting the mean
    functional
    """

    smooth = create_susan_smooth()
    
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
    Scale the median value of the run is set to 10000.
    This node is for Intensity normalization. 
    """

    meanscale = MapNode(interface=fsl.ImageMaths(suffix='_gms'),
                          iterfield=['in_file','op_string'],
                          name='meanscale')
    Intensity_Normalize = MapNode(interface=fsl.ImageMaths(op_string='-inm 10000',
                                                    suffix='_mean'),
                           iterfield=['in_file'],
                          name='intensity_norm')    

    """
    Perform temporal highpass filtering on the data
    """
    highpassfilt = MapNode(interface=fsl.ImageMaths(suffix='_tempfilt'),
                      iterfield=['in_file'],
                      name='highpassfilt')
    """
    Generate a mean functional image from the first run
    """
    meanfunc3 = MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
                                                    suffix='_mean'),
                           iterfield=['in_file'],
                          name='meanfunc3')



    add_node = MapNode(interface = fsl.ImageMaths(op_string = '-add',
                                                    suffix = 'highpassfilt_correct'),
                            iterfield = ['in_file', 'in_file2'],
                            name = 'highpassfilt_correct')

    if GSR:
        GSR_node = MapNode(util.Function(function=global_sig_regression, 
                                    input_names=['in_file','mask_file'],
                                    output_names=['glb_reg_file_name']),
                          iterfield=['in_file','mask_file'],
                          name = 'gsr')

    def whether_BETextract_or_not(featpreproc,MotionCorrection = 0, SliceTimeCorrect =0, BETextract = True):
        '''
        It selects whether BET is to be done or not and in according to that attaches the nodes. 
        In case the Motion Correction is not done, then MCFLIRT node is not attached. In case the Slice Time 
        Correction is not required, then it is not attached in the pipeline. 
        '''
        if ((MotionCorrection ==0)and (SliceTimeCorrect==0)):
            featpreproc.connect(img2float, 'out_file', meanfunc, 'in_file')
            featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')
            if BETextract:
                featpreproc.connect(img2float, 'out_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')
                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                featpreproc.connect(img2float, 'out_file', getthresh, 'in_file')
                featpreproc.connect(img2float, 'out_file', threshold, 'in_file')

        elif ((MotionCorrection ==0)and (SliceTimeCorrect!=0)):
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', meanfunc, 'in_file')
            featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')
            if BETextract:
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')
                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', getthresh, 'in_file')
                featpreproc.connect(slicetimer, 'slice_time_corrected_file', threshold, 'in_file')
        else:
            featpreproc.connect(motion_correct, 'out_file', meanfunc, 'in_file')
            featpreproc.connect(meanfunc, 'out_file', meanfuncmask, 'in_file')
            if BETextract:
                featpreproc.connect(motion_correct, 'out_file', maskfunc, 'in_file')
                featpreproc.connect(meanfuncmask, 'mask_file', maskfunc, 'in_file2')

                featpreproc.connect(maskfunc, 'out_file', getthresh, 'in_file')
                featpreproc.connect(maskfunc, 'out_file', threshold, 'in_file')
            else:
                
                featpreproc.connect(motion_correct, 'out_file', getthresh, 'in_file')
                featpreproc.connect(motion_correct, 'out_file', threshold, 'in_file')            

        return featpreproc 

    """
    Realign the functional runs to the reference (1st volume of first run)
    """
    '''
    Also, in case the Motion Correction is done, it becomes an input to the latter nodes. If the 
    slice timing is done then it is used as the input. 
    '''

    if (MotionCorrection == 1):
        if (SliceTimeCorrect != 0): 
            slicetimer = give_Slice_Timer_Node(SliceTimeCorrect,time_repeat)
            featpreproc.connect(img2float, 'out_file', slicetimer, 'in_file')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', extract_ref, 'in_file')
            featpreproc.connect(slicetimer, ('slice_time_corrected_file', pickmiddle), extract_ref, 't_min')
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', motion_correct, 'in_file')
        else:
            featpreproc.connect(img2float, 'out_file', extract_ref, 'in_file')
            featpreproc.connect(img2float, ('out_file', pickmiddle), extract_ref, 't_min')
            featpreproc.connect(img2float, 'out_file', motion_correct, 'in_file')
        
        featpreproc.connect(extract_ref, 'roi_file', outputnode, 'reference')
        featpreproc.connect(extract_ref, 'roi_file', motion_correct, 'ref_file')
        featpreproc.connect(motion_correct, 'par_file', outputnode, 'motion_parameters')
        featpreproc.connect(motion_correct, 'out_file', outputnode, 'realigned_files')

        featpreproc.connect(motion_correct, 'par_file', plot_motion, 'in_file')
        featpreproc.connect(plot_motion, 'out_file', outputnode, 'motion_plots')
        
        
        """
        SliceTimer - correct for slice wise acquisition
        """
            
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
        featpreproc.connect(motion_correct, 'out_file', smooth, 'inputnode.in_files')

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
            featpreproc.connect(slicetimer, 'slice_time_corrected_file', smooth, 'inputnode.in_files')
        
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
            featpreproc.connect(img2float, 'out_file', smooth, 'inputnode.in_files')    


   
    

    featpreproc.connect(inputnode, 'fwhm', smooth, 'inputnode.fwhm')
    featpreproc.connect(medianval, 'out_stat', smooth, 'inputnode.median_value')
    featpreproc.connect(dilatemask, 'out_file', smooth, 'inputnode.mask_file')


    featpreproc.connect(smooth, 'outputnode.smoothed_files', maskfunc3, 'in_file')

    featpreproc.connect(dilatemask, 'out_file', maskfunc3, 'in_file2')

    featpreproc.connect(maskfunc2,('out_file', tolist), concatnode, 'in1')
    featpreproc.connect(maskfunc3,('out_file', tolist), concatnode, 'in2')

    featpreproc.connect(concatnode, 'out', selectnode, 'inlist')

    featpreproc.connect(inputnode, ('fwhm', chooseindex), selectnode, 'index')
    featpreproc.connect(selectnode, 'out', outputnode, 'smoothed_files')

    datasink = Node(interface=DataSink(), name="datasink")

    if Intensity_Norm:
        featpreproc.connect(selectnode,'out', Intensity_Normalize, 'in_file')
        # To change Normalization : featpreproc.connect( inputnode, 
        #                           ('intensity_norm',operand), Intensity_Normalize,'op_string')
        """
        Define a function to get the scaling factor for intensity normalization
        """
        featpreproc.connect(Intensity_Normalize, 'out_file', outputnode, 'normalized_files')

        if highpass:

            featpreproc.connect(inputnode, ('highpass', highpass_lowpass_operand), highpassfilt, 'op_string')
            featpreproc.connect(Intensity_Normalize, 'out_file', highpassfilt, 'in_file')
            featpreproc.connect(Intensity_Normalize, 'out_file', meanfunc3, 'in_file')

            featpreproc.connect(highpassfilt, 'out_file', add_node, 'in_file')
            featpreproc.connect(meanfunc3, 'out_file', add_node, 'in_file2')
            if GSR:
                featpreproc.connect(add_node,'out_file', GSR_node, 'in_file')
                featpreproc.connect(meanfuncmask , 'mask_file', GSR_node, 'mask_file')
                featpreproc.connect(GSR_node, 'glb_reg_file_name', outputnode,'highpassed_files')
            else:
                featpreproc.connect(add_node, 'out_file', outputnode, 'highpassed_files')
            featpreproc.connect(outputnode, 'highpassed_files',
                      datasink, 'out_file')
        else:
            if GSR:
                featpreproc.connect(outputnode,'normalized_files', GSR_node, 'in_file')
                featpreproc.connect(meanfuncmask , 'mask_file', GSR_node, 'mask_file')
                featpreproc.connect(GSR_node,'glb_reg_file_name', datasink, 'out_file')
            else:
                featpreproc.connect(outputnode, 'normalized_files', datasink, 'out_file')

    else :

        featpreproc.connect(selectnode, 'out', meanscale, 'in_file')

        featpreproc.connect(medianval, ('out_stat', getmeanscale), meanscale, 'op_string')
        featpreproc.connect(meanscale, 'out_file', outputnode, 'normalized_files')
        
        if highpass:
            
            featpreproc.connect(meanscale, 'out_file', meanfunc3, 'in_file')
            featpreproc.connect(inputnode, ('highpass', highpass_lowpass_operand), highpassfilt, 'op_string')
            featpreproc.connect(meanscale, 'out_file', highpassfilt, 'in_file')
            featpreproc.connect(highpassfilt, 'out_file', add_node, 'in_file')
            featpreproc.connect(meanfunc3, 'out_file', add_node, 'in_file2')

            if GSR:
                featpreproc.connect(add_node,'out_file', GSR_node, 'in_file')
                featpreproc.connect(meanfuncmask , 'mask_file', GSR_node, 'mask_file')
                featpreproc.connect(GSR_node, 'glb_reg_file_name', outputnode,'highpassed_files')
            else:
                featpreproc.connect(add_node, 'out_file', outputnode, 'highpassed_files')
            featpreproc.connect(outputnode, 'highpassed_files',
                      datasink, 'out_file')
        else:
            if GSR:
                featpreproc.connect(outputnode,'normalized_files', GSR_node, 'in_file')
                featpreproc.connect(meanfuncmask , 'mask_file', GSR_node, 'mask_file')
                featpreproc.connect(GSR_node,'glb_reg_file_name', datasink, 'out_file')
            else:
                featpreproc.connect(outputnode, 'normalized_files', datasink, 'out_file')

    return featpreproc

def ROI_transformation(name = 'registration', applyGSR = False):
    register = Workflow(name=name)
    inputnode = Node(name = 'inputspec', interface = util.IdentityInterface(fields=['source_files','mask_files','ROI_File','func2std']))
    meanfunc = MapNode(interface=fsl.ExtractROI(t_size=1),
                             iterfield=['in_file', 't_min'],
                             name = 'meanfunc')
    inv_mat = MapNode(fsl.ConvertXFM(invert_xfm = True), 
                      iterfield = ['in_file'], 
                      name = 'inv_mat')
    transform_ROI = MapNode(fsl.ApplyXFM(interp='nearestneighbour'), 
                      iterfield = ['in_matrix_file','reference'],
                      name='transform_ROI')    
    if applyGSR:
        GSR_node = MapNode(util.Function(function=global_sig_regression, 
                                    input_names=['in_file','mask_file'],
                                    output_names=['glb_reg_file_name']),
                          iterfield=['in_file','mask_file'],
                          name = 'gsr')
        datasink_GSR_applied_ips = Node(interface=DataSink(), name="datasink_GSR_applied_ips")
        outputnode = Node(name = 'outputspec', interface= util.IdentityInterface(fields = ['transformed_ROI','GSR_applied_inputs']))
        register.connect(inputnode,'source_files', GSR_node, 'in_file')
        register.connect(inputnode, 'mask_files', GSR_node, 'mask_file')
        register.connect(GSR_node, 'glb_reg_file_name', meanfunc,'in_file')
    else:
        outputnode = Node(name = 'outputspec', interface= util.IdentityInterface(fields = ['transformed_ROI']))
        register.connect(inputnode, 'source_files', meanfunc, 'in_file')
    
    register.connect(inputnode, ('source_files', pickmiddle), meanfunc, 't_min')
    transform_ROI.inputs.padding_size = 0
    register.connect(inputnode,'func2std', inv_mat, 'in_file')
    register.connect(inputnode, 'ROI_File', transform_ROI, 'in_file')
    register.connect(inv_mat, 'out_file', transform_ROI, 'in_matrix_file')
    register.connect(meanfunc, 'roi_file', transform_ROI, 'reference')

    datasink_transformedROI = Node(interface=DataSink(), name="datasink_transformedROI")    

    """
    Assign all the output files
    """
    register.connect(transform_ROI, 'out_file', outputnode, 'transformed_ROI')
    register.connect(outputnode, 'transformed_ROI', datasink_transformedROI, 'out_file')
    if applyGSR:
        register.connect(GSR_node, 'glb_reg_file_name', outputnode, 'GSR_applied_inputs')
        register.connect(outputnode, 'GSR_applied_inputs', datasink_GSR_applied_ips, 'out_file')    
    return register



def reg_workflow_with_Anat(no_subjects,  name = 'registration'):
    """Create a FEAT preprocessing workflow
    Parameters
    ----------
    ::
        name : name of workflow (default: 'registration')
    Inputs::
        inputspec.source_files : files (filename or list of filenames to register)
        inputspec.anatomical_images : anatomical images in the same subject-wise order for coregistration.
        inputspec.target_image : registration target
    Outputs::
        outputspec.func2anat_transform : FLIRT transform
        outputspec.transformed_files : transformed files in target space
    Example
    -------
    """

    register = Workflow(name=name)
         
    inputnode = Node(interface=util.IdentityInterface(fields=['source_files',
                                                                 'anatomical_images',
                                                                 'target_image',
                                                                 'ROI_File']),
                        name='inputspec')
    outputnode = Node(interface=util.IdentityInterface(fields=['func2anat_transform',
                                                                'func2std_transform',
                                                                'transformed_ROI'
                                                                  ]),
                         name='outputspec')

    ''' 
    Pipeline:
    #1 : Calculate the mean image from the functional run.
    #2 : Do BET on the mean image.
    #3 : FAST is done for segmenting the White matter.
    #4 : Binarize is done for making the mask out of the file generated after the
         FAST segmentation.
    #5 : Now calculate the mean 2 anamtomical co-registration matrix. 
    #6 : BBR is done for refining the matrix generated from mean2anatomical co-reg.
    #7 : Now, calculate the transformation from anatomical to reference file. 
    #8 : Calculate the mask from the transformed anat2target file.
    #9 : Now , apply transformation using the target_image as reference and 
         mean2anat matrix as transformation matrix. Keep Spline interpolation.
    #10 : Now mask the output with the anat2targetmask so that the values outside 
          the brain go zero. 

    '''

    # meanfunc = MapNode(fsl.ImageMaths(op_string='-Tmean',suffix='_mean'), iterfield = ['in_file'],name = 'meanfunc')
    meanfunc = MapNode(interface=fsl.ExtractROI(t_size=1),
                             iterfield=['in_file', 't_min'],
                             name = 'meanfunc')
    stripper = MapNode(fsl.BET(frac = 0.3), 
      iterfield = ['in_file'], name='stripper')
    fast = MapNode(fsl.FAST(), 
          iterfield =['in_files'], 
          name='fast')
    selectfile = MapNode(interface=util.Select(index=[2]), 
                        iterfield = ['inlist'],
                        name='select')

    
    binarize = MapNode(fsl.ImageMaths(op_string='-nan -thr 0.5 -bin'),
                   iterfield = ['in_file'],
                   name='binarize')
    mean2anat = MapNode(fsl.FLIRT(interp = 'trilinear'), 
                        iterfield = ['in_file','reference'], 
                        name='mean2anat')
    mean2anatbbr = MapNode(fsl.FLIRT(interp = 'trilinear'), 
                    iterfield = ['in_file','wm_seg','reference','in_matrix_file'], 
                    name='mean2anatbbr')
    anat2target_affine = MapNode(fsl.FLIRT(), 
                        iterfield = ['in_file'], 
                        name='anat2target_linear')
    concat_mat = MapNode(fsl.ConvertXFM(concat_xfm = True), 
                        iterfield=['in_file','in_file2'], 
                        name='concat_mat')
    inv_mat = MapNode(fsl.ConvertXFM(invert_xfm = True), 
                      iterfield = ['in_file'], 
                      name = 'inv_mat')
    transform_ROI = MapNode(fsl.ApplyXFM(interp='nearestneighbour'), 
                      iterfield = ['in_matrix_file','reference'],
                      name='transform_ROI')

    register.connect(inputnode, 'source_files', meanfunc, 'in_file')
    register.connect(inputnode, ('source_files', pickmiddle), meanfunc, 't_min')
    """
    Estimate the tissue classes from the anatomical image. But use spm's segment
    as FSL appears to be breaking.
    """
    
    register.connect(inputnode, 'anatomical_images', stripper, 'in_file')

    register.connect(stripper, 'out_file', fast, 'in_files')
    """
    Binarize the segmentation
    """
    pickindex = lambda x, i: x[i]
    register.connect(fast, 'partial_volume_files', selectfile, 'inlist')

    register.connect(selectfile, 'out', binarize, 'in_file')
    # register.connect(fast, ('partial_volume_files', pickindex, 2),
                     # binarize, 'in_file')

    """
    Calculate rigid transform from mean image to anatomical image
    """
    mean2anat.inputs.dof = 12
    mean2anat.inputs.searchr_x = [-180, 180]
    mean2anat.inputs.searchr_y = [-180, 180]
    mean2anat.inputs.searchr_z = [-180, 180]
    register.connect(meanfunc,'roi_file', mean2anat, 'in_file')
    register.connect(stripper, 'out_file', mean2anat, 'reference')

    """
    Now use bbr cost function to improve the transform
    """

    mean2anatbbr.inputs.dof = 12
    mean2anatbbr.inputs.searchr_x = [-180, 180]
    mean2anatbbr.inputs.searchr_y = [-180, 180]
    mean2anatbbr.inputs.searchr_z = [-180, 180]
    mean2anatbbr.inputs.cost = 'bbr'
    mean2anatbbr.inputs.schedule = os.path.join(os.getenv('FSLDIR'),
                                                'etc/flirtsch/bbr.sch')
    register.connect(meanfunc,'roi_file', mean2anatbbr, 'in_file')
    register.connect(binarize, 'out_file', mean2anatbbr, 'wm_seg')
    # register.connect(inputnode, 'anatomical_images', mean2anatbbr, 'reference')
    register.connect(stripper, 'out_file', mean2anatbbr, 'reference')
    register.connect(mean2anat, 'out_matrix_file',
                     mean2anatbbr, 'in_matrix_file')
    """
    Calculate affine transform from anatomical to target
    """

    anat2target_affine.inputs.searchr_x = [-180, 180]
    anat2target_affine.inputs.searchr_y = [-180, 180]
    anat2target_affine.inputs.searchr_z = [-180, 180]
    anat2target_affine.inputs.cost = 'corratio'
    register.connect(stripper, 'out_file', anat2target_affine, 'in_file')
    register.connect(inputnode, 'target_image',
                     anat2target_affine, 'reference')
    # register.connect(anat2target_affine,'out_file', anat2targetmask, 'in_file')

    transform_ROI.inputs.padding_size = 0

    register.connect(mean2anatbbr, 'out_matrix_file', concat_mat, 'in_file')
    register.connect(anat2target_affine, 'out_matrix_file', concat_mat,'in_file2')
    register.connect(concat_mat,'out_file', inv_mat, 'in_file')
    register.connect(inputnode, 'ROI_File', transform_ROI, 'in_file')
    register.connect(inv_mat, 'out_file', transform_ROI, 'in_matrix_file')
    register.connect(meanfunc, 'roi_file', transform_ROI, 'reference')

    datasink_transformedROI = Node(interface=DataSink(), name="datasink_transformedROI")
    datasink_transformedROI.inputs.remove_dest_dir = True
    datasink_func2std = Node(interface=DataSink(), name="datasink_func2std")
    datasink_func2std.inputs.remove_dest_dir = True
    """
    Assign all the output files
    """
    register.connect(mean2anatbbr, 'out_matrix_file',
                     outputnode, 'func2anat_transform')
    register.connect(concat_mat, 'out_file', 
                      outputnode, 'func2std_transform')
    register.connect(transform_ROI, 'out_file', outputnode, 'transformed_ROI')
    register.connect(outputnode, 'transformed_ROI', datasink_transformedROI, 'out_file')
    register.connect(outputnode, 'func2std_transform', datasink_func2std,'out_file')

    return register

def reg_workflow_wo_Anat(no_subjects, name = 'registration'):
    """Create a FEAT preprocessing workflow
    Parameters
    ----------
    ::
        name : name of workflow (default: 'registration')
    Inputs::
        inputspec.source_files : files (filename or list of filenames to register)
        inputspec.anatomical_images : anatomical images in the same subject-wise order for coregistration.
        inputspec.target_image : registration target
    Outputs::
        outputspec.func2anat_transform : FLIRT transform
        outputspec.transformed_files : transformed files in target space
    Example
    -------
    """

    register = Workflow(name=name)
         
    inputnode = Node(interface=util.IdentityInterface(fields=['source_files',
                                                                 'target_image',
                                                                 'ROI_File']),
                        name='inputspec')
    outputnode = Node(interface=util.IdentityInterface(fields=['func2std_transform',
                                                                'transformed_ROI'
                                                                  ]),
                         name='outputspec')

    ''' 
    Pipeline:
    #1 : Calculate the mean image from the functional run.
    #2 : Do BET on the mean image.
    #3 : FAST is done for segmenting the White matter.
    #4 : Binarize is done for making the mask out of the file generated after the
         FAST segmentation.
    #5 : Now calculate the mean 2 anamtomical co-registration matrix. 
    #6 : BBR is done for refining the matrix generated from mean2anatomical co-reg.
    #7 : Now, calculate the transformation from anatomical to reference file. 
    #8 : Calculate the mask from the transformed anat2target file.
    #9 : Now , apply transformation using the target_image as reference and 
         mean2anat matrix as transformation matrix. Keep Spline interpolation.
    #10 : Now mask the output with the anat2targetmask so that the values outside 
          the brain go zero. 

    '''

    # meanfunc = MapNode(fsl.ImageMaths(op_string='-Tmean',suffix='_mean'), iterfield = ['in_file'],name = 'meanfunc')
    meanfunc = MapNode(interface=fsl.ExtractROI(t_size=1),
                             iterfield=['in_file', 't_min'],
                             name = 'meanfunc')
    fast = Node(fsl.FAST(), 
          name='fast')
    selectfile = Node(interface=util.Select(index=[2]),
                        name='select')

    
    binarize = Node(fsl.ImageMaths(op_string='-nan -thr 0.5 -bin'),
                   name='binarize')
    mean2ref = MapNode(fsl.FLIRT(interp = 'trilinear'), 
                        iterfield = ['in_file'], 
                        name='mean2ref')
    mean2refbbr = MapNode(fsl.FLIRT(interp = 'trilinear'), 
                    iterfield = ['in_file','in_matrix_file'], 
                    name='mean2refbbr')

    inv_mat = MapNode(fsl.ConvertXFM(invert_xfm = True), 
                      iterfield = ['in_file'], 
                      name = 'inv_mat')
    transform_ROI = MapNode(fsl.ApplyXFM(interp='nearestneighbour'), 
                      iterfield = ['in_matrix_file','reference'],
                      name='transform_ROI')

    register.connect(inputnode, 'source_files', meanfunc, 'in_file')
    register.connect(inputnode, ('source_files', pickmiddle), meanfunc, 't_min')
    """
    Estimate the tissue classes from the anatomical image. But use spm's segment
    as FSL appears to be breaking.
    """
    register.connect(inputnode, 'target_image', fast, 'in_files')
    """
    Binarize the segmentation
    """
    pickindex = lambda x, i: x[i]
    register.connect(fast, 'partial_volume_files', selectfile, 'inlist')

    register.connect(selectfile, 'out', binarize, 'in_file')
    # register.connect(fast, ('partial_volume_files', pickindex, 2),
                     # binarize, 'in_file')

    """
    Calculate rigid transform from mean image to anatomical image
    """
    mean2ref.inputs.dof = 12
    mean2ref.inputs.searchr_x = [-180, 180]
    mean2ref.inputs.searchr_y = [-180, 180]
    mean2ref.inputs.searchr_z = [-180, 180]
    register.connect(meanfunc,'roi_file', mean2ref, 'in_file')
    register.connect(inputnode, 'target_image', mean2ref, 'reference')

    """
    Now use bbr cost function to improve the transform
    """

    mean2refbbr.inputs.dof = 12
    mean2refbbr.inputs.searchr_x = [-180, 180]
    mean2refbbr.inputs.searchr_y = [-180, 180]
    mean2refbbr.inputs.searchr_z = [-180, 180]
    mean2refbbr.inputs.cost = 'bbr'
    mean2refbbr.inputs.schedule = os.path.join(os.getenv('FSLDIR'),
                                                'etc/flirtsch/bbr.sch')
    register.connect(meanfunc,'roi_file', mean2refbbr, 'in_file')
    register.connect(binarize, 'out_file', mean2refbbr, 'wm_seg')
    # register.connect(inputnode, 'anatomical_images', mean2anatbbr, 'reference')
    register.connect(inputnode, 'target_image', mean2refbbr, 'reference')
    register.connect(mean2ref, 'out_matrix_file',
                     mean2refbbr, 'in_matrix_file')

    transform_ROI.inputs.padding_size = 0

    register.connect(mean2refbbr,'out_matrix_file', inv_mat, 'in_file')
    register.connect(inputnode, 'ROI_File', transform_ROI, 'in_file')
    register.connect(inv_mat, 'out_file', transform_ROI, 'in_matrix_file')
    register.connect(meanfunc, 'roi_file', transform_ROI, 'reference')

    datasink_transformedROI = Node(interface=DataSink(), name="datasink_transformedROI")
    datasink_transformedROI.inputs.remove_dest_dir = True
    datasink_func2std = Node(interface=DataSink(), name="datasink_func2std")
    datasink_func2std.inputs.remove_dest_dir = True


    """
    Assign all the output files
    """
    register.connect(mean2refbbr, 'out_matrix_file', 
                      outputnode, 'func2std_transform')
    register.connect(transform_ROI, 'out_file', outputnode, 'transformed_ROI')
    register.connect(outputnode, 'transformed_ROI', datasink_transformedROI, 'out_file')
    register.connect(outputnode, 'func2std_transform', datasink_func2std,'out_file')

    return register



def __reg_workflow(no_subjects, name = 'registration'):

    """
    This function has been deprecated. Please no longer use it.
    Create a FEAT preprocessing workflow
    Parameters
    ----------
    ::
        name : name of workflow (default: 'registration')
    Inputs::
        inputspec.source_files : files (filename or list of filenames to register)
        inputspec.anatomical_images : anatomical images in the same subject-wise order for coregistration.
        inputspec.target_image : registration target
    Outputs::
        outputspec.func2anat_transform : FLIRT transform
        outputspec.transformed_files : transformed files in target space
    Example
    -------
    """

    register = Workflow(name=name)

    inputnode = Node(interface=util.IdentityInterface(fields=['source_files',
                                                                 'anatomical_images',
                                                                 'target_image']),
                        name='inputspec')
    outputnode = Node(interface=util.IdentityInterface(fields=['func2anat_transform',
                                                                'func2std_transform',
                                                                  'transformed_files'
                                                                  ]),
                         name='outputspec')

    ''' 
    Pipeline:
    #1 : Calculate the mean image from the functional run.
    #2 : Do BET on the mean image.
    #3 : FAST is done for segmenting the White matter.
    #4 : Binarize is done for making the mask out of the file generated after the
         FAST segmentation.
    #5 : Now calculate the mean 2 anamtomical co-registration matrix. 
    #6 : BBR is done for refining the matrix generated from mean2anatomical co-reg.
    #7 : Now, calculate the transformation from anatomical to reference file. 
    #8 : Calculate the mask from the transformed anat2target file.
    #9 : Now , apply transformation using the target_image as reference and 
         mean2anat matrix as transformation matrix. Keep Spline interpolation.
    #10 : Now mask the output with the anat2targetmask so that the values outside 
          the brain go zero. 

    '''

    # meanfunc = MapNode(fsl.ImageMaths(op_string='-Tmean',suffix='_mean'), iterfield = ['in_file'],name = 'meanfunc')
    meanfunc = MapNode(interface=fsl.ExtractROI(t_size=1),
                             iterfield=['in_file', 't_min'],
                             name = 'meanfunc')
    stripper = MapNode(fsl.BET(frac = 0.3), 
      iterfield = ['in_file'], name='stripper')
    fast = MapNode(fsl.FAST(), 
          iterfield =['in_files'], 
          name='fast')
    selectfile = MapNode(interface=util.Select(index=[2]), iterfield = ['inlist'],name='select')

    
    binarize = MapNode(fsl.ImageMaths(op_string='-nan -thr 0.5 -bin'),
                   iterfield = ['in_file'],
                   name='binarize')
    mean2anat = MapNode(fsl.FLIRT(interp = 'trilinear'), iterfield = ['in_file','reference'], name='mean2anat')
    mean2anatbbr = MapNode(fsl.FLIRT(interp = 'trilinear'), iterfield = ['in_file','wm_seg','reference','in_matrix_file'], name='mean2anatbbr')
    anat2target_affine = MapNode(fsl.FLIRT(), iterfield = ['in_file'], name='anat2target_linear')
    concat_mat = MapNode(fsl.ConvertXFM(concat_xfm = True), iterfield=['in_file','in_file2'], name='concat_mat')
    # anat2targetmask = MapNode(interface=fsl.BET(mask = True,
    #                                      no_output=True,
    #                                      frac = 0.3),
    #                       iterfield=['in_file'],
    #                       name = 'anat2targetmask')
    anat2targetmask = Node(interface=fsl.BET(mask = True,
                                     no_output=True,
                                     frac = 0.3),
                      name = 'anat2targetmask') 
    warpfile = MapNode(fsl.ApplyXFM(interp='spline'), 
                      iterfield = ['in_file', 'in_matrix_file'],
                      name='txm_registered')
    """
    Mask the functional runs with the extracted mask
    """
    maskWarpFile = MapNode(interface=fsl.ImageMaths(suffix='_masked',
                                               op_string='-mas'),
                      iterfield=['in_file'],
                      name = 'maskWarpFile')

    register.connect(inputnode, 'source_files', meanfunc, 'in_file')
    register.connect(inputnode, ('source_files', pickmiddle), meanfunc, 't_min')
    """
    Estimate the tissue classes from the anatomical image. But use spm's segment
    as FSL appears to be breaking.
    """
    
    register.connect(inputnode, 'anatomical_images', stripper, 'in_file')

    register.connect(stripper, 'out_file', fast, 'in_files')
    """
    Binarize the segmentation
    """
    pickindex = lambda x, i: x[i]
    register.connect(fast, 'partial_volume_files', selectfile, 'inlist')

    register.connect(selectfile, 'out', binarize, 'in_file')
    # register.connect(fast, ('partial_volume_files', pickindex, 2),
                     # binarize, 'in_file')

    """
    Calculate rigid transform from mean image to anatomical image
    """
    mean2anat.inputs.dof = 12
    mean2anat.inputs.searchr_x = [-180, 180]
    mean2anat.inputs.searchr_y = [-180, 180]
    mean2anat.inputs.searchr_z = [-180, 180]
    register.connect(meanfunc,'roi_file', mean2anat, 'in_file')
    register.connect(stripper, 'out_file', mean2anat, 'reference')

    """
    Now use bbr cost function to improve the transform
    """

    mean2anatbbr.inputs.dof = 12
    mean2anatbbr.inputs.searchr_x = [-180, 180]
    mean2anatbbr.inputs.searchr_y = [-180, 180]
    mean2anatbbr.inputs.searchr_z = [-180, 180]
    mean2anatbbr.inputs.cost = 'bbr'
    mean2anatbbr.inputs.schedule = os.path.join(os.getenv('FSLDIR'),
                                                'etc/flirtsch/bbr.sch')
    register.connect(meanfunc,'roi_file', mean2anatbbr, 'in_file')
    register.connect(binarize, 'out_file', mean2anatbbr, 'wm_seg')
    # register.connect(inputnode, 'anatomical_images', mean2anatbbr, 'reference')
    register.connect(stripper, 'out_file', mean2anatbbr, 'reference')
    register.connect(mean2anat, 'out_matrix_file',
                     mean2anatbbr, 'in_matrix_file')
    """
    Calculate affine transform from anatomical to target
    """

    anat2target_affine.inputs.searchr_x = [-180, 180]
    anat2target_affine.inputs.searchr_y = [-180, 180]
    anat2target_affine.inputs.searchr_z = [-180, 180]
    anat2target_affine.inputs.cost = 'corratio'
    register.connect(stripper, 'out_file', anat2target_affine, 'in_file')
    register.connect(inputnode, 'target_image',
                     anat2target_affine, 'reference')
    # register.connect(anat2target_affine,'out_file', anat2targetmask, 'in_file')
    register.connect(inputnode,'target_image', anat2targetmask, 'in_file')

    register.connect(mean2anatbbr, 'out_matrix_file', concat_mat, 'in_file')
    register.connect(anat2target_affine, 'out_matrix_file', concat_mat,'in_file2')
    warpfile.inputs.padding_size = 0

    register.connect(inputnode, 'source_files', warpfile, 'in_file')
    register.connect(concat_mat, 'out_file', warpfile, 'in_matrix_file')
    register.connect(inputnode, 'target_image', warpfile, 'reference')

    register.connect(warpfile,'out_file', maskWarpFile, 'in_file')
    register.connect(anat2targetmask, 'mask_file', maskWarpFile, 'in_file2')


    datasink = Node(interface=DataSink(), name="datasink")
    """
    Assign all the output files
    """
    register.connect(maskWarpFile, 'out_file', outputnode, 'transformed_files')
    register.connect(mean2anatbbr, 'out_matrix_file',
                     outputnode, 'func2anat_transform')
    register.connect(concat_mat, 'out_file', 
                      outputnode, 'func2std_transform')
    register.connect(maskWarpFile, 'out_file', datasink,'out_file')

    return register


