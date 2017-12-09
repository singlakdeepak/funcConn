'''
Author: Deepak Singla (singlakdeepak5@gmail.com)

The registration currently is being done just after preprocessing. 
It has to be changed and the registration of ROI to each functional
file is required. It can be done by inverse matrix of func2std.
After calculating the correlations for each subject, register them to 
standard and then mask it so that statistics can be calculated properly.

TODO
1. The reg procedure can be kept as it is, just the inverse matrix is 
needed and then apply inverse matrix to ROI for each subject.
2. After the corrs are calculated, now func2std can be done. But here I need
the corr files in the form of brains rather than matrices. Changes will be required 
in corr function because I pick up the brain voxels according to the mask. 

'''
import numpy as np
import sys
import json
import Preprocess_network as parallelPreproc
import timeit
import os
import shutil
import subprocess
from nipype.interfaces import fsl
import corr
import ttest
from os.path import join as opj
threads = 8
def get_TR(in_file):
    import nibabel
    f = nibabel.load(in_file)
    return f.get_header()['dim'][0]


def __run_Preprocessing(AnalysisParams,FunctionalFiles,StructuralFiles = None,Group = 0):
    '''
    This function has been deprecated because it was being used for registering the functional 
    images first. Rather registration shouldn't be done now. 
    '''
    ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
    B0unwarping = AnalysisParams['B0 Unwarping']
    BETextract = AnalysisParams['BET Brain Extract']
    if BETextract:
        BETextractvalue = AnalysisParams['BET Correction Value']
    else:
        BETextractvalue = 0
    FWHM = AnalysisParams['FWHM']
    TemporalFilt = AnalysisParams['Temporal Filtering']
    HPsigma = AnalysisParams['High Pass Value (in sigma)']
    LPsigma = AnalysisParams['Low Pass Value (in sigma)']
    MotionCorrection = AnalysisParams['Motion Correction']
    Registration = AnalysisParams['Registration']
    SliceTimeCorrect = AnalysisParams['Slice Time Correct']
    Intensity_Norm = AnalysisParams['Intensity Normalization']
    MelodicICA = AnalysisParams['Melodic ICA']
    # PerfusionSubtract = AnalysisParams['Perfusion Subtraction']
    OUTPUT_DIR = AnalysisParams['OutputInfo']['OutDirectory']
    TEMP_DIR_FOR_STORAGE = OUTPUT_DIR + '/tmp'
    FeatProcessName = 'featpreproc_group%s'%Group
    RegistrationName = 'registration_group%s'%Group
    RESULTS_FEAT_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink/'%FeatProcessName
    TR = AnalysisParams['Repetition Time']
    print('Using Repetition time: %s'%TR)


    if TemporalFilt:
        preproc = parallelPreproc.create_parallelfeat_preproc(name = FeatProcessName,
                                    highpass= TemporalFilt, 
                                    Intensity_Norm = Intensity_Norm,
                                    BETextract = BETextract,
                                    BETvalue = BETextractvalue,
                                    MotionCorrection = MotionCorrection,                                                
                                    SliceTimeCorrect = SliceTimeCorrect,
                                    time_repeat = TR)

        preproc.inputs.inputspec.highpass = (HPsigma,LPsigma)                                                                                                                                                                                                                                                                                                                                                       
        preproc.inputs.inputspec.func = FunctionalFiles
        preproc.inputs.inputspec.fwhm = FWHM
        preproc.base_dir = TEMP_DIR_FOR_STORAGE
        preproc.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        preproc.write_graph(graph2use='colored', format='png', simple_form=True)
        preproc.run('MultiProc', plugin_args={'n_procs': threads})          
    else:
        preproc = parallelPreproc.create_parallelfeat_preproc(name = FeatProcessName,
                                    highpass= TemporalFilt, 
                                    Intensity_Norm = Intensity_Norm,
                                    BETextract = BETextract,
                                    BETvalue = BETextractvalue,
                                    MotionCorrection = MotionCorrection, 
                                    SliceTimeCorrect = SliceTimeCorrect,
                                    time_repeat = TR)
        preproc.inputs.inputspec.func = FunctionalFiles
        preproc.inputs.inputspec.fwhm = FWHM
        preproc.base_dir = TEMP_DIR_FOR_STORAGE
        preproc.write_graph(graph2use='colored', format='png', simple_form=True)
        preproc.run('MultiProc', plugin_args={'n_procs': threads})

    datasink_results=[]
    datasink_results += [each for each in os.listdir(RESULTS_FEAT_DATASINK) if each.endswith('.json')]
    with open(RESULTS_FEAT_DATASINK + datasink_results[0]) as JSON:
        datafile = json.load(JSON)
        datafile =datafile[0][1][0][1]
    no_subjects = len(datafile)
    datasinkouts = []
    datasinkouts += [datafile[i][0] for i in range(no_subjects)]

    ProcessedFilesDIRADDRESSES = []
    dst = OUTPUT_DIR + '/Preprocessing_group%s/'%Group
    if not (os.path.exists(dst)):
        os.mkdir(dst)
    else:
        shutil.rmtree(dst)
        os.mkdir(dst)

    if Registration:
        datasinkouts_afterreg = []
        RESULTS_REG_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink/'%RegistrationName
        Reg_WorkFlow = parallelPreproc.reg_workflow(no_subjects,name = RegistrationName)
        # if (no_subjects ==1):
        #     Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts[0]
        #     Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles[0]
        # else:        
        Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts
        Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles
        Reg_WorkFlow.inputs.inputspec.target_image = ReferenceFile
        Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
        Reg_WorkFlow.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        regoutputs = Reg_WorkFlow.run('MultiProc', plugin_args={'n_procs': threads})
        datasink_results=[]
        datasink_results += [each for each in os.listdir(RESULTS_REG_DATASINK) if each.endswith('.json')]
        with open(RESULTS_REG_DATASINK + datasink_results[0]) as JSON:
            datafile = json.load(JSON)

            datafile = datafile[0][1][0][1]

        datasinkouts=[]
        # if (no_subjects==1):
        #     datasinkouts +=[datafile[0]]
        # else:
        datasinkouts += [datafile[i][0] for i in range(no_subjects)]

    for j in range(no_subjects):
        ProcessedFiles_Address = '%sProcessedFile_sub%s.nii.gz'%(dst,j)
        ProcessedFilesDIRADDRESSES += [ProcessedFiles_Address]
        shutil.copy(datasinkouts[j],ProcessedFiles_Address)
    return ProcessedFilesDIRADDRESSES

def run_Preprocessing(AnalysisParams,
                        FunctionalFiles = None,
                        StructuralFiles = None,
                        FeatFiles = None,
                        doPreprocessing= True, Group = 0):

    ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
    ROIFile = AnalysisParams['ROIFilePath']
    Registration = AnalysisParams['Registration']
    # PerfusionSubtract = AnalysisParams['Perfusion Subtraction']
    OUTPUT_DIR = AnalysisParams['OutputInfo']['OutDirectory']
    TEMP_DIR_FOR_STORAGE = OUTPUT_DIR + '/tmp'
    RegistrationName = 'registration_group%s'%Group
    TR = AnalysisParams['Repetition Time']
    print('Using Repetition time: %s'%TR)

    if doPreprocessing:
        B0unwarping = AnalysisParams['B0 Unwarping']
        BETextract = AnalysisParams['BET Brain Extract']
        if BETextract:
            BETextractvalue = AnalysisParams['BET Correction Value']
        else:
            BETextractvalue = 0
        FWHM = AnalysisParams['FWHM']
        TemporalFilt = AnalysisParams['Temporal Filtering']
        HPsigma = AnalysisParams['High Pass Value (in sigma)']
        LPsigma = AnalysisParams['Low Pass Value (in sigma)']
        MotionCorrection = AnalysisParams['Motion Correction']
        SliceTimeCorrect = AnalysisParams['Slice Time Correct']
        Intensity_Norm = AnalysisParams['Intensity Normalization']
        MelodicICA = AnalysisParams['Melodic ICA']
        FeatProcessName = 'featpreproc_group%s'%Group
        RESULTS_FEAT_DATASINK = OUTPUT_DIR +'/tmp/%s/datasink/'%FeatProcessName

        if TemporalFilt:
            preproc = parallelPreproc.create_parallelfeat_preproc(name = FeatProcessName,
                                        highpass= TemporalFilt, 
                                        Intensity_Norm = Intensity_Norm,
                                        BETextract = BETextract,
                                        BETvalue = BETextractvalue,
                                        MotionCorrection = MotionCorrection,                                                
                                        SliceTimeCorrect = SliceTimeCorrect,
                                        time_repeat = TR)

            preproc.inputs.inputspec.highpass = (HPsigma,LPsigma)                                                                                                                                                                                                                                                                                                                                                       
            preproc.inputs.inputspec.func = FunctionalFiles
            preproc.inputs.inputspec.fwhm = FWHM
            preproc.base_dir = TEMP_DIR_FOR_STORAGE
            preproc.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
            preproc.write_graph(graph2use='colored', format='png', simple_form=True)
            preproc.run('MultiProc', plugin_args={'n_procs': threads})          
        else:
            preproc = parallelPreproc.create_parallelfeat_preproc(name = FeatProcessName,
                                        highpass= TemporalFilt, 
                                        Intensity_Norm = Intensity_Norm,
                                        BETextract = BETextract,
                                        BETvalue = BETextractvalue,
                                        MotionCorrection = MotionCorrection, 
                                        SliceTimeCorrect = SliceTimeCorrect,
                                        time_repeat = TR)
            preproc.inputs.inputspec.func = FunctionalFiles
            preproc.inputs.inputspec.fwhm = FWHM
            preproc.base_dir = TEMP_DIR_FOR_STORAGE
            preproc.write_graph(graph2use='colored', format='png', simple_form=True)
            preproc.run('MultiProc', plugin_args={'n_procs': threads})

        datasink_results=[]
        datasink_results += [each for each in os.listdir(RESULTS_FEAT_DATASINK) if each.endswith('.json')]
        with open(RESULTS_FEAT_DATASINK + datasink_results[0]) as JSON:
            datafile = json.load(JSON)
            datafile =datafile[0][1][0][1]
        no_subjects = len(datafile)
        datasinkouts = []
        datasinkouts += [datafile[i][0] for i in range(no_subjects)]

        ProcessedFilesDIRADDRESSES = []
        dst = OUTPUT_DIR + '/Preprocessing_group%s/'%Group
        if not (os.path.exists(dst)):
            os.mkdir(dst)
        else:
            shutil.rmtree(dst)
            os.mkdir(dst)
        for j in range(no_subjects):
            ProcessedFiles_Address = '%sProcessedFile_sub%s.nii.gz'%(dst,j)
            ProcessedFilesDIRADDRESSES += [ProcessedFiles_Address]
            shutil.copy(datasinkouts[j],ProcessedFiles_Address)

    if Registration:
        if FeatFiles==None:
            no_subjects = len(FunctionalFiles)
            datasinkouts_afterreg = []
            ROI_REG_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink_transformedROI/'%RegistrationName
            func2std_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink_func2std/'%RegistrationName
            if StructuralFiles!=None:
                Reg_WorkFlow = parallelPreproc.reg_workflow_with_Anat(no_subjects,
                                                            name = RegistrationName)
                Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles
            else:
                Reg_WorkFlow = parallelPreproc.reg_workflow_wo_Anat(no_subjects,
                                                            name = RegistrationName)              
            Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts
            Reg_WorkFlow.inputs.inputspec.target_image = ReferenceFile
            Reg_WorkFlow.inputs.inputspec.ROI_File = ROIFile
            Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
            Reg_WorkFlow.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
            regoutputs = Reg_WorkFlow.run('MultiProc', plugin_args={'n_procs': threads})
            ROIsink_results=[]
            ROIsink_results += [each for each in os.listdir(ROI_REG_DATASINK) if each.endswith('.json')]
            func2std_results = []
            func2std_results += [each for each in os.listdir(func2std_DATASINK) if each.endswith('.json')]

            with open(ROI_REG_DATASINK + ROIsink_results[0]) as JSON:
                transROIfile = json.load(JSON)
                transROIfile = transROIfile[0][1][0][1]

            with open(func2std_DATASINK + func2std_results[0]) as JSON:
                func2stdfile = json.load(JSON)
                func2stdfile = func2stdfile[0][1][0][1]

            ROIsinkouts=[]
            ROIsinkouts += [transROIfile[i][0] for i in range(no_subjects)]
            func2stdsinkouts = []
            func2stdsinkouts += [func2stdfile[i][0] for i in range(no_subjects)]
            print('Registered ROI->func ROIs: ', ROIsinkouts)
            print('func2stdtransforms: ', func2stdsinkouts)
            return ProcessedFilesDIRADDRESSES, ROIsinkouts, func2stdsinkouts
        
        else:

            no_subjects = len(FeatFiles)
            ProcessedFuncouts = [opj(each, 'filtered_func_data.nii.gz') for each in FeatFiles]
            func2std_DATASINK = [opj(each, 'reg/example_func2standard.mat') for each in FeatFiles]
            ROI_REG_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink_transformedROI/'%RegistrationName
            Reg_WorkFlow = parallelPreproc.ROI_transformation(name = RegistrationName)
            Reg_WorkFlow.inputs.inputspec.source_files = ProcessedFuncouts
            Reg_WorkFlow.inputs.inputspec.ROI_File = ROI_File
            Reg_WorkFlow.inputs.inputspec.func2std = func2std_DATASINK
            Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
            Reg_WorkFlow.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
            regoutputs = Reg_WorkFlow.run('MultiProc', plugin_args={'n_procs': threads})
            ROIsink_results=[]
            ROIsink_results += [each for each in os.listdir(ROI_REG_DATASINK) if each.endswith('.json')]            
            with open(ROI_REG_DATASINK + ROIsink_results[0]) as JSON:
                transROIfile = json.load(JSON)
                transROIfile = transROIfile[0][1][0][1]
            ROIsinkouts=[]
            ROIsinkouts += [transROIfile[i][0] for i in range(no_subjects)]
            return ProcessedFuncouts, ROIsinkouts, func2std_DATASINK

    return ProcessedFilesDIRADDRESSES

def call_corr_wf_with_reg(Files_for_corr_dict, atlas_files_dict,
                    func2stdDict,
                     mask_file, reference,
                     TEMP_DIR_FOR_STORAGE):
    Corr_calculated_Files = {}
    for group, files in Files_for_corr_dict.items():
        datasink_dest = TEMP_DIR_FOR_STORAGE + '/' + group + '/datasink/'

        corr_wf = corr.build_correlation_wf(Registration = True,name = group)
        corr_wf.inputs.inputspec.in_files = files

        # print(atlas_files_dict[group])
        # print(func2stdDict[group])        
        
        corr_wf.inputs.inputspec.atlas_files = atlas_files_dict[group]
        corr_wf.inputs.inputspec.func2std = func2stdDict[group]
        corr_wf.inputs.inputspec.reference = reference

        corr_wf.inputs.inputspec.mask_file = mask_file
        corr_wf.base_dir = TEMP_DIR_FOR_STORAGE
        corr_wf.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        corr_wf.run('MultiProc', plugin_args = {'n_procs': threads})
        datasink_results=[]
        datasink_results += [each for each in os.listdir(datasink_dest) if each.endswith('.json')]
        with open(datasink_dest + datasink_results[0]) as JSON:
            datafile = json.load(JSON)

            datafile = datafile[0][1][0][1]

        datasinkouts=[]
        datasinkouts += [datafile[i][0] for i in range(len(files))]
        print(datasinkouts)
        Corr_calculated_Files[group] = datasinkouts
    return Corr_calculated_Files

def call_corr_wf(Files_for_corr_dict, atlas_file, mask_file, TEMP_DIR_FOR_STORAGE):
    Corr_calculated_Files = {}
    for group, files in Files_for_corr_dict.items():
        datasink_dest = TEMP_DIR_FOR_STORAGE + '/' + group + '/datasink/'

        corr_wf = corr.build_correlation_wf(Registration = False, name = group)
        corr_wf.inputs.inputspec.in_files = files
        corr_wf.inputs.inputspec.atlas_file = atlas_file
        corr_wf.inputs.inputspec.mask_file = mask_file
        corr_wf.base_dir = TEMP_DIR_FOR_STORAGE
        corr_wf.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        corr_wf.run('MultiProc', plugin_args = {'n_procs': threads})
        datasink_results=[]
        datasink_results += [each for each in os.listdir(datasink_dest) if each.endswith('.json')]
        with open(datasink_dest + datasink_results[0]) as JSON:
            datafile = json.load(JSON)

            datafile = datafile[0][1][0][1]

        datasinkouts=[]
        datasinkouts += [datafile[i][0] for i in range(len(files))]
        print(datasinkouts)
        Corr_calculated_Files[group] = datasinkouts
    return Corr_calculated_Files


def call_stat_Analysis_wt_grps(Files_for_stats_dict,
                                destination,
                                mask_file,
                                applyFDR = True):
    total_groups = len(Files_for_stats_dict)
    i = 0
    if (total_groups>1): return
    else:
        NumpyFileList = []
        FileNamesForSaving = []
        ProcName = 'CorrCalc_group%s'%i
        MeanGr1, (Tvals, Pvals) = ttest.ttest_1samp_ROIs_if_npy(
                                            Files_for_stats_dict[ProcName],
                                            save_pval_in_log_fmt = False,
                                            applyFisher = True)
        Tvals = ttest.convert_ma_to_np(Tvals)
        np.save(opj(destination,'Tvals_group_%s.npy'%i),
                                             Tvals)
        NumpyFileList.append(Tvals)
        FileNamesForSaving.append(opj(destination,'Tvals_group_%s.nii.gz'%i))
        np.save(opj(destination,'Pvals_Normal_group_%s'%i),
                                     ttest.convert_ma_to_np(Pvals))
       


        Pvalues_in_log = ttest.convert_pvals_to_log_fmt(Pvals,
                                                        Sample_mean_ArrayA = MeanGr1)
        np.save(opj(destination,'Pvals_group_%s.npy'%i),
                                     ttest.convert_ma_to_np(Pvalues_in_log))
        NumpyFileList.append(Pvalues_in_log)
        FileNamesForSaving.append(opj(destination,'Pvals_group_%s.nii.gz'%i))
        if  applyFDR :
            rejected, FDRCorrected = ttest.fdr_correction(Pvals, procs = threads, is_npy = True)
            np.save(opj(destination,'Qvals_Normal_group_%s.npy'%i),
                                    FDRCorrected)

            Qvals_in_log = ttest.convert_pvals_to_log_fmt(FDRCorrected,
                                                            Sample_mean_ArrayA = MeanGr1)
            rejected, Qvals_in_log = ttest.convert_ma_to_np(rejected),\
                                         ttest.convert_ma_to_np(Qvals_in_log)

            np.save(opj(destination,'Qvals_group_%s.npy'%i),
                                    Qvals_in_log)
            NumpyFileList.append(Qvals_in_log)
            FileNamesForSaving.append(opj(destination,'Qvals_group_%s.nii.gz'%i))
        ttest.make_brain_back_from_npy(NumpyFileList,FileNamesForSaving, mask_file)
    print('Analysis Over.')


def call_stat_Analysis_bw_grps(Files_for_stats_dict, 
                        combinations_reqd, 
                        destination, 
                        mask_file, 
                        applyFDR = True):
    tell_combs = len(combinations_reqd)
    total_groups = len(Files_for_stats_dict)
    previous = 0
    next_start = 0
    newCombs_to_take = total_groups-1
    while (previous<total_groups - 1):
        temp_combs = combinations_reqd[next_start:(next_start + newCombs_to_take)]
        for i in range(previous, total_groups-1):
            if (temp_combs[i] == 1):
                NumpyFileList = []
                FileNamesForSaving = []
                ProcName1 = 'CorrCalc_group%s'%previous
                ProcName2 = 'CorrCalc_group%s'%(previous + i + 1)
                MeanGr1 , MeanGr2, (Tvals, Pvals) = ttest.ttest_ind_samples_if_npy(
                                                                Files_for_stats_dict[ProcName1],
                                                                Files_for_stats_dict[ProcName2],
                                                                save_pval_in_log_fmt = False,
                                                                equal_var= False, 
                                                                applyFisher = True)
                Tvals = ttest.convert_ma_to_np(Tvals)
                np.save(opj(destination,'Tvals_group_{}_group_{}.npy'.format(
                                            previous,
                                            previous + i + 1)),
                                             Tvals)
                NumpyFileList.append(Tvals)
                FileNamesForSaving.append(opj(destination,'Tvals_group_{}_group_{}.nii.gz'.format(
                                            previous,
                                            previous + i + 1)))




                np.save(opj(destination,'Pvals_Normal_group_{}_group_{}.npy'.format(
                                            previous,
                                            previous + i + 1)),
                                             ttest.convert_ma_to_np(Pvals))
               


                Pvalues_in_log = ttest.convert_pvals_to_log_fmt(Pvals,
                                            Sample_mean_ArrayA = MeanGr1,
                                            Sample_mean_ArrayB = MeanGr2)
                np.save(opj(destination,'Pvals_group_{}_group_{}.npy'.format(
                                            previous,
                                            previous + i + 1)),
                                             ttest.convert_ma_to_np(Pvalues_in_log))
                NumpyFileList.append(Pvalues_in_log)
                FileNamesForSaving.append(opj(destination,'Pvals_group_{}_group_{}.nii.gz'.format(
                                            previous,
                                            previous + i + 1)))
                if  applyFDR :
                    rejected, FDRCorrected = ttest.fdr_correction(Pvals, procs = threads, is_npy = True)
                    np.save(opj(destination,'Qvals_Normal_group_{}_group_{}.npy'.format(
                                            previous,
                                            previous + i + 1)),
                                            FDRCorrected)

                    Qvals_in_log = ttest.convert_pvals_to_log_fmt(FDRCorrected,
                                            Sample_mean_ArrayA = MeanGr1,
                                            Sample_mean_ArrayB = MeanGr2)
                    rejected, Qvals_in_log = ttest.convert_ma_to_np(rejected),\
                                                 ttest.convert_ma_to_np(Qvals_in_log)

                    np.save(opj(destination,'Qvals_group_{}_group_{}.npy'.format(
                                            previous,
                                            previous + i + 1)),
                                            Qvals_in_log)
                    NumpyFileList.append(Qvals_in_log)
                    FileNamesForSaving.append(opj(destination,'Qvals_group_{}_group_{}.nii.gz'.format(
                                            previous,
                                            previous + i + 1)))
                ttest.make_brain_back_from_npy(NumpyFileList,FileNamesForSaving, mask_file)

        next_start += total_groups - previous -1
        newCombs_to_take = total_groups - previous - 2
        previous += 1
        print('Analysis Over.')


if __name__ == '__main__':

    start = timeit.default_timer()
    JSONFile = sys.argv[1]
    # JSONFile = '/home1/ee3140506/FConnectivityAnalysis/FConnectivityAnalysisDesign.json'
    # JSONFile = '/home/deepak/Desktop/FConnectivityAnalysis/FConnectivityAnalysisDesign.json'
    with open(JSONFile) as JSON:
        try :
            
            data = json.load(JSON)
            AnalysisName = data['Analysis Name']
            threads = data['No of threads']
        #     print(AnalysisName)
            WorkingDir = data['WorkingDir']
            AnalysisParams = data['AnalysisParams']
            
        except :
            print('Unexpected error:', sys.exc_info()[0])
            raise
    ROIFile = AnalysisParams['ROIFilePath']
    ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
    ProcessingWay = AnalysisParams['ProcessingType']['ProcessingWay']
    Ngroups = AnalysisParams['No of Groups']
    OUTPUT_DIR = AnalysisParams['OutputInfo']['OutDirectory']
    Registration = AnalysisParams['Registration']
    TEMP_DIR_FOR_STORAGE = opj(OUTPUT_DIR, 'tmp')


    CorrROImapFiles = {}
    transformedROIsDict = {}
    func2stdDict = {}
    if not (os.path.exists(TEMP_DIR_FOR_STORAGE)):
        os.mkdir(TEMP_DIR_FOR_STORAGE)

    mask_not_provided = AnalysisParams['MaskNotProvided']
    if  mask_not_provided:
        nonzeroportion_mask = fsl.BET(in_file = ReferenceFile,
                                  out_file = TEMP_DIR_FOR_STORAGE + '/mask_for_ttest', 
                                  mask = True, 
                                  no_output=True,
                                  frac = 0.3)
        nonzeroportion_mask.run()
        mask_file = opj(TEMP_DIR_FOR_STORAGE,'mask_for_ttest_mask.nii.gz')
        print('Made mask for the ttest: %s/mask_for_ttest.nii.gz'%TEMP_DIR_FOR_STORAGE)
    else:
        # mask_file = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
        mask_file = AnalysisParams['UseMaskFile']
    
    if (ProcessingWay==0):
        FunctionaltxtFiles = AnalysisParams['FilesInfo']['FunctionalFilePaths']
        StructuraltxtFiles = AnalysisParams['FilesInfo']['StructuralFilePaths']

        for i in range(Ngroups):
            ProcName = 'CorrCalc_group%s'%i
            # CorrROImapFiles[ProcName] = []
            with open(FunctionaltxtFiles[i]) as file:
                FunctionalFiles_in_this_group = [line.strip('\n') for line in file]
            
            if Registration:
                if (len(StructuraltxtFiles)!=0):
                    with open(StructuraltxtFiles[i]) as file:
                        StructuralFiles_in_this_group = [line.strip('\n') for line in file]
                else: StructuralFiles_in_this_group = None
                Preprocessed_Files, transformedROIs, func2stdtransforms = run_Preprocessing(AnalysisParams, 
                                                                            FunctionalFiles = FunctionalFiles_in_this_group, 
                                                                            StructuralFiles = StructuralFiles_in_this_group,
                                                                            Group= i)
                print('Preprocessed Files are: ', Preprocessed_Files)
                CorrROImapFiles[ProcName] = Preprocessed_Files
                transformedROIsDict[ProcName] = transformedROIs
                func2stdDict[ProcName] = func2stdtransforms
            else: 
                StructuralFiles_in_this_group = None
                Preprocessed_Files = run_Preprocessing(AnalysisParams, 
                                                    FunctionalFiles = FunctionalFiles_in_this_group, 
                                                    StructuralFiles = StructuralFiles_in_this_group,
                                                    Group= i)
                print('Preprocessed Files are: ', Preprocessed_Files)
                CorrROImapFiles[ProcName] = Preprocessed_Files

        # for j in range(len(Preprocessed_Files)):
        #     args = ("../../fconn.o", "-i", ProcessedFileName, "-o", dst + 'sub%d'%j, "-r", ROIFile)
        #     popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        #     popen.wait()
        #     output = popen.stdout.read()
        #     print(output)

    elif (ProcessingWay == 1):
        '''
        Already Preprocessed without the use of FSL. We are supplied the list of 
        Functional and Structural Files.
        '''
        FunctionaltxtFiles = AnalysisParams['FilesInfo']['FunctionalFilePaths']
        StructuraltxtFiles = AnalysisParams['FilesInfo']['StructuralFilePaths']
        for i in range(Ngroups):
            ProcName = 'CorrCalc_group%s'%i
            CorrROImapFiles[ProcName] = []
            with open(FunctionaltxtFiles[i]) as file:
                FunctionalFiles_in_this_group = [line.strip('\n') for line in file]
            if (len(StructuralFiles)!=0):
                with open(StructuraltxtFiles[i]) as file:
                    StructuralFiles_in_this_group = [line.strip('\n') for line in file]        
            print('Functional Files in this group: ',FunctionalFiles_in_this_group)
            CorrROImapFiles[ProcName] = FunctionalFiles_in_this_group

    elif (ProcessingWay ==2):
        '''
        Already Preprocessed with the use of FSL. I assume it has a general structure 
        and preprocessed functional files are present in folders. 
        <FEATFolder>/filtered_func_data.nii.gz
        Here I also assume that the preprocessed functional image isn't registered to 
        the standard space. So, every ROI has to be registered to the functional image. 
        '''
        FeatFiles = AnalysisParams['FeatFilesInfo']['FeatFilePaths']
        for i in range(Ngroups):
            ProcName = 'CorrCalc_group%s'%i
            with open(FeatFiles[i]) as file:
                Feat_subs_in_this_group = [line.strip('\n') for line in file]
            print('Feat Folders in this group: ',Feat_subs_in_this_group)
            Preprocessed_Files, transformedROIs, func2stdtransforms = run_Preprocessing(AnalysisParams, 
                                                    FeatFiles = Feat_subs_in_this_group,
                                                    doPreprocessing = False,
                                                    Group= i)
            CorrROImapFiles[ProcName] = Preprocessed_Files
            transformedROIsDict[ProcName] = transformedROIs
            func2stdDict[ProcName] = func2stdtransforms

    Totaltime=0
    file = open(opj(OUTPUT_DIR,'timesREADME.txt'),'w')
    stop1 = timeit.default_timer()
    file.write("Total time take for Preprocessing: %ss \n"%(stop1 - start))
    Totaltime += stop1 - start
    if (AnalysisParams['doStats']):

        stop1 = timeit.default_timer()
        if Registration:
            Corr_calculated_Files = call_corr_wf_with_reg(CorrROImapFiles,
                                                transformedROIsDict,
                                                func2stdDict,
                                                # ROI_File,
                                                mask_file,
                                                ReferenceFile,
                                                TEMP_DIR_FOR_STORAGE)
        else:
            Corr_calculated_Files = call_corr_wf(CorrROImapFiles, 
                                                ROIFile,
                                                mask_file, 
                                                TEMP_DIR_FOR_STORAGE)
        stop2 = timeit.default_timer()
        file.write("Total time taken for Correlation calculation: %ss \n"%(stop2 - stop1))
        Totaltime += stop2 - stop1

        stop2 = timeit.default_timer()

        # Corr_calculated_Files = {'CorrCalc_group0':['/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix0/ProcessedFile_sub0_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix1/ProcessedFile_sub1_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix2/ProcessedFile_sub2_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix3/ProcessedFile_sub3_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix4/ProcessedFile_sub4_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix5/ProcessedFile_sub5_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix6/ProcessedFile_sub6_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix7/ProcessedFile_sub7_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix8/ProcessedFile_sub8_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix9/ProcessedFile_sub9_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix10/ProcessedFile_sub10_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix11/ProcessedFile_sub11_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix12/ProcessedFile_sub12_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix13/ProcessedFile_sub13_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix14/ProcessedFile_sub14_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix15/ProcessedFile_sub15_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix16/ProcessedFile_sub16_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix17/ProcessedFile_sub17_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix18/ProcessedFile_sub18_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix19/ProcessedFile_sub19_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix20/ProcessedFile_sub20_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix21/ProcessedFile_sub21_fc_map.npy', '/home1/ee3140506/NYU_Cocaine_Analysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix22/ProcessedFile_sub22_fc_map.npy']}

        doAnalysiswtGrps = AnalysisParams['Stats']['Analysis within Groups']
        doAnalysisbwGrps = AnalysisParams['Stats']['Analysis between Groups']
        doNormalFisher = AnalysisParams['Stats']['doNormalFisher']
        doSeparateFDR = AnalysisParams['Stats']['Separate FDR']
        # if doAnalysiswtGrps:
            # Do Something. The Function is yet to be defined.

        # You can club the bw groups and wt groups correlations together 
        # because mean and std are already calculated in that case.
        if ((Ngroups==2)and doAnalysisbwGrps):
            GR1grGr2 = AnalysisParams['Stats']['Gr1>Gr2']
            call_stat_Analysis_bw_grps(Corr_calculated_Files,[1],OUTPUT_DIR, mask_file)
        elif ((Ngroups > 2) and doAnalysisbwGrps):
            combinations = AnalysisParams['Stats']['Combinations']
            call_stat_Analysis_bw_grps(Corr_calculated_Files, 
                                combinations, 
                                OUTPUT_DIR, 
                                mask_file)
        if doAnalysiswtGrps:
            call_stat_Analysis_wt_grps(Corr_calculated_Files,OUTPUT_DIR,mask_file)

        stop = timeit.default_timer()
        file.write("Total time taken for calculating statistics: %ss \n" %(stop - stop2))
        Totaltime += stop - stop2

    print("Total time taken for running the program: ", Totaltime)
    file.write("Total time taken for running the program: %ss \n" %(Totaltime))
    file.close()

