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


def run_Preprocessing(AnalysisParams,FunctionalFiles,StructuralFiles = None,Group = 0):

    ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
    B0unwarping = AnalysisParams['B0 Unwarping']
    BETextract = AnalysisParams['BET Brain Extract']
    FWHM = AnalysisParams['FWHM']
    TemporalFilt = AnalysisParams['Temporal Filtering']
    HPsigma = AnalysisParams['High Pass Value (in sigma)']
    LPsigma = AnalysisParams['Low Pass Value (in sigma)']
    MotionCorrection = AnalysisParams['Motion Correction']
    Registration = AnalysisParams['Registration']
    SliceTimeCorrect = AnalysisParams['Slice Time Correct']
    Intensity_Norm = AnalysisParams['Intensity Normalization']
    MelodicICA = AnalysisParams['Melodic ICA']
    PerfusionSubtract = AnalysisParams['Perfusion Subtraction']
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

def call_corr_wf(Files_for_corr_dict, atlas_file, mask_file, TEMP_DIR_FOR_STORAGE):
    Corr_calculated_Files = {}
    for group, files in Files_for_corr_dict.items():
        datasink_dest = TEMP_DIR_FOR_STORAGE + '/' + group + '/datasink/'

        corr_wf = corr.build_correlation_wf(name = group)
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





def call_stat_Analysis(Files_for_stats_dict, 
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
                    rejected, FDRCorrected = ttest.fdr_correction(Pvals, is_npy = True)
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
    # JSONFile = sys.argv[1]
    JSONFile = '/home1/ee3140506/FConnectivityAnalysis/FConnectivityAnalysisDesign.json'
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
    TEMP_DIR_FOR_STORAGE = opj(OUTPUT_DIR, 'tmp')


    CorrROImapFiles = {}
    if not (os.path.exists(TEMP_DIR_FOR_STORAGE)):
        os.mkdir(TEMP_DIR_FOR_STORAGE)

    mask_not_provided = False
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
        mask_file = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
    
    if (ProcessingWay==0):
        FunctionaltxtFiles = AnalysisParams['FilesInfo']['FunctionalFilePaths']
        StructuraltxtFiles = AnalysisParams['FilesInfo']['StructuralFilePaths']

        for i in range(Ngroups):
            ProcName = 'CorrCalc_group%s'%i
            CorrROImapFiles[ProcName] = []
            with open(FunctionaltxtFiles[i]) as file:
                FunctionalFiles_in_this_group = [line.strip('\n') for line in file]
            Registration = AnalysisParams['Registration']
            if Registration:
                with open(StructuraltxtFiles[i]) as file:
                    StructuralFiles_in_this_group = [line.strip('\n') for line in file]
            else: StructuralFiles_in_this_group = None
            Preprocessed_Files = run_Preprocessing(AnalysisParams, 
                                                    FunctionalFiles_in_this_group, 
                                                    StructuralFiles_in_this_group = StructuralFiles_in_this_group,
                                                    Group= i)
            print(Preprocessed_Files)
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
    Totaltime=0
    file = open(opj(OUTPUT_DIR,'timesREADME.txt'),'w')
    stop1 = timeit.default_timer()
    file.write("Total time take for Preprocessing: %ss \n"%(stop1 - start))
    Totaltime += stop1 - start
    if (AnalysisParams['doStats']):

        stop1 = timeit.default_timer()
        Corr_calculated_Files = call_corr_wf(CorrROImapFiles, 
                                                ROIFile,
                                                mask_file, 
                                                TEMP_DIR_FOR_STORAGE)
        stop2 = timeit.default_timer()
        file.write("Total time taken for Correlation calculation: %ss \n"%(stop2 - stop1))
        Totaltime += stop2 - stop1

        stop2 = timeit.default_timer()
        # Corr_calculated_Files = {'CorrCalc_group0':[
        # '/home/deepak/Desktop/FConnectivityAnalysis/tmp/CorrCalc_group0/'
        # 'coff_matrix/mapflow/_coff_matrix0/ProcessedFile_sub0_fc_map.npy',
        # '/home/deepak/Desktop/FConnectivityAnalysis/tmp/CorrCalc_group0/'
        # 'coff_matrix/mapflow/_coff_matrix1/ProcessedFile_sub1_fc_map.npy', 
        # '/home/deepak/Desktop/FConnectivityAnalysis/tmp/CorrCalc_group0/'
        # 'coff_matrix/mapflow/_coff_matrix2/ProcessedFile_sub2_fc_map.npy'],
        #                         'CorrCalc_group1':['/home/deepak/Desktop/'
        # 'FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/'
        # '_coff_matrix0/ProcessedFile_sub0_fc_map.npy', 
        # '/home/deepak/Desktop/FConnectivityAnalysis/tmp/CorrCalc_group1/'
        # 'coff_matrix/mapflow/_coff_matrix1/ProcessedFile_sub1_fc_map.npy']}
        #Corr_calculated_Files = {'CorrCalc_group0':['/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix0/ProcessedFile_sub0_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix1/ProcessedFile_sub1_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix2/ProcessedFile_sub2_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix3/ProcessedFile_sub3_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix4/ProcessedFile_sub4_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix5/ProcessedFile_sub5_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix6/ProcessedFile_sub6_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix7/ProcessedFile_sub7_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix8/ProcessedFile_sub8_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix9/ProcessedFile_sub9_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix10/ProcessedFile_sub10_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix11/ProcessedFile_sub11_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix12/ProcessedFile_sub12_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix13/ProcessedFile_sub13_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix14/ProcessedFile_sub14_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix15/ProcessedFile_sub15_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix16/ProcessedFile_sub16_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix17/ProcessedFile_sub17_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix18/ProcessedFile_sub18_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix19/ProcessedFile_sub19_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix20/ProcessedFile_sub20_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix21/ProcessedFile_sub21_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix22/ProcessedFile_sub22_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix23/ProcessedFile_sub23_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix24/ProcessedFile_sub24_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix25/ProcessedFile_sub25_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix26/ProcessedFile_sub26_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix27/ProcessedFile_sub27_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix28/ProcessedFile_sub28_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix29/ProcessedFile_sub29_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix30/ProcessedFile_sub30_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix31/ProcessedFile_sub31_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix32/ProcessedFile_sub32_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix33/ProcessedFile_sub33_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix34/ProcessedFile_sub34_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix35/ProcessedFile_sub35_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix36/ProcessedFile_sub36_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix37/ProcessedFile_sub37_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix38/ProcessedFile_sub38_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix39/ProcessedFile_sub39_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix40/ProcessedFile_sub40_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix41/ProcessedFile_sub41_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix42/ProcessedFile_sub42_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix43/ProcessedFile_sub43_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix44/ProcessedFile_sub44_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix45/ProcessedFile_sub45_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix46/ProcessedFile_sub46_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix47/ProcessedFile_sub47_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix48/ProcessedFile_sub48_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix49/ProcessedFile_sub49_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix50/ProcessedFile_sub50_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix51/ProcessedFile_sub51_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix52/ProcessedFile_sub52_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix53/ProcessedFile_sub53_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix54/ProcessedFile_sub54_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix55/ProcessedFile_sub55_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix56/ProcessedFile_sub56_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix57/ProcessedFile_sub57_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix58/ProcessedFile_sub58_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix59/ProcessedFile_sub59_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix60/ProcessedFile_sub60_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix61/ProcessedFile_sub61_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix62/ProcessedFile_sub62_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix63/ProcessedFile_sub63_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix64/ProcessedFile_sub64_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix65/ProcessedFile_sub65_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix66/ProcessedFile_sub66_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix67/ProcessedFile_sub67_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix68/ProcessedFile_sub68_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix69/ProcessedFile_sub69_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix70/ProcessedFile_sub70_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix71/ProcessedFile_sub71_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix72/ProcessedFile_sub72_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix73/ProcessedFile_sub73_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix74/ProcessedFile_sub74_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix75/ProcessedFile_sub75_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix76/ProcessedFile_sub76_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix77/ProcessedFile_sub77_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group0/coff_matrix/mapflow/_coff_matrix78/ProcessedFile_sub78_fc_map.npy'], 'CorrCalc_group1':['/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix0/ProcessedFile_sub0_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix1/ProcessedFile_sub1_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix2/ProcessedFile_sub2_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix3/ProcessedFile_sub3_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix4/ProcessedFile_sub4_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix5/ProcessedFile_sub5_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix6/ProcessedFile_sub6_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix7/ProcessedFile_sub7_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix8/ProcessedFile_sub8_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix9/ProcessedFile_sub9_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix10/ProcessedFile_sub10_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix11/ProcessedFile_sub11_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix12/ProcessedFile_sub12_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix13/ProcessedFile_sub13_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix14/ProcessedFile_sub14_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix15/ProcessedFile_sub15_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix16/ProcessedFile_sub16_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix17/ProcessedFile_sub17_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix18/ProcessedFile_sub18_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix19/ProcessedFile_sub19_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix20/ProcessedFile_sub20_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix21/ProcessedFile_sub21_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix22/ProcessedFile_sub22_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix23/ProcessedFile_sub23_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix24/ProcessedFile_sub24_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix25/ProcessedFile_sub25_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix26/ProcessedFile_sub26_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix27/ProcessedFile_sub27_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix28/ProcessedFile_sub28_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix29/ProcessedFile_sub29_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix30/ProcessedFile_sub30_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix31/ProcessedFile_sub31_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix32/ProcessedFile_sub32_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix33/ProcessedFile_sub33_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix34/ProcessedFile_sub34_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix35/ProcessedFile_sub35_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix36/ProcessedFile_sub36_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix37/ProcessedFile_sub37_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix38/ProcessedFile_sub38_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix39/ProcessedFile_sub39_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix40/ProcessedFile_sub40_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix41/ProcessedFile_sub41_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix42/ProcessedFile_sub42_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix43/ProcessedFile_sub43_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix44/ProcessedFile_sub44_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix45/ProcessedFile_sub45_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix46/ProcessedFile_sub46_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix47/ProcessedFile_sub47_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix48/ProcessedFile_sub48_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix49/ProcessedFile_sub49_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix50/ProcessedFile_sub50_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix51/ProcessedFile_sub51_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix52/ProcessedFile_sub52_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix53/ProcessedFile_sub53_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix54/ProcessedFile_sub54_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix55/ProcessedFile_sub55_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix56/ProcessedFile_sub56_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix57/ProcessedFile_sub57_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix58/ProcessedFile_sub58_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix59/ProcessedFile_sub59_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix60/ProcessedFile_sub60_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix61/ProcessedFile_sub61_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix62/ProcessedFile_sub62_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix63/ProcessedFile_sub63_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix64/ProcessedFile_sub64_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix65/ProcessedFile_sub65_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix66/ProcessedFile_sub66_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix67/ProcessedFile_sub67_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix68/ProcessedFile_sub68_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix69/ProcessedFile_sub69_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix70/ProcessedFile_sub70_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix71/ProcessedFile_sub71_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix72/ProcessedFile_sub72_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix73/ProcessedFile_sub73_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix74/ProcessedFile_sub74_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix75/ProcessedFile_sub75_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix76/ProcessedFile_sub76_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix77/ProcessedFile_sub77_fc_map.npy', '/home1/ee3140506/FConnectivityAnalysis/tmp/CorrCalc_group1/coff_matrix/mapflow/_coff_matrix78/ProcessedFile_sub78_fc_map.npy']}
	        

        doAnalysiswtGrps = AnalysisParams['Stats']['Analysis within Groups']
        doAnalysisbwGrps = AnalysisParams['Stats']['Analysis between Groups']
        doNormalFisher = AnalysisParams['Stats']['doNormalFisher']
        doSeparateFDR = AnalysisParams['Stats']['Separate FDR']
        # if doAnalysiswtGrps:
            # Do Something. The Function is yet to be defined.

        if ((Ngroups==2)and doAnalysisbwGrps):
            GR1grGr2 = AnalysisParams['Stats']['Gr1>Gr2']
            call_stat_Analysis(Corr_calculated_Files,[1],OUTPUT_DIR, mask_file)
        elif ((Ngroups > 2) and doAnalysisbwGrps):
            combinations = AnalysisParams['Stats']['Combinations']
            call_stat_Analysis(Corr_calculated_Files, 
                                combinations, 
                                OUTPUT_DIR, 
                                mask_file)


        stop = timeit.default_timer()
        file.write("Total time taken for calculating statistics: %ss \n" %(stop - stop2))
        Totaltime += stop - stop2

    print("Total time taken for running the program: ", Totaltime)
    file.write("Total time taken for running the program: %ss \n" %(Totaltime))
    file.close()

