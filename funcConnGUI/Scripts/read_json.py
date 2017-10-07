import numpy as np
import sys
import json
import Preprocess_network as parallelPreproc
import timeit
import os
import shutil
start = timeit.default_timer()
# JSONFile = sys.argv[1]
JSONFile = '/home/deepak/Desktop/FConnectivityAnalysis/FConnectivityAnalysisDesign.json'
with open(JSONFile) as JSON:
    try :
        
        data = json.load(JSON)
        AnalysisName = data['Analysis Name']
    #     print(AnalysisName)
        WorkingDir = data['WorkingDir']
        AnalysisParams = data['AnalysisParams']
        
    except :
        print('Unexpected error:', sys.exc_info()[0])
        raise

def get_TR(in_file):
    import nibabel
    f = nibabel.load(in_file)
    return f.get_header()['dim'][0]#['pixdim'][1:4].tolist(), f.get_header()['pixdim'][4]


def run_Preprocessing(AnalysisParams,FunctionalFiles,StructuralFiles,Group = 0):

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
    TR = get_TR(FunctionalFiles[0])
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
        preproc.run('MultiProc', plugin_args={'n_procs': 4})
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
        preproc.run('MultiProc', plugin_args={'n_procs': 4})

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
        # for i in range(no_subjects):
            # RESULTS_REG_DATASINK = OUTPUT_DIR + '/tmp/%s_%s/datasink/'%(RegistrationName,i)
        #     Reg_WorkFlow = parallelPreproc.reg_workflow(name = '%s_%s'%(RegistrationName ,i))

        #     Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts[i]
        #     Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles[i]
        #     Reg_WorkFlow.inputs.inputspec.target_image = ReferenceFile
        #     Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
        #     Reg_WorkFlow.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        #     regoutputs = Reg_WorkFlow.run()
        #     datasink_results=[]
        #     datasink_results += [each for each in os.listdir(RESULTS_REG_DATASINK) if each.endswith('.json')]
        #     with open(RESULTS_REG_DATASINK + datasink_results[0]) as JSON:
        #         datafile = json.load(JSON)
        #         datafile = datafile[0][1][0][1][0][0]
        #     ProcessedFiles_Address = '%sProcessedFile_sub%s.nii.gz'%(dst,i)
        #     ProcessedFilesDIRADDRESSES += [ProcessedFiles_Address]
        #     shutil.copy(datafile,ProcessedFiles_Address)
        RESULTS_REG_DATASINK = OUTPUT_DIR + '/tmp/%s/datasink/'%RegistrationName
        Reg_WorkFlow = parallelPreproc.reg_workflow(no_subjects,name = RegistrationName)
        if (no_subjects ==1):
            Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts[0]
            Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles[0]
        else:        
            Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts
            Reg_WorkFlow.inputs.inputspec.anatomical_images = StructuralFiles
        Reg_WorkFlow.inputs.inputspec.target_image = ReferenceFile
        Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
        Reg_WorkFlow.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}
        regoutputs = Reg_WorkFlow.run('MultiProc', plugin_args={'n_procs': 4})
        datasink_results=[]
        datasink_results += [each for each in os.listdir(RESULTS_REG_DATASINK) if each.endswith('.json')]
        with open(RESULTS_REG_DATASINK + datasink_results[0]) as JSON:
            datafile = json.load(JSON)

            datafile = datafile[0][1][0][1]
        # # ProcessedFiles_Address = '%sProcessedFile_sub%s.nii.gz'%(dst,i)
        # # ProcessedFilesDIRADDRESSES += [ProcessedFiles_Address]
        # # shutil.copy(datafile,ProcessedFiles_Address)
        datasinkouts=[]
        if (no_subjects==1):
            datasinkouts +=[datafile[0]]
        else:
            datasinkouts += [datafile[i][0] for i in range(no_subjects)]
        print(datafile)

    for j in range(no_subjects):
        ProcessedFiles_Address = '%sProcessedFile_sub%s.nii.gz'%(dst,j)
        ProcessedFilesDIRADDRESSES += [ProcessedFiles_Address]
        shutil.copy(datasinkouts[j],ProcessedFiles_Address)
    return ProcessedFilesDIRADDRESSES

ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
ProcessingWay = AnalysisParams['ProcessingType']['ProcessingWay']
Ngroups = AnalysisParams['No of Groups']
if (ProcessingWay==0):
    FunctionaltxtFiles = AnalysisParams['FilesInfo']['FunctionalFilePaths']
    StructuraltxtFiles = AnalysisParams['FilesInfo']['StructuralFilePaths']

    for i in range(0,Ngroups):
        with open(FunctionaltxtFiles[i]) as file:
            FunctionalFiles_in_this_group = [line.strip('\n') for line in file]
        with open(StructuraltxtFiles[i]) as file:
            StructuralFiles_in_this_group = [line.strip('\n') for line in file]
        Preprocessed_Files = run_Preprocessing(AnalysisParams, FunctionalFiles_in_this_group, StructuralFiles_in_this_group,Group= i)
        print(Preprocessed_Files)

stop = timeit.default_timer()
print("Total time taken for running the program: ", stop - start)
