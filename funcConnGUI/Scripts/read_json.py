import numpy as np
import sys
import json
import Preprocess_network as parallelPreproc
import timeit
import os
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
def run_Preprocessing(AnalysisParams,FunctionalFiles,StructuralFiles,Group = 0):

    ReferenceFile = AnalysisParams['ReferSummary']['ReferImgPath']
    B0unwarping = AnalysisParams['B0 Unwarping']
    BETextract = AnalysisParams['BET Brain Extract']
    FWHM = AnalysisParams['FWHM']
    HighPass = AnalysisParams['High Pass']
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
    
    if HighPass:
        preproc = parallelPreproc.create_parallelfeat_preproc(name = FeatProcessName,
                                    highpass= HighPass, 
                                    Intensity_Norm = Intensity_Norm,
                                    BETextract = BETextract,
                                    MotionCorrection = MotionCorrection, 
                                    SliceTimeCorrect = SliceTimeCorrect)
        preproc.inputs.inputspec.highpass = 128./(2*2.5)
        preproc.inputs.inputspec.func = FunctionalFiles
        preproc.inputs.inputspec.fwhm = FWHM
        preproc.base_dir = TEMP_DIR_FOR_STORAGE
        preproc.write_graph(graph2use='colored', format='png', simple_form=True)
        preproc.run('MultiProc', plugin_args={'n_procs': 4})
    else:
        preproc = parallelPreproc.create_parallelfeat_preproc(name = RegistrationName,
                                    highpass= HighPass, 
                                    Intensity_Norm = Intensity_Norm,
                                    BETextract = BETextract,
                                    MotionCorrection = MotionCorrection, 
                                    SliceTimeCorrect = SliceTimeCorrect)
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
    datasinkouts = []
    datasinkouts += [datafile[i][0] for i in range(len(datafile))]

    if Registration:
        datasinkouts_afterreg = []
        for i in range(len(datasinkouts)):
            RESULTS_REG_DATASINK = OUTPUT_DIR + '/tmp/%s_%s/datasink/'%(RegistrationName,i)
            Reg_WorkFlow = parallelPreproc.reg_workflow(name = '%s_%s'%(RegistrationName ,i))

            Reg_WorkFlow.inputs.inputspec.source_files = datasinkouts[i]
            Reg_WorkFlow.inputs.inputspec.anatomical_image = StructuralFiles[i]
            Reg_WorkFlow.inputs.inputspec.target_image = ReferenceFile
            Reg_WorkFlow.base_dir = TEMP_DIR_FOR_STORAGE
            regoutputs = Reg_WorkFlow.run()
            datasink_results=[]
            datasink_results += [each for each in os.listdir(RESULTS_REG_DATASINK) if each.endswith('.json')]
            with open(RESULTS_REG_DATASINK + datasink_results[0]) as JSON:
                datafile = json.load(JSON)
                datafile = datafile[0][1][0][1]
            
            datasinkouts_afterreg += [datafile[i][0] for i in range(len(datafile))]
        return datasinkouts_afterreg
    else:
        return datasinkouts

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
