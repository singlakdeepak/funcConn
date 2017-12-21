from os.path import join as opj
import nibabel as nib
import numpy as np
import nipype.interfaces.utility as util
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode

from nipype.interfaces import fsl
import os

def pearsonr_with_roi_mean_w_reg(in_file, atlas_file):
    import nibabel as nib
    import numpy as np
    from os.path import join as opj
    import os
    std_param = 1e-2
    sub_id = in_file.split('/')[-1].split('.')[0]
    fc_file_name = sub_id + '_fc_map.nii.gz'
    coff_matrix_file = opj(os.getcwd(),fc_file_name)
    if os.path.exists(coff_matrix_file):
        print('Saved file in : %s'%coff_matrix_file)
        
    else:

        atlasObject = nib.load(atlas_file)
        atlas = atlasObject.get_data()
        num_ROIs = int(np.max(atlas) - np.min(atlas)) 

        print('Min Index:', np.min(atlas),'Max Index', np.max(atlas))
        print('Total Number of Parcellations = ',num_ROIs)
        brain_data = nib.load(in_file)
        brain = brain_data.get_data()
        brain_affine = brain_data.affine
        x_dim, y_dim, z_dim, num_volumes = brain.shape
        
        # Initialize a matrix of ROI time series and voxel time series

        ROI_matrix = np.zeros((num_ROIs, num_volumes))
        labels = atlas.astype(int) - 1

        # mask_Obj = nib.load(mask_file)
        # mask_data = mask_Obj.get_data()
        '''
        Mask file can also be used here instead of 
        defining a threshold from standard deviation.
        '''
        std_Deviation = np.std(brain, axis = 3)
        brain_voxels_X,brain_voxels_Y,brain_voxels_Z = np.where(std_Deviation>std_param)
        num_brain_voxels = len(brain_voxels_X)

        voxel_matrix = np.zeros((num_brain_voxels, num_volumes))
        # Fill up the voxel_matrix 

        voxel_counter = 0
        num_voxels_in_ROI = np.zeros((num_ROIs,1))
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if std_Deviation[i,j,k] > std_param:
                        voxel_matrix[voxel_counter] = brain[i,j,k] 
                        voxel_counter += 1

                    currentLabel = labels[i,j,k]              
                    if currentLabel >-1:
                        # print(currentLabel)
                        # print(voxel_counter)
                        ROI_matrix[currentLabel,:] +=  brain[i,j,k]
                        num_voxels_in_ROI[currentLabel,0] += 1

        del brain
        del labels, atlas
        print('Total number of selected voxels: ', num_brain_voxels)
        ROImeans = np.divide(ROI_matrix,num_voxels_in_ROI) # Check if divide is working correctly

        Xm = np.subtract(ROImeans, np.expand_dims(np.mean(ROImeans, axis =1),axis=1))
        Ym = np.subtract(voxel_matrix, np.expand_dims(np.mean(voxel_matrix, axis =1),axis =1))

        Xstd = np.expand_dims(np.std(ROImeans, axis = 1) + 1e-7, axis =1)
        Ystd = np.expand_dims(np.std(voxel_matrix,axis = 1) + 1e-7, axis=1)

        mult1 = np.divide(Xm, Xstd)
        mult2 = np.divide(Ym , Ystd)
        coff_matrix = np.dot(mult1, mult2.T)/num_volumes

        Brainimg = np.zeros((x_dim,y_dim,z_dim,num_ROIs),dtype = np.float32)
        for i in range(num_ROIs):
            Brainimg[brain_voxels_X,brain_voxels_Y,brain_voxels_Z,i] = coff_matrix[i] 
        CoffBrain = nib.Nifti1Image(Brainimg, affine = brain_affine)
        print('Saved file in : %s'%coff_matrix_file)
        nib.save(CoffBrain, coff_matrix_file)
    return coff_matrix_file


def pearsonr_with_roi_mean(in_file, atlas_file, mask_file):
    import nibabel as nib
    import numpy as np
    from os.path import join as opj
    import os

    sub_id = in_file.split('/')[-1].split('.')[0]
    fc_file_name = sub_id + '_fc_map.npy'
    coff_matrix_file = opj(os.getcwd(),fc_file_name)
    if os.path.exists(coff_matrix_file):
        print('Saved file in : %s'%coff_matrix_file)
        
    else:

        atlasObject = nib.load(atlas_file)
        atlas = atlasObject.get_data()
        num_ROIs = int(np.max(atlas) - np.min(atlas)) 

        print('Min Index:', np.min(atlas),'Max Index', np.max(atlas))
        print('Total Number of Parcellations = ',num_ROIs)
        brain_data = nib.load(in_file)
        brain = brain_data.get_data()
        brain_affine = brain_data.affine
        x_dim, y_dim, z_dim, num_volumes = brain.shape
        
        # Initialize a matrix of ROI time series and voxel time series

        ROI_matrix = np.zeros((num_ROIs, num_volumes))
        labels = atlas.astype(int) - 1

        mask_Obj = nib.load(mask_file)
        mask_data = mask_Obj.get_data()
        brain_voxels_X,brain_voxels_Y,brain_voxels_Z = np.where(mask_data==1)
        num_brain_voxels = len(brain_voxels_X)

        voxel_matrix = np.zeros((num_brain_voxels, num_volumes))
        # Fill up the voxel_matrix 

        voxel_counter = 0
        num_voxels_in_ROI = np.zeros((num_ROIs,1))
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if mask_data[i,j,k] ==1:
                        voxel_matrix[voxel_counter] = brain[i,j,k] 
                        voxel_counter += 1

                    currentLabel = labels[i,j,k]              
                    if currentLabel >-1:
                        # print(currentLabel)
                        # print(voxel_counter)
                        ROI_matrix[currentLabel,:] +=  brain[i,j,k]
                        num_voxels_in_ROI[currentLabel,0] += 1

        del brain
        del labels, atlas

        ROImeans = np.divide(ROI_matrix,num_voxels_in_ROI) # Check if divide is working correctly

        Xm = np.subtract(ROImeans, np.expand_dims(np.mean(ROImeans, axis =1),axis=1))
        Ym = np.subtract(voxel_matrix, np.expand_dims(np.mean(voxel_matrix, axis =1),axis =1))

        Xstd = np.expand_dims(np.std(ROImeans, axis = 1) + 1e-7, axis =1)
        Ystd = np.expand_dims(np.std(voxel_matrix,axis = 1) + 1e-7, axis=1)

        mult1 = np.divide(Xm, Xstd)
        mult2 = np.divide(Ym , Ystd)
        coff_matrix = np.dot(mult1, mult2.T)/num_volumes

        np.save(fc_file_name, coff_matrix)
        print('Saved file in : %s'%coff_matrix_file)

        Brainimg = np.zeros((x_dim,y_dim,z_dim,num_ROIs),dtype = np.float32)
        for i in range(num_ROIs):
            Brainimg[brain_voxels_X,brain_voxels_Y,brain_voxels_Z,i] = coff_matrix[i] 
        CoffBrain = nib.Nifti1Image(Brainimg, affine = brain_affine)
        fc_file_nii = sub_id + '_fc_map.nii.gz'
        coff_matrix_file_in_nii = opj(os.getcwd(),fc_file_nii)
        print('Saved nii file in : %s'%coff_matrix_file_in_nii)
        nib.save(CoffBrain, coff_matrix_file_in_nii)        
    return coff_matrix_file, coff_matrix_file_in_nii


def pearson_corr_Ankita(in_file, atlas_file):
    import subprocess
    import nibabel as nib
    import os
    import numpy as np
    from os.path import join as opj
    sub_id = in_file.split('/')[-1].split('.')[0]
    coff_matrix_dir = opj(os.getcwd(),sub_id)
    if os.path.exists(coff_matrix_dir):
        shutil.rmtree(coff_matrix_dir)
    atlasObject = nib.load(atlas_file)
    atlas = atlasObject.get_data()
    N_ROIs = int(np.max(atlas) - np.min(atlas))
    # The directory for calling fconn.o needs to be changed later
    # Please change it here accordingly. I shall automate it after 
    # coming back.
    args = ("/home/deepak/Desktop/funcConn/./fconn.o", "-i", in_file, 
            "-r", atlas_file, str(N_ROIs), "-o", coff_matrix_dir)

    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    print(output)
    coff_matrix_file = opj(coff_matrix_dir, 
                            'avg_roi_time_series.nii')
    return coff_matrix_file


def make_npy_from_CorrFile(Corr_file, mask_file):
    import nibabel as nib
    import numpy as np
    from os.path import join as opj
    import os

    sub_id = Corr_file.split('/')[-1].split('.')[0]
    fc_file_name = sub_id + '_fc_map.npy'
    coff_matrix_file = opj(os.getcwd(),fc_file_name)
    if os.path.exists(coff_matrix_file):
        print('Saved file in : %s'%coff_matrix_file)
        
    else:
        brain_data = nib.load(Corr_file)
        brain = brain_data.get_data()

        x_dim, y_dim, z_dim, num_ROIs = brain.shape
        mask_Obj = nib.load(mask_file)
        mask_data = mask_Obj.get_data()
        brain_voxelsX,brain_voxelsY,brain_voxelsZ = np.where(mask_data==1)        
        num_brain_voxels = len(brain_voxelsX)
        Corr_matrix = np.zeros((num_ROIs, num_brain_voxels))
        for i in range(num_ROIs):
            Corr_matrix[i] = brain[brain_voxelsX,brain_voxelsY,brain_voxelsZ,i]
        np.save(fc_file_name, Corr_matrix)
        print('Saved file in : %s'%coff_matrix_file)
    return coff_matrix_file

def build_correlation_wf(Registration = True,
			use_Ankita_Function = False,
				name = 'pearsonCorrcalc'):
    corr_wf = Workflow(name=name)
    if Registration:
        inputnode = Node(interface=util.IdentityInterface(fields=['in_files', 
                                                                   'atlas_files',
                                                                   'func2std',
                                                                   'reference',
                                                                   'mask_file']),
                                                            name='inputspec')
        outputnode = Node(interface=util.IdentityInterface(fields=['pearsonCorr_files']),
                             name='outputspec')
        
        if use_Ankita_Function:
            coff_matrix = MapNode(util.Function(function=pearson_corr_Ankita,
                                    input_names=['in_file','atlas_file'],
                                    output_names=['coff_matrix_file']),
				iterfield = ['in_file','atlas_file'],
                          name = 'coff_matrix')
            transform_corr = MapNode(interface = fsl.ApplyXFM(interp='spline'),
                                iterfield = ['in_file','in_matrix_file'],
                                name='transform_corr')
            maskCorrFile = MapNode(interface=fsl.ImageMaths(suffix='_masked',
                                               op_string='-mas'),
       				iterfield = ['in_file'],
                      name = 'maskWarpFile')
            make_npy_from_Corr = MapNode(util.Function(function=make_npy_from_CorrFile,
                                    input_names=['Corr_file','mask_file'],
                                    output_names=['coff_matrix_file']),
                          iterfield=['Corr_file'],
                          name = 'coff_matrix_in_npy')
            
        else:
            coff_matrix = MapNode(util.Function(function=pearsonr_with_roi_mean_w_reg, 
                                    input_names=['in_file','atlas_file'],
                                    output_names=['coff_matrix_file']),
                          iterfield=['in_file','atlas_file'],
                          name = 'coff_matrix')
            transform_corr = MapNode(interface = fsl.ApplyXFM(interp='spline'), 
                                iterfield = ['in_file','in_matrix_file'],
                                name='transform_corr')
            maskCorrFile = MapNode(interface=fsl.ImageMaths(suffix='_masked',
                                               op_string='-mas'),
                      iterfield=['in_file'],
                      name = 'maskWarpFile')
            make_npy_from_Corr = MapNode(util.Function(function=make_npy_from_CorrFile, 
                                    input_names=['Corr_file','mask_file'],
                                    output_names=['coff_matrix_file']),
                          iterfield=['Corr_file'],
                          name = 'coff_matrix_in_npy')
        datasink = Node(interface=DataSink(), name='datasink')

        corr_wf.connect(inputnode, 'in_files', coff_matrix, 'in_file')
        corr_wf.connect(inputnode, 'atlas_files', coff_matrix, 'atlas_file')
        corr_wf.connect(coff_matrix,'coff_matrix_file', transform_corr, 'in_file')
        corr_wf.connect(inputnode, 'func2std', transform_corr, 'in_matrix_file')
        corr_wf.connect(inputnode, 'reference', transform_corr, 'reference')
        corr_wf.connect(transform_corr,'out_file', maskCorrFile, 'in_file')
        corr_wf.connect(inputnode, 'mask_file', maskCorrFile, 'in_file2')

        corr_wf.connect(maskCorrFile,'out_file', make_npy_from_Corr, 'Corr_file')
        corr_wf.connect(inputnode,'mask_file', make_npy_from_Corr, 'mask_file')
        corr_wf.connect(make_npy_from_Corr, 'coff_matrix_file', outputnode, 'pearsonCorr_files')
        corr_wf.connect(outputnode, 'pearsonCorr_files', datasink, 'out_file')

    else:

        inputnode = Node(interface=util.IdentityInterface(fields=['in_files', 
                                                                   'atlas_file',
                                                                   'mask_file']),
                            name='inputspec')
        outputnode = Node(interface=util.IdentityInterface(fields=['pearsonCorr_files', 
                                                                    'pearsonCorr_files_in_nii']),
                             name='outputspec')
        if use_Ankita_Function:
            coff_matrix = MapNode(util.Function(function=pearson_corr_Ankita,
                                    input_names=['in_file','atlas_file'],
                                    output_names=['coff_matrix_file']),
                                iterfield = ['in_file'],
                          name = 'coff_matrix')
            maskCorrFile = MapNode(interface=fsl.ImageMaths(suffix='_masked',
                                               op_string='-mas'),
                                iterfield = ['in_file'],
                      name = 'maskWarpFile')
            make_npy_from_Corr = MapNode(util.Function(function=make_npy_from_CorrFile,
                                    input_names=['Corr_file','mask_file'],
                                    output_names=['coff_matrix_file']),
                          iterfield=['Corr_file'],
                          name = 'coff_matrix_in_npy')
            datasink = Node(interface=DataSink(), name='datasink')

            corr_wf.connect(inputnode, 'in_files', coff_matrix, 'in_file')
            corr_wf.connect(inputnode, 'atlas_file', coff_matrix, 'atlas_file')
            corr_wf.connect(coff_matrix,'coff_matrix_file',  maskCorrFile, 'in_file')
            corr_wf.connect(inputnode, 'mask_file', maskCorrFile, 'in_file2')

            corr_wf.connect(maskCorrFile,'out_file', make_npy_from_Corr, 'Corr_file')
            corr_wf.connect(inputnode,'mask_file', make_npy_from_Corr, 'mask_file')
            corr_wf.connect(make_npy_from_Corr, 'coff_matrix_file', outputnode, 'pearsonCorr_files')
            corr_wf.connect(outputnode, 'pearsonCorr_files', datasink, 'out_file')           
        else:
            coff_matrix = MapNode(util.Function(function=pearsonr_with_roi_mean, 
                                    input_names=['in_file','atlas_file','mask_file'],
                                    output_names=['coff_matrix_file','coff_matrix_file_in_nii']),
                          iterfield=['in_file'],
                          name = 'coff_matrix')    
            datasink = Node(interface=DataSink(), name='datasink')
            # selectfile = MapNode(interface=util.Select(index=[0]), iterfield = ['inlist'],name='select')
            corr_wf.connect(inputnode, 'in_files', coff_matrix, 'in_file')
            corr_wf.connect(inputnode, 'atlas_file', coff_matrix, 'atlas_file')
            corr_wf.connect(inputnode, 'mask_file', coff_matrix, 'mask_file')

            corr_wf.connect(coff_matrix,'coff_matrix_file', outputnode, 'pearsonCorr_files')
            corr_wf.connect(coff_matrix, 'coff_matrix_file_in_nii', outputnode, 'pearsonCorr_files_in_nii')
            corr_wf.connect(outputnode, 'pearsonCorr_files', datasink, 'out_file')
        # coff_matrix = MapNode(util.Function(function=pearsonr_with_roi_mean_w_reg, 
        #                             input_names=['in_file','atlas_file'],
        #                             output_names=['coff_matrix_file']),
        #                   iterfield=['in_file'],
        #                   name = 'coff_matrix')
        # maskCorrFile = MapNode(interface=fsl.ImageMaths(suffix='_masked',
        #                                        op_string='-mas'),
        #               iterfield=['in_file'],
        #               name = 'maskWarpFile')
        # make_npy_from_Corr = MapNode(util.Function(function=make_npy_from_CorrFile, 
        #                             input_names=['Corr_file','mask_file'],
        #                             output_names=['coff_matrix_file']),
        #                   iterfield=['Corr_file'],
        #                   name = 'coff_matrix_in_npy')
        # datasink = Node(interface=DataSink(), name='datasink')

        # corr_wf.connect(inputnode, 'in_files', coff_matrix, 'in_file')
        # corr_wf.connect(inputnode, 'atlas_file', coff_matrix, 'atlas_file')
        # corr_wf.connect(coff_matrix,'coff_matrix_file', maskCorrFile, 'in_file')
        # corr_wf.connect(inputnode, 'mask_file', maskCorrFile, 'in_file2')

        # corr_wf.connect(maskCorrFile,'out_file', make_npy_from_Corr, 'Corr_file')
        # corr_wf.connect(inputnode,'mask_file', make_npy_from_Corr, 'mask_file')
        # corr_wf.connect(make_npy_from_Corr, 'coff_matrix_file', outputnode, 'pearsonCorr_files')
        # corr_wf.connect(outputnode, 'pearsonCorr_files', datasink, 'out_file')

    return corr_wf
