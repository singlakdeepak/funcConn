from os.path import join as opj
import nibabel as nib
import numpy as np
import nipype.interfaces.utility as util
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node, MapNode
import os


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
        return

    atlasObject = nib.load(atlas_file)
    atlas = atlasObject.get_data()
    num_ROIs = int(np.max(atlas) - np.min(atlas)) 

    print('Min Index:', np.min(atlas),'Max Index', np.max(atlas))
    print('Total Number of Parcellations = ',num_ROIs)
    brain_data = nib.load(in_file)
    brain = brain_data.get_data()

    x_dim, y_dim, z_dim, num_volumes = brain.shape
    
    # Initialize a matrix of ROI time series and voxel time series

    ROI_matrix = np.zeros((num_ROIs, num_volumes))
    labels = atlas.astype(int) - 1

    mask_Obj = nib.load(mask_file)
    mask_data = mask_Obj.get_data()
    brain_voxels,_,_ = np.where(mask_data==1)
    num_brain_voxels = len(brain_voxels)

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
    return coff_matrix_file


def build_correlation_wf(name = 'pearsonCorrcalc'):
    corr_wf = Workflow(name=name)

    inputnode = Node(interface=util.IdentityInterface(fields=['in_files', 
                                                               'atlas_file',
                                                               'mask_file']),
                        name='inputspec')
    outputnode = Node(interface=util.IdentityInterface(fields=['pearsonCorr_files']),
                         name='outputspec')
    coff_matrix = MapNode(util.Function(function=pearsonr_with_roi_mean, 
                                input_names=['in_file','atlas_file','mask_file'],
                                output_names=['coff_matrix_file']),
                      iterfield=['in_file'],
                      name = 'coff_matrix')    
    datasink = Node(interface=DataSink(), name='datasink')

    corr_wf.connect(inputnode, 'in_files', coff_matrix, 'in_file')
    corr_wf.connect(inputnode, 'atlas_file', coff_matrix, 'atlas_file')
    corr_wf.connect(inputnode, 'mask_file', coff_matrix, 'mask_file')

    corr_wf.connect(coff_matrix,'coff_matrix_file', outputnode, 'pearsonCorr_files')
    corr_wf.connect(outputnode, 'pearsonCorr_files', datasink, 'out_file')

    return corr_wf