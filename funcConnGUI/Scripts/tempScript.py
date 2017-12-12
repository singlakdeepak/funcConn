import nibabel as nib
import numpy as np
import os
from scipy import stats
from numpy import ma
import scipy.special as special
from statsmodels.stats import multitest
from multiprocessing import Pool


def make_brain_back_from_npy(NumpyfileList,FileListNames,mask_file):
    maskObj = nib.load(mask_file)
    maskhd = maskObj.get_header()
    maskAffine = maskObj.affine
    maskhd['data_type'] = 'FLOAT32'
    maskData = maskObj.get_data()
    maskData = maskData.astype(int)
    X,Y,Z = np.where(maskData==1)
    x_dim, y_dim, z_dim = maskData.shape
    print(x_dim)
    print(y_dim)
    print(z_dim)
    print(maskData.shape)
    ROIs = NumpyfileList[0].shape[0]
    Brainimg = np.zeros((x_dim,y_dim,z_dim,ROIs))
    print(maskData.shape)
    for i in range(len(NumpyfileList)):
        Map = NumpyfileList[i]
        for j in range(ROIs):
            ThisImg = Brainimg[:,:,:,j]
            ThisImg[X[:228483],Y[:228483],Z[:228483]] = Map[j]
            Brainimg[:,:,:,j] = ThisImg
        Brain = nib.Nifti1Image(Brainimg, affine = maskAffine, header = maskhd)
        nib.save(Brain, FileListNames[i])
