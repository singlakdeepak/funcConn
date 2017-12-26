import ipdb
import nibabel as nib
import numpy as np
import math
import copy
import pandas as pd
from pandas import DataFrame as df
import sys

QFile = sys.argv[1]
ROIfile = sys.argv[2]
qthresh = float(sys.argv[3])
Cthresh = float(sys.argv[4])


def convertROItoBin(ROImask_np):
	xdim,ydim,zdim = ROImask_np.shape
	num_roi = int(np.max(ROImask_np))

	subtractor = np.ones([xdim,ydim,zdim])
	ROI_mask_bin = np.zeros([xdim,ydim,zdim,num_roi])
	
	for i in range(num_roi):
		xcoord,ycoord,zcoord = np.where(ROImask_np == 1.0)
		ROI_mask_bin[xcoord,ycoord,zcoord,i] = 1.0
		ROImask_np = ROImask_np - subtractor

	return ROI_mask_bin

qval = nib.load(QFile)
qval_np = np.array(qval.dataobj)

ROImask = nib.load(ROIfile)
ROImask_np = np.array(ROImask.dataobj)

ROImask_bin = convertROItoBin(ROImask_np)

# qthresh = float(sys.argv[3])
qthresh_log = -1*math.log10(qthresh)

qval_np[np.abs(qval_np) < qthresh_log] = 0.0

num_roi_1 = qval_np.shape[3]
num_roi_2 = ROImask_bin.shape[3]

output=[]

for i in range(num_roi_1):
# for i in range(2):
	# # X_roi,Y_roi,Z_roi = np.where(ROImask_np == (i+1))
	# qval_pos = np.sum(np.logical_and((ROImask_np == (i+1)),(qval_np >= qthresh_log)),axis=(0,1,2))
	# qval_neg = np.sum(np.logical_and((ROImask_np == (i+1)),(qval_np <= -1.0*qthresh_log)),axis = (0,1,2))
	# roi_size = np.sum(ROImask_np == (i+1),axis=(0,1,2))

	# qval_pos_per = float(float(qval_pos)/float(roi_size))
	# qval_neg_per = float(float(qval_neg)/float(roi_size))

	# stat_matrix.append(qval_pos_per)
	# stat_matrix.append(qval_neg_per)
	# stat_matrix.append(roi_size)


	qval_i = qval_np[:,:,:,i]
	# qval_j_neg = qval[:,:,:,i]
	
	qval_rep = np.repeat(qval_i[:,:,:,np.newaxis],num_roi_2,axis=3)
	
	qval_pos = np.sum(np.logical_and((qval_rep>0),(ROImask_bin==1.0)),axis=(0,1,2))
	qval_neg = np.sum(np.logical_and((qval_rep<0),(ROImask_bin==1.0)),axis=(0,1,2))
	roi_size = np.sum(ROImask_bin,axis=(0,1,2))

	# ipdb.set_trace()
	qval_pos_per = np.divide(qval_pos,roi_size)
	qval_neg_per = np.divide(qval_neg,roi_size)
	roi_a_no = i*np.ones(num_roi_2)
	roi_b_no = np.array(range(num_roi_2))

	comb = np.vstack([roi_a_no,roi_b_no,qval_pos_per,qval_neg_per])
	comb_t = np.transpose(comb)

	output.append(comb_t)
	print("Done with ROI",i)

output_combined = np.concatenate(output)

output_df = df(output_combined)
result_df = output_df[(output_df[2]>=Cthresh)|(output_df[3]>=Cthresh)]

result_df.to_csv('ROI_connectivity_'+str(qthresh)+'_'+str(Cthresh))
# ipdb.set_trace()







# stanp.vstack()









