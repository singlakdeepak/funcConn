import numpy as np
import nibabel as nib
import math
import copy
import sys
# import ipdb

num_thresh = 3
qval = nib.load(sys.argv[1])
qval_np = np.array(qval.dataobj)

def stats(q_thresh):
	#input file is in log
	# ipdb.set_trace()
	q_thresh_log = -1*math.log10(q_thresh)

	qval_cp = copy.deepcopy(qval_np)
	#set values below threshold to 0
	thresh_indices = (np.abs(qval_cp) < q_thresh_log)
	qval_cp[thresh_indices] = 0	

	pos_count = np.sum((qval_cp>0),axis=(0,1,2))
	neg_count = np.sum((qval_cp<0),axis=(0,1,2))

	return pos_count,neg_count	

q_thresh = [0 for i in range(num_thresh)]
#take in from command line
q_thresh[0] = 0.1
q_thresh[1] = 0.05
q_thresh[2] = 0.01

q_max_log = np.amax(qval_np,axis=(0,1,2))
# q_min_pos = np.power(10,-1*q_max_log)

q_min_log  = -1*np.amin(qval_np,axis=(0,1,2))
# q_min_neg = np.power(10,q_min_log)

matrix_list = []
for index in range(num_thresh):
	pos_count,neg_count = stats(q_thresh[index])
	matrix_list.append(pos_count)
	matrix_list.append(neg_count)

matrix_list.append(q_max_log)
matrix_list.append(q_min_log)

# ipdb.set_trace()
output = np.vstack(matrix_list)
output = np.transpose(output)

np.savetxt('ROIfile.csv',output,fmt='%i,%i,%i,%i,%i,%i,%1.3f,%1.3f')


