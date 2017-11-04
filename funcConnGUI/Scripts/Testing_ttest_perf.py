import os
import numpy as np
import timeit
import ttest
from numpy import ma
start = timeit.default_timer()
#list1 = os.listdir('Group1/')
#list1 = [os.path.join(os.getcwd(), 'Group1/' + file) for file in list1]

#list2 = os.listdir('Group2/')
#list2 = [os.path.join(os.getcwd(), 'Group2/'+ file) for file in list2]
list1 = np.load('Grp1.npy')
list2 = np.load('Grp2.npy')
#Tvals, Pvals = ttest.ttest_1samp_for_all_ROIs(list1,'MNI152_T1_2mm_brain_mask.nii.gz')
Tvals, Pvals = ttest.ttest_ind_samples_if_npy(list1,list2,'MNI152_T1_2mm_brain.nii.gz')
Tvals, Pvals = ma.filled(Tvals), ma.filled(Pvals)
stop = timeit.default_timer()
print( Tvals)

print (Pvals)

np.save('Tvals.npy',Tvals)
np.save('Pvals.npy',Pvals)
print ('Time taken = %s'%(stop - start))
