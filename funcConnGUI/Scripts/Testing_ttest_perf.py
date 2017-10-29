import os
import numpy as np
import timeit
import ttest

start = timeit.default_timer()
list1 = os.listdir('Group1/')
list1 = [os.path.join(os.getcwd(), 'Group1/' + file) for file in list1]

list2 = os.listdir('Group2/')
list2 = [os.path.join(os.getcwd(), 'Group2/'+ file) for file in list2]

Tvals, Pvals = ttest.ttest_1samp_for_all_ROIs(list1,'MNI152_T1_2mm_brain_mask.nii.gz')

stop = timeit.default_timer()
print Tvals

print Pvals

print 'Time taken = %s'%(stop - start)
