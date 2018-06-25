# Statistics

The stats window looks like this

![Stats](images/statsWindow.png)

Currently, a user can perform a simple one sample or two sample Student&rsquo;s t-test. By selecting,

* **Within Groups**, one can perform a one sampled t-test with zero group mean and the population standard deviation (estimated from the number of subjects).
* **Between Groups**, one can perform a two-sampled t-test with either of the two alternate hypothesis, **G1>G2** or **G2>G1**. Another option to choose the combination of groups for comparison is also available here. 

The normalization method currently available is Fisher method (introduced by Ronald Fisher).

While doing False Discovery Rate (FDR) correction, a user can choose, either to perform all the hypothesis tests together (**Joint FDR Correction**) or treat the N region of interest (ROI) maps independently (**Separate FDR correction for each ROI**), where N is the number of ROIs. The first method largely makes the corrected maps insignificant, thus the results might not represent the truth. If a large number of data samples are available, then joint FDR correction can be used for the analysis.