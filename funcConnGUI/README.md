# Functional Connectivity GUI & Statistics Calculation

## GUI Architecture

* Made in **Qt Creator**
* Option for Registration of the data is also available if already not registered to the reference image. It used FLIRT for the registration.
* FSL is must for the working of the GUI.
* Accepts the groups upto 5 for now because the statistical analysis for ROI correlation maps has been done only for difference between two groups for now.
* An ROI atlas needs to be supplied where each voxel is labeled from 1 to N where N is total number of ROIs. The voxels outside the brain should be numbered 0.
* Statistical analysis within a group is also allowed. It is assumed that the mean correlation for the population is 0 as the Null Hypothesis. 
* FDR Correction has also been implemented. If you want, you can also check Fisher Transformation for further normalization of the data.
* The statistical tests currently available are:
    * Simple 1 sample ttest
    * 2 sample t-test
    * Paired t-test

### Preprocessing of fMRI data
* User can also preprocess the data accordingly. It also has options for proceeding without preprocessing or if it has been done already in FSL. By the way, here the GUI contains the basic preprocessing pipeline with the methods:
    * Slice Timing Correction
    * Motion Correction
    * BET Brain Extraction
    * Spatial Smoothing using SUSAN
    * Temporal filtering with both High Pass as well as low pass or BandPass
    * Intensity Normalization
* All the Preprocessing methods have been implemented using **NiPype**:Pipeline for Neuroimaging analysis.

