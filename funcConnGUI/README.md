# Functional Connectivity GUI & Statistics Calculation

## ABOUT
This directory contains the graphical user interface for the toolkit on **Region of Interest based Functional Connectivity Analysis**. The main window looks like this.
![MainWindow](docs/images/MainWindow.png)

For more information about the fields in the main window, please look into [this file](docs/MainWindow.md).
## GUI Architecture

* Made in **Qt Creator**
* Option for Registration of the data is also available if already not registered to the reference image. It used FLIRT for the registration.
* FSL is must for the working of the GUI.
* Accepts the groups upto 5 for now. The statistical analysis for ROI correlation maps has been implemented only for the various combinations of comparison between groups.
* An ROI atlas needs to be supplied where each voxel is labeled from 1 to N where N is total number of ROIs. The voxels outside the brain should be numbered 0.
* Statistical analysis within a group is also allowed. It is assumed that the mean correlation for the population is 0 as the Null Hypothesis. 
* FDR Correction has also been implemented. If you want, you can also check Fisher Transformation for further normalization of the data.
* The statistical tests currently available are:
    * Simple 1 sample Student&rsquo;s t-test
    * 2 sample t-test


## Preprocessing of fMRI data
A user can also preprocess the data if he wants. All the preprocessing methods have been implemented using [NiPype](http://miykael.github.io/nipype-beginner-s-guide/installation.html): a pipeline for neuroimaging analysis. Nipype creates parallel workflows for preprocessing and finding ROI based correlation in different subjects. It also has options for calculating the correlations and statistics without preprocessing or if it has been done already in FSL. The GUI contains the basic preprocessing pipeline with the methods written below:

* **Slice Time Correction**: 
Ascending, Descending or Interleaved order. You can also provide a Slice Time Correction file in the same format as that of FSL. 

* **Motion Correction**: This option corrects the motion of the brain during the scan and calls FSL&rsquo;s MCFLIRT command at the backend.
* **BET Brain Extraction**: This option extracts the brain from the skull and the threshold for brain extraction can be varied from 0 to 1. It also contains the options for doing robust brain extraction. 
* **Spatial Smoothing using SUSAN**: It uses the SUSAN method for spatial smoothing, the same way as implemented in FSL. The value can be varied in mm.
* **Temporal Filtering**: It has options for doing both high pass as well as low pass filtering. The user can also try to do band-pass. The value in GUI is put in seconds while in the json, the value is in units of sigma.
* **Intensity Normalization**: It has the same intensity normalization function as used in FSL. 
