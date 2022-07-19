# RS-FetMRI
RS-FetfMRI for processing Fetal resting-state functional MRI scans

### Introduction

Resting-state functional magnetic resonance imaging (rs-fMRI) has most recently proved to open a measureless window on functional neurodevelopment in utero. We focused the Fetal resting state functional MRI processing pipeline (RS-FetMRI)  development on creating an effective, user-friendly and easy to use package completely integrated in SPM (Statistical Parametric Mapping).  This processing pipeline is suitable for both **single subject** and **group-based analyses** and can be used for both **1.5T** and **3T scanner** acquisitions.This pipeline does not need any structural scans (**structural-free**).
The RS-FetMRI is a semi-automatic and standardized pipeline composed of six (M1 to M6) modules for processing fetal resting-state functional MRI (Fetal rs-fMRI) **Nifti data**. While the first three modules, from M1 to M3, work Within Session (WS) the last four, from M4 to M6, work Between Session (BS) (Figure 1).
The RS-FetMR is capable of 1) detecting and correcting fetal-specific motion effects and signal intensity changes, especially through 2) accuracy of time-series spatial normalization to a standardized gestational-week specific fetal template space via the synergetic action of each module. Furthermore, the whole processing does not need a structural fetal scan, which can be difficult to acquire/process. This RS-FetMRI protocol is suitable for a large pool of users, from **beginners** to **experts**, although a basic technical knowledge of fetal functional image processing is required. A detailed explanation of how to deal with this pipeline is presented in the **User-Manual** in the RS-FetMRI folder.

<p align="center">
<img src="https://github.com/NicoloPecco/RS-FetMRI/blob/main/Images/Flowchart.png" width="800" height="800">
</p>
<p align="center">
Flowchart of the RS-FetMRI pipeline.
</p>

# Usage 

The entire RS-FetMRI package is composed of a main script, called ‘RSFetfMRI.m’, which is subdivided into six modules (Figure 1), a MATLAB custom-built function ‘Create_template.m’ and an automatic modified version of the Artifact Detection Tool (https://www.nitrc.org/projects/artifact_detect) called ‘art.m’ accompanied with relative configuration files. The RS-FetMRI script uses a wide range of  SPM functions for image processing and image visualization. It also contains a folder called ‘Templates’ that includes three subfolders (‘Template_for_session’, ‘Template_Priors_Seg’ and ‘Template_orig’) containing different Nifti file useful for the 1st-pass masking, segmentation and visualization. Particularly, as explained in the ‘Installation and Requirements’ paragraph, the ‘Template_orig’ files need to be downloaded from a different website.

# Requirement

The RS-FetMRI package can be downloaded from GitHub (https://github/NicoloPecco/RS-FetMRI). The RS-FetMRI  package can be used on any computer operating system with installed:

- Matlab2013 (or above);
- SPM12 (or above);
- The user must download the CRL Fetal Brain Atlas images  – GW 21 through 37– from <http://crl.med.harvard.edu/research/fetal_brain_atlas/> and then unzip the 'STA*.nii.gz' file (i.e. STA21.nii.gz --> STA21.nii or STA37exp.nii.gz --> STA37.nii). These file must be inserted in the subfolder ‘Template_orig’ which is found within the ‘Templates’ folder;
- Initialization.

# Initialization

We strongly encourage the users to read the user manual and to see the 'Tutorial_RSFetMRI.mp4' viedo, expecially for the initialization step, before to run any analyses.  
Enter the ‘art’ folder within the ‘RS-FetMRI’ package. This folder contains **three configuration files** (‘.cfg’). Within each of the configuration file, **three paths (‘image_dir’, ‘motion_dir’ and ‘mask_file’) need to be modified** before any analysis can be exploited (See supplementary Initialization in the User Manual).
If the ART toolbox was previously installed, it should be removed from the MATLAB paths and then updated with the new script path instead.
Once the configuration files were modified, **open Matlab** software and **set the current folder** within the ‘RS-FetMRI’ folder and press the run button.

# TestSet

The RS-FetMRI package also includes a complete Test Set. The subject under examination was at the 35 GW.All the images present in this manual were selected from this subject. If the user would like to try the analysis on this subject, simply download the material and then the user should move all of the folders regarding M2 to M6 into a different repository. Please note that initialization is required to try the TestSet. All of the parameters inserted during the processing are reported in the supplementary TestSet parameters. 


### Final Data Movie
SPM movie display of the final Dataset. This movie is related to the TestSet. See Manual.
<p align="center">
<img src="https://github.com/NicoloPecco/RS-FetMRI/blob/main/Images/Screen%20Recording%202021-07-06%20at%2013.02.57.gif">
</p>

# Citation

Pecco N, Canini M, Mosser KHH, Caglioni M, Scifo P, Castellano A, Cavoretto P, Candiani M, Baldoli C, Falini A, Rosa PAD. RS-FetMRI: a MATLAB-SPM Based Tool for Pre-processing Fetal Resting-State fMRI Data. Neuroinformatics. 2022 Jul 14. doi: 10.1007/s12021-022-09592-5. Epub ahead of print. PMID: 35834105.
