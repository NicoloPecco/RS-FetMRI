# RS-FetfMRI
RS-FetfMRI for processing Fetal resting-state functional MRI scans

### Introduction

Resting-state functional magnetic resonance imaging (rs-fMRI) has most recently proved to open a measureless window on functional neurodevelopment in utero. We focused the Fetal resting state functional MRI processing pipeline (RS-FetMRI)  development on creating an effective, user-friendly and easy to use package completely integrated in SPM (Statistical Parametric Mapping).  This processing pipeline is suitable for both **single subject** and **group-based analyses** and can be used for both **1.5T** and **3T scanner** acquisitions.This pipeline does not need any structural scans (**structural-free**).
The RS-FetMRI is a semi-automatic and standardized pipeline composed of six (M1 to M6) modules for processing fetal resting-state functional MRI (Fetal rs-fMRI) **Nifti data**. While the first three modules, from M1 to M3, work Within Session (WS) the last four, from M4 to M6, work Between Session (BS) (Figure 1).
The RS-FetMR is capable of 1) detecting and correcting fetal-specific motion effects and signal intensity changes, especially through 2) accuracy of time-series spatial normalization to a standardized gestational-week specific fetal template space via the synergetic action of each module. Furthermore, the whole processing does not need a structural fetal scan, which can be difficult to acquire/process. This RS-FetMRI protocol is suitable for a large pool of users, from **beginners** to **experts**, although a basic technical knowledge of fetal functional image processing is required. A detailed explanation of how to deal with this pipeline is presented in the **User-Manual** in the RS-FetMRI folder.

![](https://github.com/NicoloPecco/RS-FetMRI/blob/main/Images/Flowchart.png)
> Flowchart of the RS-FetMRI pipeline.

### Usage ...Continued

# END
