# High Throughput Optical Physiology
Code companion to https://www.biorxiv.org/content/10.1101/2024.12.23.629904v1. Includes _peakFinder.R_, an automated peak detection algorithm for fluorescence transient analysis. Copyright 2024, University of Maryland, Baltimore. All rights reserved. Patent pending.


**What is the purpose of this repository?** <br/>
The goal of this repository is leverage the intensity-based glutamate sensing fluorescent reporter (iGluSnFR3) for synaptic physiology, at scale. <br/>  
We achieve this with three key innovations. <br/>

1) Automated data acquisition by coordinating hardware via Python.
2) Automated segmentation of synapses via iGluSnFR3 activity or a fluorescent marker (using SynQuant).
3) Automated correction and normalization of intensity-time traces and extraction of fluorescence transients for downstream analysis.

NOTE: Our approach is also compatible with calcium imaging, but this may require additional development. We will strive to be responsive to users to help enable the pipeline for other fluorescent sensors. <br/>

**What is contained in this repository?**
1) _Automated segmentation via activity or marker_:  <br/>
A pipeline of macros in ImageJ macro language and MATLAB. See iGlu_activitySegmentation/ and iGlu_markerSegmentation/.  <br/>

2) _Automated trace normalization and peak extraction_:  <br/>
A standalone _R_ script which accepts an arbitrary number of .csv files across an arbitrary number of subdirectories. So long as the traces can be discriminated on the basis of sufficient string variables, _peakFinder.R_ automatically extracts iGluSnFR3 transients and normalizes traces via an iterative baseline identification algorithm predicated on a median filter. _peakFinder.R_ is also viable for calcium sensors, but performs baseline correction using a percentile filter. This repository will be actively maintained according to the needs of the user base and we will strive to be responsive to user inquiries.  <br/>  

3) _Demo scripts to produce the Figures in the manuscript_:  <br/>
_R_ has a number of packages which make publication-quality data visualization very accessible. All data in this manuscript was analyzed and visualized using custom _R_ scripts. We provide access to the original datasets and the bespoke analysis scripts to produce all the figures in the manuscript. As advanced users will observe, there are a number of sub-functions which carry out small actions on the dataset. The terminal scripts for each "peak annotation module" are master scripts which manipulate the dataset, extract summary statistics, and plot the data. If you appreciate the visualizations in this manuscript and are feeling brave, we invite you to inspect the scripts and take the visualizations you like. We anticipate that certain scripts, like _tracePlotter.R_, will be handy for a variety of users for getting high-quality graphs of activity for each ROI and debugging.  <br/>

**How do I use the code in this repository?**
1) _Automated segmentation_:  <br/>
To use activity-based segmentation, you will need a MATLAB instance, as the code for frequency-domain filtering of iGluSnFR recordings was originally written by Phillipe Mendonça in MATLAB (see Mendonça, et al. _Nature Communications_, 2021). 
To use marker-based segmentation, you will need to install SynQuant and Fast4DReg in ImageJ (or FIJI).

2) _Automated trace normalization and peak extraction_:  <br/>
To use _peakFinder.R_, you will need to download R and RStudio (https://posit.co/download/rstudio-desktop/). The code takes advantage of a few dozen packages for data manipulation, curve-fitting, and data visualization which you can install in RStudio. The command for installing these packages is written into _run_peakFinder.R_.  

3) _Demo scripts to produce the Figures in the manuscript_:  <br/>
Follow the same instructions as for _peakFinder.R_.


**How long will it take me to get started using this code repository?**
Downloading and installing R, RStudio, and installing dependcy libraries: <br/>
Running _peakFinder.R_ on the demo dataset: 2 minutes <br/>
Running the analysis to generate Figure 5: 3.5 minutes <br/>
Running the analysis to generate Figure 6: 7.5 minutes <br/>

**Help! I can't find the demo datasets!** <br/>
_Data for peakFinder.R_: https://drive.google.com/drive/folders/1_Mm8QNhFO6zKyGGvtbehgSkhrEo5kW_H?usp=drive_link  <br/>
_Data for Figure 5_: https://drive.google.com/drive/folders/1ENJNOxWQrlXi2OHbr8fbtA3porr2fKp0?usp=drive_link <br/>
_Data for Figure 6_: https://drive.google.com/drive/folders/19ac_NfB_gWiH5fPiM7eQ7aH5aw4GKAwd?usp=drive_link <br/>

**COMING SOON**
1) Detailed README files to help users leverage the full pipeline.
2) _fileManipulate.R_, for scalable manipulation of directories and filenames.
3) Demo scripts for Figure 2, 3, 4.


