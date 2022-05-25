# ued-tools
Matlab code for processing and interpreting data from Ultrafast Electron Diffraction (UED) experiments

# Introduction
UED is a type of pump-probe experiment where a sample of interest is excited by a laser pulse (pump) and probed with a similarly short high-energy electron pulse.  The time delay between the pump laser and the electron probe is varied precisely and the experiment is repeated many times to produce a series of diffraction patterns with a range of time delays between the pump and probe. This can be thought of as a reciprocal- or Fourier-space movie of the structural dynamics, however, for chemical insight into the atomic-level motion we need to develop data analysis protocols. This can be divided into two main steps: data reduction and structure fitting.

# 1. Data reduction
These experiments typically have low signal-to-noise ratios, so the measurement is performed many times at each specific timepoint rendering a 4d dataset with dimensions `[nx, ny, nt, nframes]` (where nx and ny are the dimensions in pixels of the detector, nt is the number of timepoints measured, and nframes is the number of frames at each time point). This is typically larger than the RAM of a desktop computer, but contains much redundancy of information. For structural analysis what is required is the diffracted intensity for each region of interest corresponding to one or more (h,k,l) indicies as a function of time, ie a dataset with dimensions `[nroi, nt]`. The following steps are necessary to reduce the full dataset to this form:
* Masking problematic detector regions
* Outlier removal due to rare experimental errors
* ROI detection for useful data
* Averaging and uncertainty estimation
* Image distortion correction

Functions and classes for these steps are located in ./lib/Data_Processing
The most important classes are the following:
* ImageProcessingClass2 - Tool for creating a beam-stop mask, locating diffraction spots, and generating ROIs
* ScanProcessingClass2/3 - Takes an instance of ImageProcessingClass2 and experimental data and returns a reduced data set
* ImageDistortionCalculator - Calculates the distortion of a diffraction image based on an ideal image

# 2. Structure refinement
In this stage we want to infer the most likely trajectory of the atoms from the reduced diffraction data. This involves the following steps:
* Indexing: Assigning contributions of (h,k,l) indices to ROIs
* Generating a model based on prior knowledge of the system
* Optimization of atom positions to mimimize loss function (disagreement between experimental and calculated intensities and prior function)

Functions and classes for these steps are located in ./lib/Ediff_functions
The most important of these are
* refineStructure_v14 - Performs the optimization using experimental data, the crystal structure and the model as input
* DiffractionCalculator6  - A tool for calculating realistic electron diffraction patterns.  This can be used for indexing or cross-checking, but is rather slow for fitting purposes


# Example usage
Livescripts showing the correct usage of the functions and classes above will be provided in the examples folder


### Note
The files contained in Crystallography_fileIO are not my own. See licence.txt for credits.