# ued-tools
Matlab code for processing and interpreting data from Ultrafast Electron Diffraction (UED) experiments

## 1. Introduction
UED is a type of pump-probe experiment where a sample of interest is excited by a laser pulse (pump) and probed with a similarly short high-energy electron pulse.  The time delay between the pump laser and the electron probe is varied precisely and the experiment is repeated many times to produce a series of diffraction patterns with a range of time delays between the pump and probe. These experiments typically have low signal-to-noise ratios, so the measurement is typically performed many times at each specific timepoint and the resulting data is typically a 4d dataset with dimensions [nx, ny, nt, nframes], where nx,ny are the number of pixels on the detector, nt is the number of timepoints measured and nframes is the number of nframes at each time point. This is typically larger than t
