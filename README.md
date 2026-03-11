# EEG-Metastability

MATLAB code for calculating EEG-derived metastability-related measures, including channel-wise **Synchrony Coalition Entropy (SCE)** and a phase-dispersion summary, from preprocessed continuous EEG recordings.

This repository is intended for researchers who want to quantify frequency-specific phase-synchrony dynamics in EEG data using a relatively simple and transparent pipeline.

---

## Overview

This pipeline performs the following steps for each subject and frequency bin:

1. Load preprocessed EEG data from `.mat` files
2. Optionally apply **Current Source Density (CSD)** transformation
3. Compute complex wavelet coefficients using `izmy_gbweeg.m`
4. Extract instantaneous phase
5. Build binary synchrony coalitions using a phase-difference threshold
6. Compute **channel-wise Synchrony Coalition Entropy (SCE)**
7. Save per-subject and per-frequency results

---

## Features

- Frequency-resolved analysis
- Channel-wise SCE output
- Optional CSD preprocessing
- Batch processing of multiple subjects
- Parallel processing support (`parfor`)
- Structured output suitable for downstream statistics

---

## Requirements

### MATLAB
Tested in MATLAB with standard numeric and parallel computing functionality.

### Required external dependencies
This repository does **not** bundle all third-party dependencies. You need:

- **EEGLAB**  
  Required for `eegfilt`

- **CSD toolbox / CSDconvert function**  
  Required only if `UseCSD = true`

- **CSD basis file**  
  Example parameter: `CSDbasis.mat`

### Included in this repository
- `calcMSISCE.m`
- `izmy_gbweeg.m`

---

## Input data format

Each input file must be a `.mat` file containing an EEG structure:

```matlab
EEG.data
