# EEG-Metastability

# This script calculates the MSI (Metastability Index) and SCE (Synchrony Coalition Entropy) from EEG data for Linux. 

It performs the following steps:
  1. Load EEG data for multiple subjects.
  2. Apply Current Source Density (CSD) transformation to preprocess EEG data.
  3. Perform Wavelet Transform to calculate the instantaneous phase of the EEG signals.
  4. Calculate MSI, which quantifies phase dispersion across electrodes.
  5. Calculate SCE, which measures phase synchrony patterns between electrodes and over frequency bands.

Terms:
 - EEG data to be loaded is after noise rejection.
 - CSD: A preprocessing technique to highlight local activity and reduce volume conduction effects.
 - Wavelet Transform: A method to decompose a signal into frequency components, used here to calculate the phase.
 - MSI: A measure of the variability of phase differences across electrodes.
 - SCE: A measure of phase synchrony patterns across electrodes and their entropy.

Data:
 - EEG data is loaded from files named 'sub_*_r.mat'.
 - Results are saved as .mat files in the specified output directory.
