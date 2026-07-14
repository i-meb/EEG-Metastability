% make_example_eeg:  Generate a synthetic continuous EEG dataset for testing.
%
%   This script generates a reproducible synthetic EEG dataset compatible
%   with calcMSISCE.m. The generated data are intended solely for testing,
%   demonstration, and software validation. 
%
% GENERATED DATA
%   EEG.data
%       Synthetic continuous EEG signals with dimensions:
%
%           [channels x time samples] = [63 x 10000]
%
%   The dataset has the following properties:
%
%       Number of channels:       63
%       Sampling frequency:       1000 Hz
%       Number of time samples:   10000
%       Signal duration:          10 seconds
%
% SIGNAL MODEL
%   Each channel is generated independently as:
%
%       x_ch(t) = sin(2*pi*10*t + phi1_ch) ...
%               + 0.5*sin(2*pi*18*t + phi2_ch) ...
%               + 0.3*epsilon_ch(t)
%
%   where:
%
%       phi1_ch, phi2_ch
%           Independent random phases drawn uniformly from [0, 2*pi).
%
%       epsilon_ch(t)
%           Independent Gaussian white noise with zero mean and unit
%           variance.
%
%   The random-number generator is initialized using RNG(1), allowing the
%   same dataset to be reproduced across runs under compatible MATLAB
%   versions.
%
% OUTPUT FILE
%   The script saves the generated structure as:
%
%       examples/dummy_001.mat
%
%   The output MAT-file contains one variable:
%
%       EEG
%
%   with the signal stored in:
%
%       EEG.data
%
%   The directory "example_data" must exist before the current
%   implementation is executed.
%
% USAGE
%   Run this script from MATLAB:
%
%       run('examples/make_example_eeg.m')
%
%   The resulting file can then be analyzed using:
%
%       results = calcMSISCE('example_data/sub_001.mat', ...
%           'SampleRate', 1000, ...
%           'Channels', 63);
%
% NOTES
%   - This dataset is synthetic and contains no personal or participant
%     information.
%   - It is not designed to reproduce realistic EEG spatial covariance,
%     volume conduction, nonstationarity, artifacts, or physiological
%     spectral structure.
%   - Random phases are generated independently for each channel.
%   - The dataset should therefore be used only to confirm that the
%     analysis pipeline executes correctly.
%
% COPYRIGHT
%   Copyright (c) 2026 Mebuki Izumiya
%
% LICENSE
%   SPDX-License-Identifier: BSD-3-Clause

rng(1);

nCh = 63;
nT = 10000;
fs = 1000;
t = (0:nT-1) / fs;

EEG = struct();
EEG.data = zeros(nCh, nT);

for ch = 1:nCh
    phase1 = 2 * pi * rand();
    phase2 = 2 * pi * rand();
    EEG.data(ch,:) = ...
        sin(2*pi*10*t + phase1) + ...
        0.5 * sin(2*pi*18*t + phase2) + ...
        0.3 * randn(1, nT);
end

save(fullfile('examples', 'dummy_001.mat'), 'EEG');
disp('Example dataset saved to examples/dummy_001.mat .');
