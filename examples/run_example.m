% run_example  Minimal executable example for EEG-Metastability.
%
% This script:
%   1. Locates the repository root.
%   2. Adds the repository functions to the MATLAB path.
%   3. Checks that EEGLAB/eegfilt is available.
%   4. Loads examples/sub_001.mat.
%   5. Runs calcMSISCE using a small frequency range.
%   6. Displays the resulting MSI and mean SCE.
%
% The example dataset is synthetic and contains no human participant data.
%
% SPDX-License-Identifier: BSD-3-Clause
% Copyright (c) 2026 Mebuki Izumiya


scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(repoRoot);
exampleFile = fullfile(scriptDir, 'sub_001.mat');

% If you want to work with a single sample files:
results = calcMSISCE(exampleFile, ...
    'FilePattern', 'sub_*.mat', ...
    'SampleRate', 1000, ...
    'Channels', 63, ...
    'FrequencyRange', 1:30, ...
    'Threshold', 1.2, ...
    'SaveResults', false, ...
    'Verbose', true);

% If you want to work with multiple files:
results = calcMSISCE(scriptDir, ...
    'FilePattern', 'sub_*.mat', ...
    'SampleRate', 1000, ...
    'Channels', 63, ...
    'FrequencyRange', 1:30, ...
    'Threshold', 1.2, ...
    'SaveResults', false, ...
    'Verbose', true);

% Visualization

figure;
plot(results.frequencies, squeeze(results.MSI(1,1,:)), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('MSI');
title('Example MSI spectrum');

figure;
plot(results.frequencies, squeeze(results.meanSCE(1,1,:)), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Mean normalized SCE');
title('Example SCE spectrum');
