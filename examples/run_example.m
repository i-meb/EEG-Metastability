% run_example  Minimal executable example for EEG-Metastability.
%
% This script:
%   1. Locates the repository root.
%   2. Adds the repository functions to the MATLAB path.
%   3. Calls calcMSISCE, which checks that EEGLAB/eegfilt is available.
%   4. Loads examples/sub_001.mat.
%   5. Runs calcMSISCE using a small frequency range.
%   6. Displays the resulting MSI and mean SCE.
%
% SPDX-License-Identifier: BSD-3-Clause
% Copyright (c) 2026 Mebuki Izumiya

if exist('eegfilt', 'file') ~= 2
    error(['EEGLAB function eegfilt.m was not found. ' ...
           'Start EEGLAB or add EEGLAB to the MATLAB path before running this example.']);
end

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
    'SaveResults', true, ...
    'Verbose', true);

% If you want to work with multiple files:
results = calcMSISCE(scriptDir, ...
    'FilePattern', 'sub_*.mat', ...
    'SampleRate', 1000, ...
    'Channels', 63, ...
    'FrequencyRange', 1:30, ...
    'Threshold', 1.2, ...
    'SaveResults', true, ...
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
