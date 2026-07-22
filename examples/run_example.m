% run_example  Minimal executable example for EEG-Metastability.
%
% This script analyzes the included sample EEG file at selected
% frequencies, saves the results, and plots MSI and mean SCE.
%
% License: BSD-3-Clause
% Copyright (c) 2026 Mebuki Izumiya

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);

addpath(repoRoot);

if exist('eegfilt', 'file') ~= 2
    error(['EEGLAB function eegfilt.m was not found. ', ...
        'Start EEGLAB or add EEGLAB to the MATLAB path.']);
end

exampleFile = fullfile(scriptDir, 'sub_001.mat');
outputDir   = fullfile(repoRoot, 'output');

results = calcMSISCE(exampleFile, ...
    'SampleRate', 1000, ...
    'Channels', 63, ...
    'FrequencyRange', [10, 18], ...
    'Threshold', 1.2, ...
    'WaveletCycles', 1, ...
    'UseParallel', false, ...
    'SaveResults', true, ...
    'OutputDir', outputDir, ...
    'OutputFileName', 'example_results.mat', ...
    'Verbose', true);

figure;
plot(results.frequencies, squeeze(results.MSI(1,1,:)), ...
    '-o', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('MSI');
title('Example MSI');
grid on;

figure;
plot(results.frequencies, squeeze(results.meanSCE(1,1,:)), ...
    '-o', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Mean normalized SCE');
title('Example mean SCE');
grid on;
