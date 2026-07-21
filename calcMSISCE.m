function results = calcMSISCE(inputPath, varargin)
% calcMSISCE  Compute EEG Metastability Index and Synchrony Coalition Entropy.
%
%   RESULTS = calcMSISCE(INPUTPATH) calculates frequency-resolved
%   Metastability Index (MSI) and channel-wise Synchrony Coalition Entropy
%   (SCE) from preprocessed continuous EEG data.
%
%   -----
%   INPUT
%
%       Path to either:
%       (1) a directory containing MATLAB EEG files, or
%       (2) a single MATLAB EEG file.
%
%       By default, each file must contain:
%
%           EEG.data
%
%       with dimensions:
%
%           [channels x time samples]
%
%   ---------------------
%   NAME-VALUE PARAMETERS
%
%   'FilePattern'       Filename pattern. Default: 'sub_*.mat'
%   'OutputDir'         Directory for saved results. Default: ''
%   'SaveResults'       Save output to disk. Default: false
%   'OutputFileName'    Output MAT filename.
%                       Default: 'msi_sce_results.mat'
%   'DataVariable'      Top-level MAT variable. Default: 'EEG'
%   'DataField'         Signal field in DataVariable. Default: 'data'
%   'SampleRate'        Sampling frequency in Hz. Default: 1000
%   'Channels'          Number of channels analyzed. Default: 63
%   'TimeIndices'       Time-sample indices. Default: [] (all samples)
%   'FrequencyRange'    Wavelet center frequencies in Hz. Default: 1:47
%   'BandWidth'         Pre-wavelet filter bandwidth in Hz. Default: 1
%   'Threshold'         Absolute phase-difference threshold in radians.
%                       Default: 1.2
%   'WaveletCycles'     Wavelet duration parameter passed to izmy_gbweeg.
%                       Default: 1
%   'UseParallel'       Enable PARFOR. Default: false
%   'Verbose'           Display progress messages. Default: true
%
%   ------
%   OUTPUT
%
%   results >
%       Structure containing:
%
%       .MSI
%           Metastability Index, defined as the temporal variance of the
%           Kuramoto order parameter. Dimensions:
%           [subjects x 1 x frequencies]
%
%       .nSCE
%           Normalized channel-wise Synchrony Coalition Entropy.
%           Dimensions:
%           [subjects x channels x frequencies]
%
%       .meanSCE
%           Mean normalized SCE across channels.
%           Dimensions:
%           [subjects x 1 x frequencies]
%
%       .patternValues
%           Observed synchrony-coalition pattern identifiers.
%
%       .patternCounts
%           Occurrence counts of the observed patterns.
%
%   --------------
%   SCE DEFINITION
%
%   For each reference channel, binary synchronization states are computed
%   against all other channels:
%
%       synchronized = abs(phase difference) < Threshold
%
%   The self-channel is excluded. For M total channels, each coalition
%   therefore consists of M-1 binary elements and has 2^(M-1) possible
%   states. The normalized channel-wise SCE is:
%
%       SCE_i = -sum_s p_i(s) log2(p_i(s)) / (M - 1)
%
%   Mean SCE is the arithmetic mean of SCE_i across the M channels.
%
%   ------------
%   DEPENDENCIES
%
%   - izmy_gbweeg      Included in this repository
%
%   -----
%   NOTES
%
%   - Input data must already be preprocessed and artifact-cleaned.
%   - CSD transformation is not performed by this function.
%   - Coalition patterns are encoded using uint64; therefore the current
%     implementation supports at most 64 partner channels.
%
%   ---------------------
%   AUTHORS AND COPYRIGHT
%
%   Copyright (c) 2026 Mebuki Izumiya
%
%   Maintainer:
%       Mebuki Izumiya
%       Division of Neural Dynamics
%       National Institute for Physiological Sciences
%
%   -------
%   LICENSE
%
%   SPDX-License-Identifier: BSD-3-Clause
%   - This file is part of EEG-Metastability and is distributed under the
%     BSD 3-Clause License. See the LICENSE file in the repository root.
%
%   See also izmy_gbweeg.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'inputPath', @(x) ischar(x) || isstring(x));

addParameter(ip, 'FilePattern', 'sub_*.mat', @(x) ischar(x) || isstring(x));
addParameter(ip, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'SaveResults', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'OutputFileName', 'msi_sce_results.mat', @(x) ischar(x) || isstring(x));

addParameter(ip, 'DataVariable', 'EEG', @(x) ischar(x) || isstring(x));
addParameter(ip, 'DataField', 'data', @(x) ischar(x) || isstring(x));

addParameter(ip, 'SampleRate', 1000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'Channels', 63, @(x) isnumeric(x) && isscalar(x) && x >= 2);
addParameter(ip, 'TimeIndices', [], @(x) isnumeric(x) && isvector(x));

addParameter(ip, 'FrequencyRange', [1 47], @(x) isnumeric(x) && numel(x) >= 2);
addParameter(ip, 'BandWidth', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'Threshold', 1.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'WaveletCycles', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'UseParallel', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(ip, inputPath, varargin{:});
opt = ip.Results;

inputPath = char(inputPath);
filePattern = char(opt.FilePattern);
outputDir = char(opt.OutputDir);
dataVariable = char(opt.DataVariable);
dataField = char(opt.DataField);
outputFileName = char(opt.OutputFileName);

% -------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------
if isfolder(inputPath)
    inputDir = inputPath;
    files = dir(fullfile(inputDir, filePattern));

elseif isfile(inputPath)
    [inputDir, fileName, fileExt] = fileparts(inputPath);
    files = dir(fullfile(inputDir, [fileName, fileExt]));

else
    error(['Input path does not exist. Please provide either a directory ', ...
           'containing EEG .mat files or a single EEG .mat file:\n%s'], inputPath);
end

if isempty(files)
    error('No files matching pattern "%s" were found in: %s', filePattern, inputDir);
end

if exist('izmy_gbweeg', 'file') ~= 2
    error('Required function "izmy_gbweeg.m" was not found on the MATLAB path.');
end

nSubjects = numel(files);
freqs = opt.FrequencyRange(:)';
nFreq = numel(freqs);
nCh = opt.Channels;

if opt.Verbose
    fprintf('Found %d file(s) in %s\n', nSubjects, inputDir);
    fprintf('Processing %d frequency bins: %s ... %s Hz\n', ...
        nFreq, num2str(freqs(1)), num2str(freqs(end)));
end

% -------------------------------------------------------------------------
% Preallocate outputs
% -------------------------------------------------------------------------
subjectIds    = nan(nSubjects, 1);
fileNames     = cell(nSubjects, 1);
MSI           = nan(nSubjects, 1, nFreq);
nSCE          = nan(nSubjects, nCh, nFreq);
meanSCE       = nan(nSubjects, 1, nFreq);

patternValues = cell(nSubjects, nFreq, nCh);
patternCounts = cell(nSubjects, nFreq, nCh);

% -------------------------------------------------------------------------
% Start parallel pool if requested
% -------------------------------------------------------------------------
if opt.UseParallel
    pool = gcp('nocreate');
    if isempty(pool)
        parpool;
    end
end

% -------------------------------------------------------------------------
% Subject loop
% -------------------------------------------------------------------------
if opt.UseParallel
    parfor s = 1:nSubjects
        [subjectIds(s), fileNames{s}, MSI_subj, nSCE_subj, ...
            meanSCE_subj, patternValues_subj, patternCounts_subj] = processOneFile( ...
            files(s), inputDir, dataVariable, dataField, opt, freqs, nCh);

        MSI(s,1,:) = MSI_subj;
        nSCE(s,:,:) = nSCE_subj;
        meanSCE(s,1,:) = meanSCE_subj;

        for f = 1:nFreq
            for ch = 1:nCh
                patternValues{s,f,ch} = patternValues_subj{f,ch};
                patternCounts{s,f,ch} = patternCounts_subj{f,ch};
            end
        end
    end
else
    for s = 1:nSubjects
        [subjectIds(s), fileNames{s}, MSI_subj, nSCE_subj, ...
            meanSCE_subj, patternValues_subj, patternCounts_subj] = processOneFile( ...
            files(s), inputDir, dataVariable, dataField, opt, freqs, nCh);

        MSI(s,1,:) = MSI_subj;
        nSCE(s,:,:) = nSCE_subj;
        meanSCE(s,1,:) = meanSCE_subj;

        for f = 1:nFreq
            for ch = 1:nCh
                patternValues{s,f,ch} = patternValues_subj{f,ch};
                patternCounts{s,f,ch} = patternCounts_subj{f,ch};
            end
        end
    end
end

% -------------------------------------------------------------------------
% Build results struct
% -------------------------------------------------------------------------
results = struct();
results.subjectIds = subjectIds;
results.fileNames = fileNames;
results.frequencies = freqs;
results.threshold = opt.Threshold;
results.sampleRate = opt.SampleRate;
results.channels = nCh;
results.timeIndices = opt.TimeIndices;

results.MSI = MSI;
results.nSCE = nSCE;
results.meanSCE = meanSCE;
results.patternValues = patternValues;
results.patternCounts = patternCounts;

results.metadata = struct();
results.metadata.inputPath = inputPath;
results.metadata.inputDir = inputDir;
results.metadata.filePattern = filePattern;
results.metadata.outputDir = outputDir;
results.metadata.dataVariable = dataVariable;
results.metadata.dataField = dataField;
results.metadata.bandWidth = opt.BandWidth;
results.metadata.waveletCycles = opt.WaveletCycles;
results.metadata.useParallel = logical(opt.UseParallel);
results.metadata.generatedAt = datestr(now, 30);
results.metadata.msiDefinition = 'Temporal variance of network-wide phase synchrony';
results.metadata.sceDefinition = 'Normalized Shannon entropy of channel-wise binary synchrony coalition patterns';

% -------------------------------------------------------------------------
% Optional save
% -------------------------------------------------------------------------
if opt.SaveResults
    if isempty(outputDir)
        error('OutputDir must be specified when SaveResults=true.');
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    save(fullfile(outputDir, outputFileName), 'results', '-v7.3');

    if opt.Verbose
        fprintf('Saved results to: %s\n', fullfile(outputDir, outputFileName));
    end
end

if opt.Verbose
    fprintf('Done.\n');
end

end

% =========================================================================
% Helper function
% =========================================================================
function [subjectId, fileName, MSI_subj, nSCE_subj, ...
          meanSCE_subj, patternValues_subj, patternCounts_subj] = processOneFile( ...
          fileInfo, inputDir, dataVariable, dataField, opt, freqs, nCh)

fileName = fileInfo.name;
filePath = fullfile(inputDir, fileName);

tmp = load(filePath);

if ~isfield(tmp, dataVariable)
    error('File "%s" does not contain variable "%s".', fileName, dataVariable);
end

S = tmp.(dataVariable);
if ~isstruct(S) || ~isfield(S, dataField)
    error('Variable "%s" in file "%s" does not contain field "%s".', ...
        dataVariable, fileName, dataField);
end

data = S.(dataField);

if ~isnumeric(data) || ndims(data) ~= 2
    error('Expected %s.%s in "%s" to be a numeric 2D matrix [channels x time].', ...
        dataVariable, dataField, fileName);
end

if size(data,1) < nCh
    error('File "%s" has only %d channels, but Channels=%d was requested.', ...
        fileName, size(data,1), nCh);
end

data = data(1:nCh, :);

if isempty(opt.TimeIndices)
    timeIdx = 1:size(data,2);
else
    timeIdx = opt.TimeIndices(:)';
    if max(timeIdx) > size(data,2)
        error('TimeIndices exceed data length in file "%s".', fileName);
    end
end

data = data(:, timeIdx);

tok = regexp(fileName, '\d+', 'match');
if isempty(tok)
    subjectId = NaN;
else
    subjectId = str2double(tok{1});
end

nFreq = numel(freqs);
MSI_subj = nan(1, nFreq);
nSCE_subj = nan(nCh, nFreq);
meanSCE_subj = nan(1, nFreq);

patternValues_subj = cell(nFreq, nCh);
patternCounts_subj = cell(nFreq, nCh);

for f = 1:nFreq
    cf = freqs(f);
    % complex wavelet transform
    cwtData = izmy_gbweeg(data, cf, opt.SampleRate, opt.WaveletCycles);

    if size(cwtData,1) ~= nCh
        error('izmy_gbweeg returned unexpected channel dimension for file "%s".', fileName);
    end

    % ---------------------------------------------------------------------
    % MSI: temporal variability of network-wide phase synchrony
    % ---------------------------------------------------------------------
    phaseData = angle(cwtData);                         % [channels x time]
    synchrony_t = abs(mean(exp(1i * phaseData), 1));    % [1 x time]
    MSI_subj(f) = var(synchrony_t, 0, 2, 'omitnan');

    % ---------------------------------------------------------------------
    % SCE: channel-wise synchrony coalition entropy
    % ---------------------------------------------------------------------
    numDataPoints = size(cwtData, 2);
    cwt3d = reshape(cwtData, [nCh, 1, numDataPoints]);
    cwt3dRep = repmat(cwt3d, [1, nCh, 1]);
    phaseDiff = angle(cwt3dRep ./ permute(cwt3dRep, [2, 1, 3]));

    binaryData = abs(phaseDiff) < opt.Threshold;
    scePerChannel = nan(nCh, 1);

    for ch = 1:nCh
        coalition = squeeze(binaryData(ch,:,:));   % [nCh x time]
        coalition(ch,:) = [];                      % remove self-channel
        nPartners = size(coalition, 1);

        if nPartners > 64
            error(['Current implementation uses uint64 pattern encoding and supports ', ...
                   'at most 64 coalition bits. Requested channel count is too large.']);
        end

        decdata = binaryRowsToUint64(coalition');  % [time x 1]
        [patVals, ~, idxPat] = unique(decdata, 'stable');
        patCounts = accumarray(idxPat, 1);

        P = patCounts ./ size(coalition, 2);
        sceBits = -sum(P .* log2(P), 'omitnan');

        % number of possible coalition states = 2^(nPartners)
        % therefore max entropy in bits = nPartners
        scePerChannel(ch) = sceBits / nPartners;   

        patternValues_subj{f, ch} = patVals;
        patternCounts_subj{f, ch} = patCounts;
    end

    nSCE_subj(:, f) = scePerChannel;
    meanSCE_subj(f) = mean(scePerChannel, 'omitnan');

    if opt.Verbose
        fprintf('[%s] %g Hz done\n', fileName, cf);
    end
end

end

% =========================================================================
% Convert binary rows to uint64
% rows: [nRows x nBits], values must be 0/1 or logical
% =========================================================================
function decdata = binaryRowsToUint64(rows)
rows = logical(rows);
[nRows, nBits] = size(rows);

if nBits > 64
    error('binaryRowsToUint64 supports at most 64 bits.');
end

powers = uint64(2) .^ uint64((nBits - 1):-1:0);
decdata = sum(uint64(rows) .* powers, 2, 'native');

if numel(decdata) ~= nRows
    error('Unexpected binary-to-decimal conversion error.');
end
end
