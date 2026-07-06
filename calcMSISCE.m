function results = calcMSISCE(inputPath, varargin)
% calcMSISCE
%
% Calculate Metastability Index (MSI) and channel-wise
% Synchrony Coalition Entropy (SCE) from EEG data.
%
% In this implementation, MSI is defined as the temporal variance
% of network-wide phase synchrony across channels at each frequency.
%
% SCE is computed for each channel by thresholding phase differences
% between the reference channel and all other channels, constructing
% binary synchrony coalition patterns across time, and evaluating the
% Shannon entropy of the pattern distribution.
%
% This function is intended for open-source use and avoids hard-coded paths.
% It processes MATLAB .mat files containing an EEG structure with data in
% EEG.data (channels x timepoints).
%
% REQUIRED EXTERNAL FUNCTIONS / TOOLBOXES
%   - eegfilt      (EEGLAB)
%   - izmy_gbweeg  (included in this repository)
%
% INPUT
%   inputPath : directory containing EEG .mat files, or a single EEG .mat file
%
% PARAMETERS
%   'FilePattern'        : filename pattern (default: 'sub_*.mat')
%   'OutputDir'          : directory to save output (default: '')
%   'SaveResults'        : true/false (default: false)
%   'OutputFileName'     : .mat filename if saving (default: 'msi_sce_results.mat')
%   'DataVariable'       : top-level variable name in .mat (default: 'EEG')
%   'DataField'          : field inside structure containing signal (default: 'data')
%   'SampleRate'         : sampling rate in Hz (default: 1000)
%   'Channels'           : number of channels to use (default: 63)
%   'TimeIndices'        : indices to analyze (default: [])
%   'FrequencyRange'     : vector of center frequencies (default: 1:47)
%   'BandWidth'          : band-pass width in Hz for eegfilt (default: 1)
%   'Threshold'          : phase-difference threshold in radians (default: 1.2)
%   'WaveletCycles'      : nco parameter for izmy_gbweeg (default: 1)
%   'UseParallel'        : true/false (default: true)
%   'Verbose'            : true/false (default: true)
%
% OUTPUT
%   results : struct with fields
%       .subjectIds            [nSubjects x 1]
%       .fileNames             {nSubjects x 1}
%       .frequencies           [1 x nFreq]
%       .threshold             scalar
%       .sampleRate            scalar
%       .channels              scalar
%       .timeIndices           vector
%       .MSI                   [nSubjects x 1 x nFreq]
%       .nSCE                  [nSubjects x nChannels x nFreq]
%       .meanSCE               [nSubjects x 1 x nFreq]
%       .patternValues         {nSubjects, nFreq, nChannels}
%       .patternCounts         {nSubjects, nFreq, nChannels}
%       .metadata              struct
%
% NOTE
%   - This function assumes each input file contains one continuous recording
%     in EEG.data with shape [channels x timepoints].
%   - CSD preprocessing is NOT included in this function.
%     If required, it must be applied externally.
%
% EXAMPLE
%   results = calcMSISCE('./data', ...
%       'OutputDir', './results', ...
%       'SaveResults', true);

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

addParameter(ip, 'FrequencyRange', 1:47, @(x) isnumeric(x) && isvector(x) && all(x > 0));
addParameter(ip, 'BandWidth', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'Threshold', 1.2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'WaveletCycles', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'UseParallel', true, @(x) islogical(x) || isnumeric(x));
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

if exist('eegfilt', 'file') ~= 2
    error(['Required function "eegfilt" was not found. ', ...
           'Please add EEGLAB to your MATLAB path.']);
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
    lowCut = cf;
    highCut = cf + opt.BandWidth;

    % band-pass
    filteredData = eegfilt(data, opt.SampleRate, lowCut, highCut);

    % complex wavelet transform
    cwtData = izmy_gbweeg(filteredData, cf, opt.SampleRate, opt.WaveletCycles);

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
        % P = patCounts ./ numel(opt.TimeIndices);  % 削除してよい
        sceBits = -sum(P .* log2(P), 'omitnan');

        % number of possible coalition states = 2^(nPartners)
        % therefore max entropy in bits = nPartners
        maxEntropy = nPartners;
        scePerChannel(ch) = sceBits / nCh;

        patternValues_subj{f, ch} = patVals;
        patternCounts_subj{f, ch} = patCounts;
    end

    nSCE_subj(:, f) = scePerChannel;
    meanSCE_subj(f) = mean(scePerChannel, 'omitnan');

    if opt.Verbose
        fprintf('[%s] freq=%g Hz done\n', fileName, cf);
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
