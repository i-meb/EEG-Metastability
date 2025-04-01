%% 250326

clear
tic

addpath('/fileserver/user/mebuki/eeglab2022.1')
addpath('/fileserver/user/mebuki/OSF/CSD_transform_for_suetani')
eeglab; close all hidden;

%% Define parameters
subj = 1; 
tt = 5001:15000;
ch = 63;
splrate = 1000;
freqnum = 47;

% Preallocate matrices to store results
nSCE_all = zeros(subj, ch+1, freqnum);
extractPat_all = zeros(subj, length(tt), freqnum, 'uint64');  
countPat_all   = zeros(subj, length(tt), freqnum, 'uint64');  

%% Load EEG data files
File = dir('sub_*.mat');
% File = File(~contains({File.name}, '_r.mat'));
%File = File(subj*4+1:subj*5, :);

% Start parallel pool if not already started
if isempty(gcp('nocreate'))
    parpool;
end

%% Loop through threshold values for phase synchronization
for th = 1.2
    % Preallocate local variables to gather results from each iteration
    local_nSCE_all = zeros(subj, ch, freqnum); 
    local_extractPat_all = zeros(subj, length(tt), freqnum, 'uint64'); 
    local_countPat_all   = zeros(subj, length(tt), freqnum, 'uint64'); 
    local_subjNo_all     = zeros(subj, 1); 

    parfor Filenum = 1:subj
        eegdata = load(File(Filenum).name, 'EEG');
        eegdata = eegdata.EEG.data(:,tt);
        subjNo = str2double(regexp(File(Filenum).name, '\d+', 'match'));

        local_nSCE = zeros(ch, freqnum);     
        local_extractPat = zeros(length(tt), freqnum);  
        local_countPat = zeros(length(tt), freqnum);    
        local_subjNo_all(Filenum, :) = subjNo;

        for freq = 1:freqnum 
            eegdataa = eegfilt(eegdata, splrate, freq, freq+1);

            % CSD transformation
            csdData = CSDconvert(eegdataa,'CSDbasis.mat');  
            eegdataa = []; 

            % Wavelet transform
            cwtData = izmy_gbweeg(csdData, freq, splrate, 1);
            [numChannels, numDataPoints] = size(cwtData);  
            csdData = [];

            %% SCE calculation
            cwt_3d = reshape(cwtData, [ch, 1, numDataPoints]);
            cwt_3drep = repmat(cwt_3d, [1, ch, 1]);  
            phaseData = angle(cwt_3drep ./ permute(cwt_3drep, [2, 1, 3]));  
            cwtData = []; cwt_3d = []; cwt_3drep =[];

            % Binary matrix based on phase difference threshold
            binaryData = zeros(ch, ch, numDataPoints);
            binaryData(abs(phaseData) <  th) = 1;  
            binaryData(abs(phaseData) >= th) = 0;  
            phaseData = [];

            % Temporary storage for this frequency
            extractPat_temp = zeros(length(tt), ch, 'uint64');
            countPat_temp   = zeros(length(tt), ch, 'uint64');

            for nch = 1:ch
                sq1 = squeeze(binaryData(nch,:,:));
                sq2 = [sq1(1:nch-1,:); sq1(nch+1:ch,:)]; 

                % Convert binary matrix to decimal for unique pattern identification
                binaryString = char(sq2' + '0');  % Convert sq2 to a char array of '0'/'1'

                % bitごとに2のべき乗を計算し、binaryStringに対応する2のべき乗行列を生成
                powersOfTwo = uint64(2) .^ uint64(61:-1:0);  
                
                % binaryStringを数値変換し、bit shift処理でuint64に一括変換
                decdata = sum(uint64(binaryString == '1') .* powersOfTwo, 2, 'native');

                % Extract unique patterns and count occurrences
                [extractPat, ~, patternIndices] = unique(decdata, 'stable');
                countPat = accumarray(patternIndices, 1);  %, [], @sum, int64(0)

                % Store patterns and counts temporarily
                extractPat_temp(1:length(extractPat),nch) = extractPat;
                countPat_temp(1:length(countPat),nch) = countPat;

                % Calculate SCE
                P = countPat ./ numel(tt);  
                SCE = sum(-P .* log2(P), 'omitnan');  
                nSCE = SCE / ch;  

                % Store result locally
                local_nSCE(nch, freq) = nSCE;
                disp(['th = ', num2str(th), '  subj', num2str(subjNo), '  ', num2str(freq), 'Hz', '  sce(', num2str(nch), ') = ', num2str(nSCE)]);
            end

            % Store temporary results for current frequency
            local_extractPat(1:size(extractPat_temp, 1),freq) = extractPat_temp(:,nch);
            local_countPat(1:size(countPat_temp, 1),freq) = countPat_temp(:,nch);
        end

        % Store local variables into global variables after loop
        local_nSCE_all(Filenum, :, :) = local_nSCE; 
        local_extractPat_all(Filenum, :, :) = local_extractPat;
        local_countPat_all(Filenum, :, :)   = local_countPat;
    end

    % Gather results into main matrices
    nSCE_all(:,1,:) = repmat(local_subjNo_all, [1, 1, freqnum]);
    nSCE_all(:,2:end,:) = local_nSCE_all;
    extractPat_all = local_extractPat_all;
    countPat_all = local_countPat_all;

    clear local_nSCE_all local_extractPat_all local_countPat_all

    %% Save results
    % thname = strrep(num2str(th), '.', ''); 
    % addpath(genpath('\\133.48.97.35\user\mebuki\OSF\Results\freqone\th12\120subj\sce\10sec'))
    % cd('\\133.48.97.35\user\mebuki\OSF\Results\freqone\th12\120subj\sce\10sec');
    % save(fullfile(pwd, ['sce_all', thname]), "extractPat_all", "countPat_all", "nSCE_all", '-v7.3', '-nocompression');
    %clear  nSCE_all extractPat_all countPat_all
end

toc
