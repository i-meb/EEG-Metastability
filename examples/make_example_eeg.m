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
