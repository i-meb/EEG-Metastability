function h = izmy_gbweeg(data, freq, sample_rate, nco)

% izmy_gbweeg:  Compute complex Gabor-wavelet coefficients for EEG signals.
%
%   h = izmy_gbweeg(data, freq, sample_rate, nco) convolves each EEG
%   channel with a complex Gaussian-windowed sinusoid centered at 'freq'.
%
%   INPUTS
%   ------
%   data
%       EEG signal matrix with dimensions:
%
%           [channels x time samples]
%
%   freq
%       Wavelet center frequency in Hz. Positive scalar.
%
%   sampleRate
%       Sampling frequency in Hz. Positive scalar.
%
%   nCycles
%       Wavelet-duration parameter. The Gaussian standard deviation is:
%
%           sigma = nCycles / (6 * freq)
%
%       and the half-window duration is:
%
%           nCycles / (2 * freq)
%
%   OUTPUT
%   ------
%   h
%       Complex wavelet coefficients with dimensions:
%
%           [channels x valid time samples]
%
%       Because convolution is performed using the 'valid' option, the
%       output length is:
%
%           size(data,2) - length(motherWavelet) + 1
%
%   METHOD
%   ------
%   The signal mean is removed separately from each channel. The complex
%   mother wavelet is then defined as:
%
%       sqrt(freq) ...
%       .* exp(-t.^2 / (2*sigma^2)) ...
%       .* exp(1i*2*pi*freq*t)
%
%   Each EEG channel is convolved with this wavelet using CONV(...,'valid').
%
%   PROVENANCE
%   ----------
%   Based on an earlier Gabor-wavelet implementation developed by
%   Keiichi Kitajo. Adapted, documented, and integrated into the
%   EEG-Metastability pipeline by Mebuki Izumiya.
%
%   AUTHORS AND COPYRIGHT
%   ---------------------
%   Copyright (c) 2026 Keiichi Kitajo
%
%   Contact:
%       Keiichi Kitajo
%       Mebuki Izumiya
%
%   LICENSE
%   -------
%   SPDX-License-Identifier: BSD-3-Clause
%
%   This file is part of EEG-Metastability and is distributed under the
%   BSD 3-Clause License. See the LICENSE file in the repository root.
%
%   See also calcMSISCE, conv.

[channel, nTime] = size(data);

sigma = nco / (6 * freq);
windowlength = fix(((nco / freq) * sample_rate) / 2);

% demean across time for each channel
me = mean(data, 2, 'omitnan');
data = data - repmat(me, 1, nTime);

% wavelet time vector
t = -(windowlength / sample_rate):(1 / sample_rate):(windowlength / sample_rate);

mother = sqrt(freq) ...
    .* exp(-(t .^ 2) ./ (2 * sigma ^ 2)) ...
    .* exp(1i * 2 * pi * freq * t);

% conv(...,'valid') length = nTime - length(mother) + 1
validLength = nTime - length(mother) + 1;

if validLength <= 0
    error(['Wavelet window is longer than the input data. ', ...
           'Increase data length or reduce WaveletCycles.']);
end

h = zeros(channel, validLength);

for ch = 1:channel
    temp = conv(data(ch, :), mother, 'valid');
    h(ch, :) = temp;
end

end
