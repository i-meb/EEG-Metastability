%% 

function h = izmy_gbweeg(data, freq, sample_rate, nco)

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
