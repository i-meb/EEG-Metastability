function h=izmy_gbweeg(data,freq,sample_rate,nco)

[channel, time] = size(data); 

    sigma = nco/(6*freq);
    windowlength = fix(((nco/freq) * sample_rate) / 2);

    me = squeeze(mean(data, 2));
    %data(:,:,1) = data(:,:,1)-me(:,1);
    data = data - repmat(me, 1, time, 1);
    
    time = -(windowlength/sample_rate):(1/sample_rate):(windowlength/sample_rate);
    mother=sqrt(freq)*exp(-(time.^2)/(2*(sigma^2))).*exp(1i*2*pi*freq*time);
    validLength = 5000 - 2 * windowlength;  % align tt(in Calc code)
    h = zeros(channel, validLength);

    for channels = 1:channel
        temp=conv(data(channels,:),mother,'valid');
        h(channels,:) = temp;
    end
end
%     convdata=data;
%     
%     ch=1;
%     trial=1;
% 
%     temp=conv(data(ch,:,trial),mother);
%     convdata(ch,:,trial)=temp(windowlength+1:(end-windowlength));
%     h(i,:)=convdata;
% end
