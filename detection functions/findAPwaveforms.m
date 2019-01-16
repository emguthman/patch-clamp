function [waveforms] = findAPwaveforms(data,centerTimes,sweep,sample_rate)
%gets trace for each AP waveform for plotting purposes

waveforms=zeros(length(centerTimes),66);
for ii = 1:length(centerTimes)
    waveforms(ii,:)=data(centerTimes(ii)-0.0015*sample_rate:centerTimes(ii)+0.005*sample_rate,1,sweep);
end
end
