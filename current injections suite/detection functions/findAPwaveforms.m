function [waveforms] = findAPwaveforms(data,peak_times,sweep,sample_rate)
%gets trace for each AP waveform for plotting purposes

waveforms=zeros(length(peak_times),41);
for ii = 1:length(peak_times)
    waveforms(ii,:)=data(peak_times(ii)-0.0015*sample_rate:peak_times(ii)+0.0025*sample_rate,1,sweep);
end
end