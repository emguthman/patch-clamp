function [hw,halfamp1,t_half1,halfamp2,t_half2] = findAPhalfwidth(data,sweep,sample_rate,threshtimes,peaktime,amp)
%defined as the width (in msec) of the AP waveform at the half-amplitude
halfamp1=zeros(length(threshtimes),1);
halfamp2=zeros(length(threshtimes),1);
hw=zeros(length(threshtimes),1);
t_half1=zeros(length(threshtimes),1);
t_half2=zeros(length(threshtimes),1);
for ii = 1:length(threshtimes)
    dfrom_halfamp1=[];
    dfrom_halfamp2=[]; %#ok<*NASGU>
    halfamp=amp(ii)/2;
    dfrom_halfamp1=abs(data(threshtimes(ii):peaktime(ii)-1,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist1]=min(dfrom_halfamp1);
    dfrom_halfamp2=abs(data(peaktime(ii):peaktime(ii)+0.0025*sample_rate,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist2]=min(dfrom_halfamp2);
    t_half1(ii)=threshtimes(ii)+t_mindist1-1;
    halfamp1(ii)=data(t_half1(ii),1,sweep);
    t_half2(ii)=peaktime(ii)+t_mindist2-1;
    halfamp2(ii)=data(t_half2(ii),1,sweep);
    hw(ii)=(t_half2(ii)-t_half1(ii))/(sample_rate/1000);
end
end