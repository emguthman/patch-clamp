function [freq] = instFreq(times)
%calculates instantaneous frequency of events
%assumes times are in seconds
%created 08-25-2016; modified 08-25-2016

freq=zeros(length(times)-1,1);
for ii = 1:length(times)-1 %does not get instFreq for last event (no event to compare it to)
    freq(ii) = 1/(times(ii+1)-times(ii));
end
end