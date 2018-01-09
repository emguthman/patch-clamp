function [peaks, times, number] = find_rspks(data,sweep,t1,stop)
%gets rebound spikes, define as voltage peaks above -10mV

[allpeaks,alltimes]=findpeaks(data(t1:stop,1,sweep));

peaks=allpeaks(allpeaks>-15);
times=alltimes(allpeaks>-15)-1+t1;
number=length(peaks);
end