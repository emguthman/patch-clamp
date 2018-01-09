function [peakval, tpeak, nopeaks] = findspkpeaks(data,thresh_vals,thresh_times,sweep)
%detects action potential peaks

nopeaks=sum(thresh_vals(sweep,:)~=0);
peakval=zeros(nopeaks,1);
tpeak=zeros(nopeaks,1);
if nopeaks~=0 %if APs in current sweep
    for ii=1:nopeaks
        if ii ~= nopeaks
            [peakval(ii),tpeak(ii)]=max(data(thresh_times(sweep,ii):thresh_times(sweep,ii+1),1,sweep)); %finds max point between threshold_n and threshold_n+1
        elseif ii == nopeaks
            [peakval(ii),tpeak(ii)]=max(data(thresh_times(sweep,ii):thresh_times(sweep,ii)+100,1,sweep)); %finds max point between threshold_n and threshold_n+1
        end
        tpeak(ii)=tpeak(ii)+thresh_times(sweep,ii)-1; %fixes time
    end
end
end