function [FR] = findFR(thresh_times,sweep,window_start,window_end,sample_rate)
%finds FR within a specified window

spks=(thresh_times(sweep,:)>=window_start)&(thresh_times(sweep,:)<=window_end);
tspks=thresh_times(sweep,spks); %gets time for spk threshold

%get ISI in sec
for ii = 1:length(tspks)
    if ii ~=length(tspks)
        ISI(ii)=(tspks(ii+1)-tspks(ii))/sample_rate;
    else
        if ii ~= length(thresh_times(sweep,:))
            ISI(ii)=(thresh_times(sweep,ii+1)-tspks(ii))/sample_rate;
        end
    end
end
meanISI=mean(ISI);

%get FR in Hz
FR=1/meanISI;

end