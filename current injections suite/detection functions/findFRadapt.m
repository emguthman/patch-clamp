function [ratio] = findFRadapt(threshtimes,sweep,nospks)
%ratio of last ISI to first ISI in sweep

tspks=threshtimes(sweep,1:nospks-1); %gets time for spk threshold
%get ISI in dps
for ii = 1:length(tspks)
    if ii ~=length(tspks)
        ISI(ii)=(tspks(ii+1)-tspks(ii));
    else
        ISI(ii)=(threshtimes(sweep,ii+1)-tspks(ii));
    end
end

if length(ISI) >= 3
    ratio=ISI(1)/mean([ISI(end) ISI(end-1) ISI(end-2)]);
else
    ratio=NaN;
end
end