function vsag = findvsag(data,sweep,injstart,injend,blstart,blend,sample_rate)
%returns voltage sag
% 100*(minV-steadystateV)/(minV-baselineV)  
% minV is the mean value of the minumum point in the hyperpolarizing injection
% steadystateV is the mean value over the last 200ms of the hyperpolarizing injection
% baselineV is the mean value over the first 500ms of the sweep

minV=min(data(injstart:injend,1,sweep));
ssV=mean(data((injend-0.2*sample_rate):injend,1,sweep));
blV=mean(data(blstart:blend,1,sweep));

vsag=100*((minV-ssV)/(minV-blV));

end