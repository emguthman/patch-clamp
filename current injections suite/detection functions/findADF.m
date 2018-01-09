function [ADF] = findADF(data,blstart,blend,injend,sweep,sample_rate)
%find after deflection as calculated in Cea-del Rio et al 2010
ADF=mean(data(injend+.2*sample_rate:injend+.3*sample_rate,1,sweep))-mean(data(blstart:blend,1,sweep));

end