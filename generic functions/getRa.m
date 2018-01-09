function [Ra, meanRa, keep] = getRa(pClamp_data,injstart,injend,sample_rate)
% Get Access Resistane Function, created 7-26-2016, modified 08-15-2016
% takes sweeps and outputs access resistance
% also returns the mean Ra
% keep is a logical variable, if max(Ra) <= 25 and does not vary > +/-20% from
% meanRa, keep == 1; if Ra violates these conditions, keep == 0
% assumes that dV is 10mV -- since this sometimes is recorded incorrectly
% all my protocols use a 10mV step for checking Ra. (EMG 08-15-2016)

%init vars
if length(size(pClamp_data)) == 3
    no_sweeps=size(pClamp_data,3);
    no_files=1;
elseif length(size(pClamp_data)) == 4
    no_sweeps=size(pClamp_data,3);
    no_files=size(pClamp_data,4);
elseif length(size(pClamp_data)) < 3
    error('Not enough input arguments in pClamp Data file -- it may be filtered')
end

if no_files == 1 %calc Ra if only one .abf file
    dV=10.*ones(no_sweeps,1);
    dI=zeros(no_sweeps,1);
    Ra=zeros(no_sweeps,1);
    for ii = 1:no_sweeps
        %dV(ii) = (mean(pClamp_data((injstart-.51*sample_rate):(injstart-.01*sample_rate),2,ii)) - mean(pClamp_data(injstart:injend,2,ii))); %difference bw holding voltage and voltage step (mV)
        dI(ii) = mean(pClamp_data((injstart-.051*sample_rate):(injstart-.01*sample_rate),1,ii)) - min(pClamp_data(injstart:injend,1,ii)); %difference bw baseline I and Vstep evoked peak I (pA)
        Ra(ii) = 1000*(dV(ii)/dI(ii)); %in MOhms
    end    
elseif no_files > 1 %calc Ra if multiple .abf files
    dV=10.*ones(no_sweeps*no_files,1);
    dI=zeros(no_sweeps*no_files,1);
    Ra=zeros(no_sweeps*no_files,1);
    for jj = 1:no_files
        for ii = 1:no_sweeps
            %dV(ii+jj*no_sweeps-no_sweeps) = (mean(pClamp_data((injstart-.51*sample_rate):(injstart-.01*sample_rate),2,ii,jj)) - mean(pClamp_data(injstart:injend,2,ii,jj))); %difference bw holding voltage and voltage step (mV)
            dI(ii+jj*no_sweeps-no_sweeps) = mean(pClamp_data((injstart-.51*sample_rate):(injstart-.01*sample_rate),1,ii,jj)) - min(pClamp_data(injstart:injend,1,ii,jj)); %difference bw baseline I and Vstep evoked peak I (pA)
            Ra(ii+jj*no_sweeps-no_sweeps) = 1000*(dV(ii+jj*no_sweeps-no_sweeps)/dI(ii+jj*no_sweeps-no_sweeps)); %in MOhms
        end
    end
end

meanRa=mean(Ra);

%exclude if it doesn't meet the above criteria
if max(Ra) > 25
    keep = false;
elseif (max(Ra) > meanRa * 1.2) || (min(Ra) < meanRa * .8)
    keep = false;
else
    keep = true;
end

%plot Ra to view
Rafig=figure;
plot(Ra,'-o')
line([0 no_sweeps*no_files+1],[meanRa * 1.2, meanRa * 1.2],'LineStyle','--','Color','k')
line([0 no_sweeps*no_files+1],[meanRa * .8, meanRa * .8],'LineStyle','--','Color','k')
ylim([min(Ra)-5 max(Ra)+5])
xlim([0 no_sweeps*no_files+1])
close(Rafig)

end