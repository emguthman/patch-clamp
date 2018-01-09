%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Action Potential Threshold %%
%%%%% for current injection protocol %%%%%
%%%%%%%%%% Created: 11-23-2015 %%%%%%%%%%%
%%%%%%%%%%% Edited: 07-08-2017 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thresh_vals,thresh_times]=find_apthresh(CS_data,blstart,blend,injstart,injend,samplerate,dI)
no_sweeps=size(CS_data,3); %gives number of sweeps per file

%init vars
thresh_vals=[];
thresh_times=[];
thresh = 50; %can also set to be a function of the sd dVm
dVm = diff(CS_data); %takes approx derivative (lose element n = 1)
dVm(:,2,:) = [];
dVm = dVm.*(samplerate/1000); %converts to V/s
CS_data(1,:,:) = []; %remove element n = 1 from Vm data
t = samplerate^-1:samplerate^-1:size(CS_data,1)/samplerate;
sddVm = std(dVm(blstart:blend)); %get standard deviation of the derivative of baseline V-mem
dVm_thresh = 20*ones(1,no_sweeps); %manually set to 20V/s
apthresh = -1.*ones(250,2,no_sweeps);

% Set up CS figure
for ii = 1:no_sweeps
    %find current step
    t_thresh=[];
    over_thresh=[];
    if dI(ii) > 0
        %find spikes
        noSpkSweep=1;
        over_thresh=find(dVm((.002*samplerate+injstart):injend,1,ii)>dVm_thresh(ii)); %find points above the threshold of dmV during current step (2ms buffer from start of step)
        for jj = 1:length(over_thresh)
            if jj == 1
                %get point before and after threshold
                p0 = over_thresh(jj)+injstart+.002*samplerate-1; %after
                p1 = over_thresh(jj)+injstart+.002*samplerate-2; %before
                afterVal = abs(20-dVm(p0,1,ii));
                beforeVal = abs(20-dVm(p1,1,ii));
                if afterVal <= beforeVal %if value after threshold is closer to threshold...
                    t_thresh(noSpkSweep) = p0;
                else %if value before threshold is closer to threshold...
                    t_thresh(noSpkSweep) = p1;
                end
                noSpkSweep=noSpkSweep+1;
            else
                if over_thresh(jj)>.001*samplerate+t_thresh(noSpkSweep-1)-injstart %if time over threshold is at least 1ms past the last AP (refractory period)
                    p0 = over_thresh(jj)+injstart+.002*samplerate-1; %after
                    p1 = over_thresh(jj)+injstart+.002*samplerate-2; %before
                    afterVal = abs(20-dVm(p0,1,ii));
                    beforeVal = abs(20-dVm(p1,1,ii));
                    if afterVal <= beforeVal %if value after threshold is closer to threshold...
                        t_thresh(noSpkSweep) = p0;
                    else %if value before threshold is closer to threshold...
                        t_thresh(noSpkSweep) = p1;
                    end
                    noSpkSweep=noSpkSweep+1;
                end
            end
        end
        noSpkSweep = noSpkSweep -1
        
        %Plots AP traces and thresholding
        figure(1)
        subplot(1,2,1)
        hold on
        plot(t,CS_data(:,1,ii))
        ylim([min(min(CS_data(:,1,:)))-0.1*max(range(CS_data(:,1,:))) max(max(CS_data(:,1,:)))+0.1*max(range(CS_data(:,1,:)))])
        ylabel('Membrane Potential (mV)','fontweight','bold')
        xlabel('Time (ms)','fontweight','bold')
        xlim([.5 1.75])
        if isempty(t_thresh) ~= 1
            scatter(t_thresh/samplerate,CS_data(t_thresh,1,ii),'markeredgecolor',[0 .4 0])
        end
        subplot(1,2,2)
        hold on
        plot(t,dVm(:,1,ii))
        line([t(1) t(end)],[dVm_thresh(ii) dVm_thresh(ii)],'color','k','linewidth',1.25);
        ylabel('dV/dt (V/sec)')
        xlabel('Time (ms)')
        xlim([.5 1.75])
        clf
    end
    
    %Fill Matrix
    if isempty(t_thresh)~=1
        if isempty(thresh_vals) == 1 && isempty(thresh_times) ==1
            thresh_vals(ii,:)=CS_data(t_thresh,1,ii);
            thresh_times(ii,:)=t_thresh;
        elseif length(t_thresh)==size(thresh_vals,2)
            thresh_vals(ii,:)=CS_data(t_thresh,1,ii);
            thresh_times(ii,:)=t_thresh;
        elseif length(t_thresh)<size(thresh_vals,2)
            no_rows = size(thresh_vals,2)-length(t_thresh);
            addrows=zeros(1,no_rows);
            thresh_vals(ii,:)=vertcat(CS_data(t_thresh,1,ii),addrows');
            t_thresh=horzcat(t_thresh,addrows);
            thresh_times(ii,:)=t_thresh;
        elseif length(t_thresh)>size(thresh_vals,2)
            no_col = length(t_thresh)-size(thresh_vals,2);
            addcol=zeros(size(thresh_vals,1),no_col);
            thresh_vals=horzcat(thresh_vals,addcol);
            thresh_times=horzcat(thresh_times,addcol);
            thresh_vals(ii,:)=CS_data(t_thresh,1,ii);
            thresh_times(ii,:)=t_thresh;
        end
    elseif isempty(t_thresh)==1
        thresh_vals(ii,:)=0;
        thresh_times(ii,:)=0;
    end
end
