%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculates Membrane Resistance %%%%
%%%%%%%%%% Created: 09-10-2015 %%%%%%%%%%
%%%%%%%%%%% Edited: 01-28-2016 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rm, dp]=MemRes(CS_data,blstart,blend,injend,sample_rate,threshold)
%%% takes CS Step data imported from pClamp with abfload.m
%all sweeps prior to rheobase method, last 200ms of injection
if threshold > 3
    
    limit = 1; %must have more data points than this to obtain Rm
    
    %IV plot to get membrane resistance
    %get dI and dV
    dI=zeros(1,threshold-1);
    dV=zeros(1,threshold-1);
    for ii = 1:threshold-1
        dI(ii) = roundn(mean(CS_data(injend-0.2*sample_rate:injend,2,ii))-mean(CS_data(blstart:blend,2,ii)),1); %difference between current step and holding current
        dV(ii) = mean(CS_data(injend-0.2*sample_rate:injend,1,ii))-mean(CS_data(blstart:blend,1,ii)); %difference between evoked voltage and baseline voltage -- last 200ms of 600ms current injection
    end
    
    %get no data points
    dp = length(dV); %number of data points in V-I approximation
    
    %get linear fit of V-I curve
    VIfit=nlinfit(dI,dV,@dr_fitline,[0 1]);
    plot_VIfit=VIfit(1).*dI+VIfit(2);
    if dp <=limit
        Rm=NaN;
    else
        Rm=1000*VIfit(1); %gives you Rm in MOhms
    end
    
    %Plot V-I
    figure(99)
    hold on
    scatter(dI,dV,50,[1 .5 .5],'fill')
    plot(dI,plot_VIfit)
    ylabel('Voltage (mV)','fontweight','bold')
    xlabel('Current (pA)','fontweight','bold')
    title('V-I plot, all sweeps prior to rheobase method','fontweight','bold')
    close(figure(99))
elseif threshold == 1
    Rm=NaN
    dp=NaN
end