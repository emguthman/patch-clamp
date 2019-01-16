%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% sEvents for sEvent suite %%%%%%%%%%
%%%%%%%%%%% Created: 03-23-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 10-18-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decrease buffer for EPSCs and change to .33 cuttoff

function [outputdata] = sponPSCFxn(inputs)

%% INIT DATA
close all

%sample rate
samplerate = inputs.SamplerateHzEditField.Value;

%% LOAD DATA
iData = -1.*ones(inputs.LengthofsweepsecEditField.Value*samplerate,inputs.NumberofsweepsEditField.Value);
rawData = -1.*ones(inputs.LengthofsweepsecEditField.Value*samplerate,2,inputs.NumberofsweepsEditField.Value);

%load
cd(inputs.spscDir);
contents = dir('*.abf');
fn = {contents.name}';
spscFile = fullfile(cd,fn);
pCFile = abfload(spscFile{1},'sweeps','a');
iData(:,:) = sgolayfilt(pCFile(:,1,:),3,11,[],1); %Savitsky-Golay 3rd Order Filter w +/- 0.5ms window on recorded current
rawData(:,:,:) = pCFile;

%% CHECK ACCESS
%injection for Ra
injstart = inputs.RainjectionstartmsEditField.Value*(samplerate/1000); %voltage injection start
injend = inputs.RainjectionendmsEditField.Value*(samplerate/1000); %voltage injection end

%Calculate access resistance
[Ra, meanRa, keep] = getRa(rawData,injstart,injend,samplerate) %#ok<*NOPRT>


%% ANALYSIS
if keep == 1
    
    %remove V-steps for Ra, with before and after for buffer
    iDataNoVstep = iData;
    iDataNoVstep((injstart-.001*samplerate):(injend+.099*samplerate),:)=[]; %removes 1.1sec of data from 1ms before Vstep to 99ms after
    
    %collate sweeps
    outputdata.sponData = reshape(iDataNoVstep,[],1); %reshapes the matrix of trials into a single vector
    
    %time vector
    t = samplerate^-1:samplerate^-1:length(outputdata.sponData)/samplerate;
    
    %find spontaneous events
    noEvents = 0;
    bufferEnd = 0;
    
    %get current template
    iTemplate = findTemplate(outputdata.sponData,samplerate,inputs);
    
    %truncate template to rise time and tau decay
    if inputs.CurrentSwitch.Value == 'EPSC'
        %20-80 rt
        [peakVal, peakInd] = min(iTemplate);
        [rt, ~, tTwenty] = findRise(iTemplate,1,peakInd,samplerate);
        
        %decay tau
        tau = findTau(iTemplate,1,samplerate,'Single');
        
    elseif inputs.CurrentSwitch.Value == 'IPSC'
        %20-80 rt
        [peakVal, peakInd] = max(iTemplate);
        [rt, ~, tTwenty] = findRise(iTemplate,1,peakInd,samplerate);
        
        %decay tau
        tau = findTau(iTemplate,1,samplerate,'Double');
        
    end
    
    %truncate
    iTempTruncated = iTemplate(tTwenty:peakInd+ceil(samplerate*tau/1000));
    preTime = peakInd-tTwenty;
    postTime = ceil(samplerate*tau/1000);
    
    %get charge integral
    chargeIntegral = trapz(iTempTruncated).*(1000/samplerate);
    cumCI = cumtrapz(iTempTruncated).*(1000/samplerate);
    
    %plot on template figure
    figure(1)
    plot(tTwenty:length(iTempTruncated)+tTwenty-1,iTempTruncated,'linewidth',5,'color','k');
    
    
    %get all potential spontaneous events
    [eventTimes, eventBaseline] = findEvents(outputdata.sponData,samplerate,inputs);
    
    %plot
    rawTraceFig = figure(2);
    rawTraceFig.Position = [300 15 750 355];
    hold on
    plot(t,outputdata.sponData,'color',inputs.plotColor(7,:))
    scatter(t(eventTimes),outputdata.sponData(eventTimes),'markeredgecolor',inputs.plotColor(3,:)) %all possible events, pre template search
    xlim([t(1) t(end)])
    xlabel('Time (s)')
    ylabel('current (pA)')
    title('raw trace')
    thisAx = gca;
    setAx(thisAx);
    
    %loop thru events to find those that match template
    %     normRT = -1.*ones(length(eventTimes),1);
    %     normTau = -1.*ones(length(eventTimes),1);
    normCI = -1.*ones(length(eventTimes),1);
    for ii = 1:length(eventTimes)
        %get event
        %         thisEventTrace =
        %         outputdata.sponData(eventTimes(ii)-2*preTime:eventTimes(ii)+3*postTime)-eventBaseline(ii); %for kinetic fitting
        thisEventTrace = outputdata.sponData(eventTimes(ii)-preTime:eventTimes(ii)+postTime)-eventBaseline(ii); %for integral fitting
        
        %fit by kinetics -- not all individual traces can be fit, makes sense as not all are synaptic currents
        %         if inputs.CurrentSwitch.Value == 'EPSC'
        %             %20-80 rt
        %             [thisEventRT, ~, ~] = findRise(thisEventTrace,1,1+0.01*samplerate,samplerate);
        %
        %             %decay tau
        %             thisEventTau = findTau(thisEventTrace,1+2*preTime,samplerate,'Single');
        %
        %         elseif inputs.CurrentSwitch.Value == 'IPSC'
        %             %20-80 rt
        %             [thisEventRT, ~, ~] = findRise(thisEventTrace,1,1+2*preTime,samplerate);
        %
        %             %decay tau
        %             thisEventTau = findTau(thisEventTrace,1+2*preTime,samplerate,'Double');
        %
        %         end
        %
        %         %normalize to template
        %         normRT(ii) = thisEventRT./rt; thisEventRT = [];
        %         normTau(ii) = thisEventTau./tau; thisEvenTau = [];
        
        %scale for template
        thisScale = [];
        thisScale = (outputdata.sponData(eventTimes(ii))-eventBaseline(ii))./peakVal;
        
        %fit by charge integral
        thisEventChargeIntegral = trapz(thisEventTrace).*(1000/samplerate); %in ms
        thisTemplateCharge = trapz(thisScale.*iTempTruncated).*(1000/samplerate); %in ms
        cumChargeDiff = 100.*((cumtrapz(thisScale.*iTempTruncated).*(1000/samplerate)) - (cumtrapz(thisEventTrace).*(1000/samplerate)))./(cumtrapz(thisScale.*iTempTruncated).*(1000/samplerate)); %normalized to charge integral, in percent
        normCI(ii) = thisEventChargeIntegral./thisTemplateCharge; thisEventChargeIntegral = []; thisTemplateCharge = [];
        chargeComputationProgress = 100*(ii/length(eventTimes))
        
        %if start of putative event is of greater amplitude than peak, "event" is actually just a blip on downslope of another event
        if inputs.CurrentSwitch.Value == 'EPSC'
            if mean(thisEventTrace(1:.001*samplerate)) <= (outputdata.sponData(eventTimes(ii))-eventBaseline(ii))
                normCI(ii) = -1;
            end
        elseif inputs.CurrentSwitch.Value == 'IPSC'
            if mean(thisEventTrace(1:.001*samplerate)) >= (outputdata.sponData(eventTimes(ii))-eventBaseline(ii))
                normCI(ii) = -1;
            end
        end
        
        %plot to check fit
        fitFig = figure(3);
        fitFig.Position = [475 200 500 500];
        subplot(4,2,1:6)
        hold on
        %         plot(outputdata.sponData(eventTimes(ii)-preTime:eventTimes(ii)+postTime)-eventBaseline(ii),'color',inputs.plotColor(1,:),'linewidth',3) % for kinetics data
        plot(thisEventTrace,'color',inputs.plotColor(3,:),'linewidth',3)
        plot(iTempTruncated.*thisScale,'color','k','linestyle','--','linewidth',2)
        setAx(gca);
        ylabel('current (pA)')
        xlabel('time (data points)')
        title('rise and decay of template and current event')
        subplot(4,2,7:8)
        plot(cumChargeDiff,'color',inputs.plotColor(7,:),'linewidth',3)
        line([0 length(cumChargeDiff)],[0 0],'linewidth',1,'linestyle','--','color','k')
        setAx(gca);
        xlabel('time (data points)')
        ylabel('\Delta Charge Integral (%)')
        title('Charge Integral, cumulative change from template')
        clf
    end
    
    %look at histogram
    histFig = figure(3);
    histFig.Position = [1060 600 350 180];
    nbins = ceil((max(normCI)-min(normCI))/.25); %bins at 25% for zoomed out
    subplot(1,2,1)
    histogram(normCI,nbins,'EdgeColor',[1 1 1],'FaceColor',inputs.plotColor(7,:))
    wholeHistAx = gca;
    setAx(wholeHistAx);
    ylabel('count')
    xlabel('normalized charge integral')
    subplot(1,2,2)
    nbins = ceil((max(normCI)-min(normCI))/.025); %bins at 2.5% for zoomed out
    histogram(normCI,nbins,'EdgeColor',[1 1 1],'FaceColor',inputs.plotColor(7,:))
    setAx(gca);
    ylabel('count')
    xlabel('normalized charge integral')
    xlim([0 1])
    
    % get putative events
    ciCutoff = inputdlg('charge integral cut-off (in normalized units):');
    ciCutoff = str2double(ciCutoff);
    outputdata.putativeEvents = find(normCI >= ciCutoff);
    outputdata.sponEventTraces = -1.*ones(.1*samplerate+1,length(outputdata.putativeEvents));
    for ii = 1:length(outputdata.putativeEvents)
        gettingEvents  = 100*(ii./length(outputdata.putativeEvents))
        outputdata.sponEventTraces(:,ii) = outputdata.sponData(eventTimes(outputdata.putativeEvents(ii))-.025*samplerate:eventTimes(outputdata.putativeEvents(ii))+.075*samplerate)-eventBaseline(outputdata.putativeEvents(ii));
    end
    
    %align and get mean
    if inputs.CurrentSwitch.Value == 'EPSC'
        [outputdata.alignedEvents, outputdata.meanEvent] = meanTraceMaxRise(outputdata.sponEventTraces,.025*samplerate+1,samplerate,.06,-1);
    elseif inputs.CurrentSwitch.Value == 'IPSC'
        [outputdata.alignedEvents, outputdata.meanEvent] = meanTraceMaxRise(outputdata.sponEventTraces,.025*samplerate+1,samplerate,.06,1);
    end
    
    %get putative low noise
    lowNoise = find(normCI < ciCutoff)
    lowNoiseTraces = -1.*ones(.1*samplerate+1,length(lowNoise));
    for ii = 1:length(lowNoise)
        gettinglowNoise  = 100*(ii./length(lowNoise))
        lowNoiseTraces(:,ii) = outputdata.sponData(eventTimes(lowNoise(ii))-.025*samplerate:eventTimes(lowNoise(ii))+.075*samplerate)-eventBaseline(lowNoise(ii));
    end
    
    %align and get mean
    if inputs.CurrentSwitch.Value == 'EPSC'
        [alignedLowNoise, meanLowNoise] = meanTraceMaxRise(lowNoiseTraces,.025*samplerate+1,samplerate,.06,-1);
    elseif inputs.CurrentSwitch.Value == 'IPSC'
        [alignedLowNoise, meanLowNoise] = meanTraceMaxRise(lowNoiseTraces,.025*samplerate+1,samplerate,.06,1);
    end
    
    %high noise were just all large events
    %     %get putative high noise
    %     highNoise = find(normCI > 3)
    %     highNoiseTraces = -1.*ones(.1*samplerate+1,length(highNoise));
    %     for ii = 1:length(highNoise)
    %         gettinghighNoise  = 100*(ii./length(highNoise))
    %         highNoiseTraces(:,ii) = outputdata.sponData(eventTimes(highNoise(ii))-.025*samplerate:eventTimes(highNoise(ii))+.075*samplerate)-eventBaseline(highNoise(ii));
    %     end
    %
    %      %align and get mean
    %     if inputs.CurrentSwitch.Value == 'EPSC'
    %         [alignedHighNoise, meanHighNoise] = meanTraceMaxRise(highNoiseTraces,.025*samplerate+1,samplerate,.06,-1);
    %     elseif inputs.CurrentSwitch.Value == 'IPSC'
    %         [alignedHighNoise, meanHighNoise] = meanTraceMaxRise(highNoiseTraces,.025*samplerate+1,samplerate,.06,1);
    %     end
    
    %plot events vs noise
    figure(1);
    clf
    %events
    subplot(1,4,1:3)
    hold on
    plot(outputdata.alignedEvents,'color',inputs.plotColor(9,:))
    plot(outputdata.meanEvent,'color',inputs.plotColor(3,:),'linewidth',5)
    eventAx = gca;
    setAx(eventAx);
    title('detected events and mean')
    xlim([0 length(outputdata.meanEvent)]);
    ylabel('current (pA)')
    xlabel('time (data points)')
    %noise
    subplot(1,4,4)
    hold on
    plot(alignedLowNoise,'color',inputs.plotColor(9,:))
    plot(meanLowNoise,'color',inputs.plotColor(3,:),'linewidth',5)
    setAx(gca);
    title('rejected noise and mean')
    xlim([0 length(outputdata.meanEvent)]);
    ylim(eventAx.YLim);
    ylabel('current (pA)')
    xlabel('time (data points)')

    %replot raw trace with events shown
    figure(2)
    clf
    subplot(3,1,1:2)
    hold on
    plot(t,outputdata.sponData,'color',inputs.plotColor(9,:))
    scatter(t(eventTimes(outputdata.putativeEvents)),outputdata.sponData(eventTimes(outputdata.putativeEvents)),'markeredgecolor',inputs.plotColor(1,:)) %all possible events, pre template search
    xlim([t(1) t(end)])
    xlabel('Time (s)')
    ylabel('current (pA)')
    title('raw trace')
    setAx(gca);
    subplot(3,1,3)
    hold on
    plot(t,outputdata.sponData,'color',inputs.plotColor(9,:))
    scatter(t(eventTimes(outputdata.putativeEvents)),outputdata.sponData(eventTimes(outputdata.putativeEvents)),'markeredgecolor',inputs.plotColor(1,:)) %all possible events, pre template search
    tWindow(1) = randperm(length(t),1);
    tWindow(2) = tWindow(1) + samplerate;
    xlim([t(tWindow(1)) t(tWindow(2))])
    xlabel('Time (s)')
    ylabel('current (pA)')
    title('raw trace')
    
    %get iei for frequency
    outputdata.eventiei = diff(eventTimes(outputdata.putativeEvents))./samplerate; %converts to s iei
    outputdata.eventFreq = outputdata.eventiei.^(-1); %inverse of iei is frequency in Hz
    outputdata.eventiei = outputdata.eventiei.*1000; %converts to ms
    
    %get amplitude for event
    outputdata.eventAmp = -1.*ones(1,length(outputdata.putativeEvents));
    for ii = 1:length(outputdata.putativeEvents)
        outputdata.eventAmp(ii) = outputdata.sponData(eventTimes(outputdata.putativeEvents(ii)))-eventBaseline(outputdata.putativeEvents(ii));
    end
    
    %get kinetics
    if inputs.CurrentSwitch.Value == 'EPSC'
        %20-80 rt
        [~, meanTracePeak] = min(outputdata.meanEvent);
        [outputdata.sEPSCrt, ~, ~] = findRise(outputdata.meanEvent,1,meanTracePeak,samplerate);
        
        %decay tau
        outputdata.sEPSCtau = findTau(outputdata.meanEvent,1,samplerate,'Single');
        
    elseif inputs.CurrentSwitch.Value == 'IPSC'
        %20-80 rt
        [~, meanTracePeak] = max(outputdata.meanEvent);
        [outputdata.sIPSCrt, ~, ~] = findRise(outputdata.meanEvent,1,meanTracePeak,samplerate);
        
        %decay tau
        outputdata.sIPSCtau = findTau(outputdata.meanEvent,1,samplerate,'Double');
        
    end
    
    %plot iei, amplitude, display mean values in command window
    outputdata.mAmp = median(outputdata.eventAmp);
    outputdata.mFreq = median(outputdata.eventFreq)
    
    nBins.amp = ceil(max(abs(outputdata.eventAmp)));
    nBins.iei = ceil(max(abs(outputdata.eventiei))./5);
    
    histFig = figure(4);
    histFig.Position = [1050 15 350 500];
    subplot(2,1,1)
    title('iei histogram')
    hold on
    histogram(outputdata.eventiei,nBins.iei,'EdgeColor','w','FaceColor',inputs.plotColor(7,:))
    ylabel('count')
    xlabel('iei (ms)')
    subplot(2,1,2)
    title('amplitude histogram')
    hold on
    histogram(outputdata.eventAmp,nBins.amp,'EdgeColor','w','FaceColor',inputs.plotColor(7,:))
    ylabel('count')
    xlabel('amplitude (pA)')
    
end
end

%%subfunctions
function [currentTemplate] = findTemplate(data,samplerate,inputs)
%inits
noTempEvents = 0;
bufferEnd = 0;

%find currents
for ii = 0.025*samplerate:length(data)-0.075*samplerate %get from first 25ms to 75ms before end of trace (want only whole currents)
    if ii > bufferEnd
        %init
        noise_std = [];
        bsData = [];
        
        if ii <= 0.05*samplerate %for first 50ms
            %use 50ms rolling baseline, and use first the 50ms of data as baseline for all potential events in first 50ms of data
            %get noise sd
            noise_std = getMAD(data(1:0.05*samplerate));
            
            %get spontaneous events
            if inputs.CurrentSwitch.Value == 'EPSC'
                thr = median(data(1:0.05*samplerate)) - inputs.TemplateFactorEditField.Value*noise_std;
                if data(ii) < thr
                    %baseline subtract trace
                    thisBaseline = mean(data(1:0.05*samplerate));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = min(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    %if local maximum/minimum is at either end of the vector, then it's not acutally a local maximum/minimum but rather on a slope of a real event, so move on
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noTempEvents = noTempEvents +1
                            iTempTrace(noTempEvents,:) = bsData(loc+ii-0.005*samplerate-1-0.025*samplerate:loc+ii-0.005*samplerate-1+0.075*samplerate);
                            tTempEvents(noTempEvents) = ii+loc-0.005*samplerate-1;
                            bufferEnd = ii + loc + 0.003*samplerate; %start next window 3ms after EPSC peak
                        end
                    end
                end
            elseif inputs.CurrentSwitch.Value == 'IPSC'
                thr = median(data(1:0.05*samplerate)) + inputs.TemplateFactorEditField.Value*noise_std;
                if data(ii) > thr
                    %baseline subtract trace
                    thisBaseline = mean(data(1:0.05*samplerate));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = max(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noTempEvents = noTempEvents +1
                            iTempTrace(noTempEvents,:) = bsData(loc+ii-0.005*samplerate-1-0.025*samplerate:loc+ii-0.005*samplerate-1+0.075*samplerate);
                            tTempEvents(noTempEvents) = ii+loc-0.005*samplerate-1;
                            bufferEnd = ii + loc + 0.006*samplerate; %start next window 6ms after IPSC peak
                        end
                    end
                end
            end
        else %for the rest of the trace, use the 50ms pre-data point ii as the baseline
            %get noise_std
            noise_std = getMAD(data(1:0.05*samplerate));
            
            %get spontaneous events
            if inputs.CurrentSwitch.Value == 'EPSC'
                thr = median(data(ii-0.05*samplerate:ii-1)) - inputs.TemplateFactorEditField.Value*noise_std;
                if data(ii) < thr
                    %get baseline subtracted trace
                    thisBaseline = mean(data(ii-0.05*samplerate:ii-1));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = min(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noTempEvents = noTempEvents +1
                            iTempTrace(noTempEvents,:) = bsData(loc+ii-0.005*samplerate-1-0.025*samplerate:loc+ii-0.005*samplerate-1+0.075*samplerate);
                            tTempEvents(noTempEvents) = ii+loc-0.005*samplerate-1;
                            bufferEnd = ii + loc + 0.003*samplerate; %start next window 3ms after EPSC peak
                        end
                    end
                end
            elseif inputs.CurrentSwitch.Value == 'IPSC'
                thr = median(data(ii-0.05*samplerate:ii-1)) + inputs.TemplateFactorEditField.Value*noise_std;
                if data(ii) > thr
                    %get baseline subtracted trace
                    thisBaseline = mean(data(ii-0.05*samplerate:ii-1));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = max(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noTempEvents = noTempEvents +1
                            iTempTrace(noTempEvents,:) = bsData(loc+ii-0.005*samplerate-1-0.025*samplerate:loc+ii-0.005*samplerate-1+0.075*samplerate);
                            tTempEvents(noTempEvents) = ii+loc-0.005*samplerate-1;
                            bufferEnd = ii + loc + 0.006*samplerate;
                        end
                    end
                end
            end
        end
    end
end

%get random events
% if noTempEvents >= 100
%     theseEvents = iTempTrace(randperm(noTempEvents,50),:);
%     theseEvents = theseEvents';
% else
theseEvents = iTempTrace(randperm(noTempEvents,floor(noTempEvents.*.1)),:); %use 10% of potential events for template
theseEvents = theseEvents';
% end

if inputs.CurrentSwitch.Value == 'EPSC'
    [alignedTraces, currentTemplate] = meanTraceMaxRise(theseEvents,.025*samplerate+1,samplerate,.06,-1);
elseif inputs.CurrentSwitch.Value == 'IPSC'
    [alignedTraces, currentTemplate] = meanTraceMaxRise(theseEvents,.025*samplerate+1,samplerate,.06,1);
end

%plot template
templateFig = figure(1);
templateFig.Position = [300 600 750 355];
hold on
plot(alignedTraces,'color',inputs.plotColor(ceil(size(inputs.plotColor,1)/2),:),'linewidth',1)
plot(currentTemplate,'color',inputs.plotColor(1,:),'linewidth',5)
title('template current')
ylabel('current (pA)')
xlabel('time (data points, 10kHz sample)')
xlim([1 length(currentTemplate)])
tempAx = gca;
setAx(tempAx);

end

function [tEvents,eventBaseline] = findEvents(data,samplerate,inputs)
%inits
noEvents = 0;
bufferEnd = 0;

%find currents
for ii = 0.025*samplerate:length(data)-0.075*samplerate %get from first 25ms to 75ms before end of trace (want only whole currents)
    if ii > bufferEnd
        %init
        noise_std = [];
        bsData = [];
        
        if ii <= 0.05*samplerate %for first 50ms
            %use 50ms rolling baseline, and use first the 50ms of data as baseline for all potential events in first 50ms of data
            %get noise sd
            noise_std = getMAD(data(1:0.05*samplerate));
            
            %get spontaneous events
            if inputs.CurrentSwitch.Value == 'EPSC'
                thr = median(data(1:0.05*samplerate)) - inputs.NoiseSDEditField.Value*noise_std;
                if data(ii) < thr
                    %baseline subtract trace
                    thisBaseline = mean(data(1:0.05*samplerate));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = min(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 5ms buffer for iei
                    %if local maximum/minimum is at either end of the vector, then it's not acutally a local maximum/minimum but rather on a slope of a real event, so move on
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noEvents = noEvents +1
                            tEvents(noEvents) = ii+loc-0.005*samplerate-1;
                            eventBaseline(noEvents) = thisBaseline;
                            bufferEnd = ii + loc + 0.003*samplerate; %start next window 3ms after EPSC peak
                        end
                    end
                end
            elseif inputs.CurrentSwitch.Value == 'IPSC'
                thr = median(data(1:0.05*samplerate)) + inputs.NoiseSDEditField.Value*noise_std;
                if data(ii) > thr
                    %baseline subtract trace
                    thisBaseline = mean(data(1:0.05*samplerate));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = max(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noEvents = noEvents +1
                            tEvents(noEvents) = ii+loc-0.005*samplerate-1;
                            eventBaseline(noEvents) = thisBaseline;
                            bufferEnd = ii + loc + 0.006*samplerate;
                        end
                    end
                end
            end
        else %for the rest of the trace, use the 50ms pre-data point ii as the baseline
            %get noise_std
            noise_std = getMAD(data(ii-0.05*samplerate:ii-1));
            
            %get spontaneous events
            if inputs.CurrentSwitch.Value == 'EPSC'
                thr = median(data(ii-0.05*samplerate:ii-1)) - inputs.NoiseSDEditField.Value*noise_std;
                if data(ii) < thr
                    %get baseline subtracted trace
                    thisBaseline = mean(data(ii-0.05*samplerate:ii-1));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = min(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noEvents = noEvents +1
                            tEvents(noEvents) = ii+loc-0.005*samplerate-1;
                            eventBaseline(noEvents) = thisBaseline;
                            bufferEnd = ii + loc + 0.003*samplerate; %start next window 3ms after EPSC peak
                        end
                    end
                end
            elseif inputs.CurrentSwitch.Value == 'IPSC'
                thr = median(data(ii-0.05*samplerate:ii-1)) + inputs.NoiseSDEditField.Value*noise_std;
                if data(ii) > thr
                    %get baseline subtracted trace
                    thisBaseline = mean(data(ii-0.05*samplerate:ii-1));
                    bsData = data - thisBaseline;
                    %find peak
                    [~,loc] = max(bsData(ii-0.005*samplerate:ii+0.005*samplerate)); %give 10ms buffer for spc iei
                    if loc ~= 1
                        if loc ~= .01*samplerate+1
                            noEvents = noEvents +1
                            tEvents(noEvents) = ii+loc-0.005*samplerate-1;
                            eventBaseline(noEvents) = thisBaseline;
                            bufferEnd = ii + loc + 0.006*samplerate;
                        end
                    end
                end
            end
        end
    end
end
end
