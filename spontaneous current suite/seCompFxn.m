%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% spon event Grp Comp %%%%%%%%%%%%
%%%%%%%%%%% Created: 10-20-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 01-06-2019 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputdata] = seCompFxn(inputs)
close all

%% inits
samplerate = inputs.SamplerateHzEditField.Value;
ctrlColor = [.3333 .8039 1];
expColor = inputs.plotColor(end,:);
lightgray=[.75 .75 .75]; %light gray for individual data points

%% LOAD DATA
%ctrl
if strcmp(inputs.CurrentSwitch_2.Value,'EPSC')
    cd([inputs.spscDir,'/',inputs.ControlGroupblueEditField.Value,'/sEPSC']);
elseif strcmp(inputs.CurrentSwitch_2.Value,'IPSC')
    cd([inputs.spscDir,'/',inputs.ControlGroupblueEditField.Value,'/sIPSC']);
end
contents = dir('*.mat');
filenames = {contents.name}';
spscFiles.(inputs.ControlGroupblueEditField.Value) = fullfile(cd,filenames);

%load
for ii = 1:length(spscFiles.(inputs.ControlGroupblueEditField.Value))
    load(spscFiles.(inputs.ControlGroupblueEditField.Value){ii})
    outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)(ii) = abs(app.spscData.mAmp);
    outputdata.frequency.(inputs.ControlGroupblueEditField.Value)(ii) = app.spscData.mFreq;
    if strcmp(inputs.CurrentSwitch_2.Value,'EPSC')
        outputdata.rise.(inputs.ControlGroupblueEditField.Value)(ii) = app.spscData.sEPSCrt;
        outputdata.decay.(inputs.ControlGroupblueEditField.Value)(ii) = app.spscData.sEPSCtau;
    elseif strcmp(inputs.CurrentSwitch_2.Value,'IPSC')
        outputdata.rise.(inputs.ControlGroupblueEditField.Value)(ii) = app.spscData.sIPSCrt;
        outputdata.decay.(inputs.ControlGroupblueEditField.Value)(ii) = app.spscData.sIPSCtau;
    end
    outputdata.meanEventTrace.(inputs.ControlGroupblueEditField.Value)(ii,:) = app.spscData.meanEvent;
    
    %example trace
    if ii == inputs.CtrlgroupexamplecellEditField.Value
        outputdata.eventTraces.(inputs.ControlGroupblueEditField.Value) = app.spscData.alignedEvents;
        outputdata.ieiDist.(inputs.ControlGroupblueEditField.Value) = app.spscData.eventiei;
        outputdata.frequencyDist.(inputs.ControlGroupblueEditField.Value) = app.spscData.eventFreq;
        outputdata.rawTrace.(inputs.ControlGroupblueEditField.Value) = app.spscData.sponData;
    end
    clear app
end

%Experimental Group
if strcmp(inputs.CurrentSwitch_2.Value,'EPSC')
    cd([inputs.spscDir,'/',inputs.ExperimentalGroupEditField.Value,'/sEPSC']);
elseif strcmp(inputs.CurrentSwitch_2.Value,'IPSC')
    cd([inputs.spscDir,'/',inputs.ExperimentalGroupEditField.Value,'/sIPSC']);
end
contents = dir('*.mat');
filenames = {contents.name}';
spscFiles.(inputs.ExperimentalGroupEditField.Value) = fullfile(cd,filenames);

%load
for ii = 1:length(spscFiles.(inputs.ExperimentalGroupEditField.Value))
    load(spscFiles.(inputs.ExperimentalGroupEditField.Value){ii})
    outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)(ii) = abs(app.spscData.mAmp);
    outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)(ii) = app.spscData.mFreq;
    if strcmp(inputs.CurrentSwitch_2.Value,'EPSC')
        outputdata.rise.(inputs.ExperimentalGroupEditField.Value)(ii) = app.spscData.sEPSCrt;
        outputdata.decay.(inputs.ExperimentalGroupEditField.Value)(ii) = app.spscData.sEPSCtau;
    elseif strcmp(inputs.CurrentSwitch_2.Value,'IPSC')
        outputdata.rise.(inputs.ExperimentalGroupEditField.Value)(ii) = app.spscData.sIPSCrt;
        outputdata.decay.(inputs.ExperimentalGroupEditField.Value)(ii) = app.spscData.sIPSCtau;
    end
    outputdata.meanEventTrace.(inputs.ExperimentalGroupEditField.Value)(ii,:) = app.spscData.meanEvent;
    
    %example trace
    if ii == inputs.CtrlgroupexamplecellEditField.Value
        outputdata.eventTraces.(inputs.ExperimentalGroupEditField.Value) = app.spscData.alignedEvents;
        outputdata.ieiDist.(inputs.ExperimentalGroupEditField.Value) = app.spscData.eventiei;
        outputdata.frequencyDist.(inputs.ExperimentalGroupEditField.Value) = app.spscData.eventFreq;
        outputdata.rawTrace.(inputs.ExperimentalGroupEditField.Value) = app.spscData.sponData;
    end
    clear app
end

%time vector
t = samplerate^-1:samplerate^-1:length(outputdata.rawTrace.(inputs.ControlGroupblueEditField.Value))/samplerate;
%% ANALYSIS
%normality check
if length(spscFiles.(inputs.ControlGroupblueEditField.Value)) > 3
    [normalityAmp.(inputs.ControlGroupblueEditField.Value),~] = adtest(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value));
    [normalityFreq.(inputs.ControlGroupblueEditField.Value),~] = adtest(outputdata.frequency.(inputs.ControlGroupblueEditField.Value));
    [normalityRise.(inputs.ControlGroupblueEditField.Value),~] = adtest(outputdata.rise.(inputs.ControlGroupblueEditField.Value));
    [normalityDecay.(inputs.ControlGroupblueEditField.Value),~] = adtest(outputdata.decay.(inputs.ControlGroupblueEditField.Value));
else %assume normal for testing purposes so code can run
    normalityAmp.(inputs.ControlGroupblueEditField.Value) = 0;
    normalityFreq.(inputs.ControlGroupblueEditField.Value) = 0;
    normalityRise.(inputs.ControlGroupblueEditField.Value) = 0;
    normalityDecay.(inputs.ControlGroupblueEditField.Value) = 0;
end
if length(spscFiles.(inputs.ExperimentalGroupEditField.Value)) > 3
    [normalityAmp.(inputs.ExperimentalGroupEditField.Value),~] = adtest(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value));
    [normalityFreq.(inputs.ExperimentalGroupEditField.Value),~] = adtest(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value));
    [normalityRise.(inputs.ExperimentalGroupEditField.Value),~] = adtest(outputdata.rise.(inputs.ExperimentalGroupEditField.Value));
    [normalityDecay.(inputs.ExperimentalGroupEditField.Value),~] = adtest(outputdata.decay.(inputs.ExperimentalGroupEditField.Value));
else
    normalityAmp.(inputs.ExperimentalGroupEditField.Value) = 0;
    normalityFreq.(inputs.ExperimentalGroupEditField.Value) = 0;
    normalityRise.(inputs.ExperimentalGroupEditField.Value) = 0;
    normalityDecay.(inputs.ExperimentalGroupEditField.Value) = 0;
end


%effect, amplitude
if sum([normalityAmp.(inputs.ControlGroupblueEditField.Value) normalityAmp.(inputs.ExperimentalGroupEditField.Value)]) == 0
    %if data are normal
    [~,pAmp_tt] = ttest2(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),...
        outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value))
else
    %if data are non-normal
    [pAmp_mwu,~] = ranksum(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),...
        outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)) %#ok<*ASGLU,*NOPRT>
end
%effect, frequency
if sum([normalityFreq.(inputs.ControlGroupblueEditField.Value) normalityFreq.(inputs.ExperimentalGroupEditField.Value)]) == 0
    %if data are normal
    [~,pFreq_tt] = ttest2(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),...
        outputdata.frequency.(inputs.ExperimentalGroupEditField.Value))
else
    %if data are non-normal
    [pFreq_mwu,~] = ranksum(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),...
        outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)) %#ok<*ASGLU,*NOPRT>
end
%effect, rise
if sum([normalityRise.(inputs.ControlGroupblueEditField.Value) normalityRise.(inputs.ExperimentalGroupEditField.Value)]) == 0
    %if data are normal
    [~,pRise_tt] = ttest2(outputdata.rise.(inputs.ControlGroupblueEditField.Value),...
        outputdata.rise.(inputs.ExperimentalGroupEditField.Value))
else
    %if data are non-normal
    [pRise_mwu,~] = ranksum(outputdata.rise.(inputs.ControlGroupblueEditField.Value),...
        outputdata.rise.(inputs.ExperimentalGroupEditField.Value)) %#ok<*ASGLU,*NOPRT>
end
%effect, decay
if sum([normalityDecay.(inputs.ControlGroupblueEditField.Value) normalityDecay.(inputs.ExperimentalGroupEditField.Value)]) == 0
    %if data are normal
    [~,pDecay_tt] = ttest2(outputdata.decay.(inputs.ControlGroupblueEditField.Value),...
        outputdata.decay.(inputs.ExperimentalGroupEditField.Value))
else
    %if data are non-normal
    [pDecay_mwu,~] = ranksum(outputdata.decay.(inputs.ControlGroupblueEditField.Value),...
        outputdata.decay.(inputs.ExperimentalGroupEditField.Value)) %#ok<*ASGLU,*NOPRT>
end

%get central tendency values for amplitude
if sum([normalityAmp.(inputs.ControlGroupblueEditField.Value) normalityAmp.(inputs.ExperimentalGroupEditField.Value)]) == 0
    outputdata.mAmp.(inputs.ControlGroupblueEditField.Value) = mean(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value));
    outputdata.errAmp.(inputs.ControlGroupblueEditField.Value) = sem(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),1);
    outputdata.mAmp.(inputs.ExperimentalGroupEditField.Value) = mean(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value));
    outputdata.errAmp.(inputs.ExperimentalGroupEditField.Value) = sem(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),1);
    outputdata.ampNormality = 0; %data are normal
else
    outputdata.mAmp.(inputs.ControlGroupblueEditField.Value) = median(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value));
    outputdata.negErrAmp.(inputs.ControlGroupblueEditField.Value) = median(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)) - quantile(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),.25);
    outputdata.posErrAmp.(inputs.ControlGroupblueEditField.Value) = quantile(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),.75) - median(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value));
    outputdata.mAmp.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value));
    outputdata.negErrAmp.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)) - quantile(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),.25);
    outputdata.posErrAmp.(inputs.ExperimentalGroupEditField.Value) = quantile(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),.75) - median(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value));
    outputdata.ampNormality = 1; %data are non-normal
end

%get central tendency values for frequency
if sum([normalityFreq.(inputs.ControlGroupblueEditField.Value) normalityFreq.(inputs.ExperimentalGroupEditField.Value)]) == 0
    outputdata.mFreq.(inputs.ControlGroupblueEditField.Value) = mean(outputdata.frequency.(inputs.ControlGroupblueEditField.Value));
    outputdata.errFreq.(inputs.ControlGroupblueEditField.Value) = sem(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),1);
    outputdata.mFreq.(inputs.ExperimentalGroupEditField.Value) = mean(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value));
    outputdata.errFreq.(inputs.ExperimentalGroupEditField.Value) = sem(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),1);
    outputdata.freqNormality = 0; %data are normal
else
    outputdata.mFreq.(inputs.ControlGroupblueEditField.Value) = median(outputdata.frequency.(inputs.ControlGroupblueEditField.Value));
    outputdata.negErrFreq.(inputs.ControlGroupblueEditField.Value) = median(outputdata.frequency.(inputs.ControlGroupblueEditField.Value)) - quantile(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),.25);
    outputdata.posErrFreq.(inputs.ControlGroupblueEditField.Value) = quantile(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),.75) - median(outputdata.frequency.(inputs.ControlGroupblueEditField.Value));
    outputdata.mFreq.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value));
    outputdata.negErrFreq.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)) - quantile(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),.25);
    outputdata.posErrFreq.(inputs.ExperimentalGroupEditField.Value) = quantile(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),.75) - median(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value));
    outputdata.freqNormality = 1; %data are non-normal
end

%get central tendency values for rise
if sum([normalityRise.(inputs.ControlGroupblueEditField.Value) normalityRise.(inputs.ExperimentalGroupEditField.Value)]) == 0
    outputdata.mRise.(inputs.ControlGroupblueEditField.Value) = mean(outputdata.rise.(inputs.ControlGroupblueEditField.Value));
    outputdata.errRise.(inputs.ControlGroupblueEditField.Value) = sem(outputdata.rise.(inputs.ControlGroupblueEditField.Value),1);
    outputdata.mRise.(inputs.ExperimentalGroupEditField.Value) = mean(outputdata.rise.(inputs.ExperimentalGroupEditField.Value));
    outputdata.errRise.(inputs.ExperimentalGroupEditField.Value) = sem(outputdata.rise.(inputs.ExperimentalGroupEditField.Value),1);
    outputdata.riseNormality = 0; %data are normal
else
    outputdata.mRise.(inputs.ControlGroupblueEditField.Value) = median(outputdata.rise.(inputs.ControlGroupblueEditField.Value));
    outputdata.negErrRise.(inputs.ControlGroupblueEditField.Value) = median(outputdata.rise.(inputs.ControlGroupblueEditField.Value)) - quantile(outputdata.rise.(inputs.ControlGroupblueEditField.Value),.25);
    outputdata.posErrRise.(inputs.ControlGroupblueEditField.Value) = quantile(outputdata.rise.(inputs.ControlGroupblueEditField.Value),.75) - median(outputdata.rise.(inputs.ControlGroupblueEditField.Value));
    outputdata.mRise.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.rise.(inputs.ExperimentalGroupEditField.Value));
    outputdata.negErrRise.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.rise.(inputs.ExperimentalGroupEditField.Value)) - quantile(outputdata.rise.(inputs.ExperimentalGroupEditField.Value),.25);
    outputdata.posErrRise.(inputs.ExperimentalGroupEditField.Value) = quantile(outputdata.rise.(inputs.ExperimentalGroupEditField.Value),.75) - median(outputdata.rise.(inputs.ExperimentalGroupEditField.Value));
    outputdata.riseNormality = 1; %data are non-normal
end

%get central tendency values for decay
if sum([normalityDecay.(inputs.ControlGroupblueEditField.Value) normalityDecay.(inputs.ExperimentalGroupEditField.Value)]) == 0
    outputdata.mDecay.(inputs.ControlGroupblueEditField.Value) = mean(outputdata.decay.(inputs.ControlGroupblueEditField.Value));
    outputdata.errDecay.(inputs.ControlGroupblueEditField.Value) = sem(outputdata.decay.(inputs.ControlGroupblueEditField.Value),1);
    outputdata.mDecay.(inputs.ExperimentalGroupEditField.Value) = mean(outputdata.decay.(inputs.ExperimentalGroupEditField.Value));
    outputdata.errDecay.(inputs.ExperimentalGroupEditField.Value) = sem(outputdata.decay.(inputs.ExperimentalGroupEditField.Value),1);
    outputdata.decayNormality = 0; %data are normal
else
    outputdata.mDecay.(inputs.ControlGroupblueEditField.Value) = median(outputdata.decay.(inputs.ControlGroupblueEditField.Value));
    outputdata.negErrDecay.(inputs.ControlGroupblueEditField.Value) = median(outputdata.decay.(inputs.ControlGroupblueEditField.Value)) - quantile(outputdata.decay.(inputs.ControlGroupblueEditField.Value),.25);
    outputdata.posErrDecay.(inputs.ControlGroupblueEditField.Value) = quantile(outputdata.decay.(inputs.ControlGroupblueEditField.Value),.75) - median(outputdata.decay.(inputs.ControlGroupblueEditField.Value));
    outputdata.mDecay.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.decay.(inputs.ExperimentalGroupEditField.Value));
    outputdata.negErrDecay.(inputs.ExperimentalGroupEditField.Value) = median(outputdata.decay.(inputs.ExperimentalGroupEditField.Value)) - quantile(outputdata.decay.(inputs.ExperimentalGroupEditField.Value),.25);
    outputdata.posErrDecay.(inputs.ExperimentalGroupEditField.Value) = quantile(outputdata.decay.(inputs.ExperimentalGroupEditField.Value),.75) - median(outputdata.decay.(inputs.ExperimentalGroupEditField.Value));
    outputdata.decayNormality = 1; %data are non-normal
end

%% DISPLAY DATA
%amplitude
ampFig = figure(1);
ampFig.Position = [500 485 150 325];
hold on
scatter(ones(length(spscFiles.(inputs.ControlGroupblueEditField.Value)),1),...
    outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),75,lightgray,'filled');
scatter(2.*ones(length(spscFiles.(inputs.ExperimentalGroupEditField.Value)),1),...
    outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),75,lightgray,'filled');
if sum([normalityAmp.(inputs.ControlGroupblueEditField.Value) normalityAmp.(inputs.ExperimentalGroupEditField.Value)]) == 0
    errorbar(1,mean(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)),...
        sem(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),2),...
        'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,mean(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,mean(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)),...
        sem(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),2),...
        'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,mean(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
else
    errorbar(1,median(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)),...
        outputdata.negErrAmp.(inputs.ControlGroupblueEditField.Value),outputdata.posErrAmp.(inputs.ControlGroupblueEditField.Value),'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,median(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,median(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)),...
        outputdata.negErrAmp.(inputs.ExperimentalGroupEditField.Value),outputdata.posErrAmp.(inputs.ExperimentalGroupEditField.Value),'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,median(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
end
ampAx = gca;
setAx(ampAx);
xlim([.5 2.5]);
ylim([0 max([outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value) ...
    outputdata.amplitude.(inputs.ControlGroupblueEditField.Value)])*1.15])
ylabel('amplitude (pA)')
ampAx.XTick = [1 2];
ampAx.XTickLabel{1} = (inputs.ControlGroupblueEditField.Value);
ampAx.XTickLabel{2} = (inputs.ExperimentalGroupEditField.Value);

%frequency
freqFig = figure(2);
freqFig.Position = [675 485 150 325];
hold on
scatter(ones(length(spscFiles.(inputs.ControlGroupblueEditField.Value)),1),...
    outputdata.frequency.(inputs.ControlGroupblueEditField.Value),75,lightgray,'filled');
scatter(2.*ones(length(spscFiles.(inputs.ExperimentalGroupEditField.Value)),1),...
    outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),75,lightgray,'filled');
if sum([normalityFreq.(inputs.ControlGroupblueEditField.Value) normalityFreq.(inputs.ExperimentalGroupEditField.Value)]) == 0
    errorbar(1,mean(outputdata.frequency.(inputs.ControlGroupblueEditField.Value)),...
        sem(outputdata.frequency.(inputs.ControlGroupblueEditField.Value),2),...
        'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,mean(outputdata.frequency.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,mean(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)),...
        sem(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),2),...
        'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,mean(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
else
    errorbar(1,median(outputdata.frequency.(inputs.ControlGroupblueEditField.Value)),...
        outputdata.negErrFreq.(inputs.ControlGroupblueEditField.Value),outputdata.posErrFreq.(inputs.ControlGroupblueEditField.Value),'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,median(outputdata.frequency.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,median(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)),...
        outputdata.negErrFreq.(inputs.ExperimentalGroupEditField.Value),outputdata.posErrFreq.(inputs.ExperimentalGroupEditField.Value),'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,median(outputdata.frequency.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
end
freqAx = gca;
setAx(freqAx);
xlim([.5 2.5]);
ylim([0 max([outputdata.frequency.(inputs.ExperimentalGroupEditField.Value) ...
    outputdata.frequency.(inputs.ControlGroupblueEditField.Value)])*1.15])
ylabel('frequency (Hz)')
freqAx.XTick = [1 2];
freqAx.XTickLabel{1} = (inputs.ControlGroupblueEditField.Value);
freqAx.XTickLabel{2} = (inputs.ExperimentalGroupEditField.Value);

%freq vs amplitude plot
faFig = figure(5);
faFig.Position = [500 85 325 300];
hold on
scatter(outputdata.amplitude.(inputs.ControlGroupblueEditField.Value),...
    outputdata.frequency.(inputs.ControlGroupblueEditField.Value),75,ctrlColor,'filled')
scatter(outputdata.amplitude.(inputs.ExperimentalGroupEditField.Value),...
    outputdata.frequency.(inputs.ExperimentalGroupEditField.Value),75,expColor,'filled')
faAx = gca;
setAx(faAx);
if faAx.XLim(2) >= faAx.YLim(2)
    faAx.XLim(1) = 0;
    ylim([0 faAx.XLim(2)])
else
    faAx.YLim(1) = 0;
    xlim([0 faAx.YLim(2)])
end
xlabel('Amplitude (pA)')
ylabel('Frequency (Hz)')

%rise
riseFig = figure(3);
riseFig.Position = [850 485 150 325];
hold on
scatter(ones(length(spscFiles.(inputs.ControlGroupblueEditField.Value)),1),...
    outputdata.rise.(inputs.ControlGroupblueEditField.Value),75,lightgray,'filled');
scatter(2.*ones(length(spscFiles.(inputs.ExperimentalGroupEditField.Value)),1),...
    outputdata.rise.(inputs.ExperimentalGroupEditField.Value),75,lightgray,'filled');
if sum([normalityRise.(inputs.ControlGroupblueEditField.Value) normalityRise.(inputs.ExperimentalGroupEditField.Value)]) == 0
    errorbar(1,mean(outputdata.rise.(inputs.ControlGroupblueEditField.Value)),...
        sem(outputdata.rise.(inputs.ControlGroupblueEditField.Value),2),...
        'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,mean(outputdata.rise.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,mean(outputdata.rise.(inputs.ExperimentalGroupEditField.Value)),...
        sem(outputdata.rise.(inputs.ExperimentalGroupEditField.Value),2),...
        'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,mean(outputdata.rise.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
else
    errorbar(1,median(outputdata.rise.(inputs.ControlGroupblueEditField.Value)),...
        outputdata.negErrRise.(inputs.ControlGroupblueEditField.Value),outputdata.posErrRise.(inputs.ControlGroupblueEditField.Value),'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,median(outputdata.rise.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,median(outputdata.rise.(inputs.ExperimentalGroupEditField.Value)),...
        outputdata.negErrRise.(inputs.ExperimentalGroupEditField.Value),outputdata.posErrRise.(inputs.ExperimentalGroupEditField.Value),'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,median(outputdata.rise.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
end
riseAx = gca;
setAx(riseAx);
xlim([.5 2.5]);
ylim([0 max([outputdata.rise.(inputs.ExperimentalGroupEditField.Value) ...
    outputdata.rise.(inputs.ControlGroupblueEditField.Value)])*1.15])
ylabel('20-80% risetime (ms)')
riseAx.XTick = [1 2];
riseAx.XTickLabel{1} = (inputs.ControlGroupblueEditField.Value);
riseAx.XTickLabel{2} = (inputs.ExperimentalGroupEditField.Value);

%decay
decayFig = figure(4);
decayFig.Position = [1025 485 150 325];
hold on
scatter(ones(length(spscFiles.(inputs.ControlGroupblueEditField.Value)),1),...
    outputdata.decay.(inputs.ControlGroupblueEditField.Value),75,lightgray,'filled');
scatter(2.*ones(length(spscFiles.(inputs.ExperimentalGroupEditField.Value)),1),...
    outputdata.decay.(inputs.ExperimentalGroupEditField.Value),75,lightgray,'filled');
if sum([normalityDecay.(inputs.ControlGroupblueEditField.Value) normalityDecay.(inputs.ExperimentalGroupEditField.Value)]) == 0
    errorbar(1,mean(outputdata.decay.(inputs.ControlGroupblueEditField.Value)),...
        sem(outputdata.decay.(inputs.ControlGroupblueEditField.Value),2),...
        'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,mean(outputdata.decay.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,mean(outputdata.decay.(inputs.ExperimentalGroupEditField.Value)),...
        sem(outputdata.decay.(inputs.ExperimentalGroupEditField.Value),2),...
        'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,mean(outputdata.decay.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
else
    errorbar(1,median(outputdata.decay.(inputs.ControlGroupblueEditField.Value)),...
        outputdata.negErrDecay.(inputs.ControlGroupblueEditField.Value),outputdata.posErrDecay.(inputs.ControlGroupblueEditField.Value),'color',ctrlColor,'linewidth',5,'CapSize',0)
    scatter(1,median(outputdata.decay.(inputs.ControlGroupblueEditField.Value)),250,ctrlColor,'filled')
    errorbar(2,median(outputdata.decay.(inputs.ExperimentalGroupEditField.Value)),...
        outputdata.negErrDecay.(inputs.ExperimentalGroupEditField.Value),outputdata.posErrDecay.(inputs.ExperimentalGroupEditField.Value),'color',expColor,'linewidth',5,'CapSize',0)
    scatter(2,median(outputdata.decay.(inputs.ExperimentalGroupEditField.Value)),250,expColor,'filled')
end
decayAx = gca;
setAx(decayAx);
xlim([.5 2.5]);
ylim([0 max([outputdata.decay.(inputs.ExperimentalGroupEditField.Value) ...
    outputdata.decay.(inputs.ControlGroupblueEditField.Value)])*1.15])
ylabel('\tau-decay (ms)')
decayAx.XTick = [1 2];
decayAx.XTickLabel{1} = (inputs.ControlGroupblueEditField.Value);
decayAx.XTickLabel{2} = (inputs.ExperimentalGroupEditField.Value);

%example traces
exTrace;
