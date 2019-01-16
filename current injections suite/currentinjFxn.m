function [output] = currentinjFxn(inputs)
%for analysis of biophysics differences bw responding and nonresponding neurons in minstim data
%created 11-27-18 %modified 12-03-18
close all

%% INIT VARS
%sample rate
samplerate=inputs.SamplerateHzEditField.Value;

%stimulus start in ms
injstart = inputs.PulsestartmsEditField.Value*(samplerate/1000);
injend = inputs.PulseendmsEditField.Value*(samplerate/1000);

%baseline
blstart=injstart-inputs.BaselinedurationmsEditField.Value*(samplerate/1000)-1*(samplerate/1000);
blend=injstart-1*(samplerate/1000);

%% LOAD DATA
cd(inputs.ciDir);

%data
contents = dir('*.abf');
filenames = {contents.name}';
ciFiles = fullfile(cd,filenames);

%number of files
noRuns = str2double(inputs.NumberofrunsDropDown.Value);
if noRuns == 2
    if inputs.InjectiondeltapAEditField.Value <= inputs.InjectiondeltapAEditField_2.Value
        runOrder = [1 2];
    else
        runOrder = [2 1];
    end
end

%% PREPARE DATA
if noRuns == 1
    output.ccFile = abfload(ciFiles{1},'sweeps','a');
    %creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
    %M is dependent on sample rate
    %first data point in pCFile is at t=0
    noSweeps=inputs.NumberofsweepsEditField.Value; %gives number of sweeps per file
    t = 1/samplerate:1/samplerate:size(output.ccFile,1)/samplerate;
    
    %find dI
    sweepDelta=zeros(noSweeps,2);
    sweepDelta(:,1)=1:1:noSweeps;
    sweepDelta(:,2)=(inputs.StartinginjectionpAEditField.Value:inputs.InjectiondeltapAEditField.Value:...
        (inputs.InjectiondeltapAEditField.Value*(noSweeps-1)+inputs.StartinginjectionpAEditField.Value));
    [~,sweepTau]=min(abs(sweepDelta(:,2)+50)); %finds sweep index closest to -50pA sweep -- for tau-m
    
    %find AP Threshold
    [thVals,thTimes]=findSpkThresh(output.ccFile,injstart,injend,samplerate,sweepDelta(:,2),inputs.plotColor);
    thTimes(thTimes~=-1) = thTimes(thTimes~=-1) +1;
    rheoSweep=find(thVals(:,1)~=-1);
    if length(rheoSweep) >1
        rheoSweep(2:end)=[];
        output.rheo = sweepDelta(rheoSweep,2);
    elseif isempty(rheoSweep)
        rheoSweep=noSweeps+1;
    end
elseif noRuns == 2
    output.ccFile.RunOne = abfload(ciFiles{1},'sweeps','a');
    output.ccFile.RunTwo = abfload(ciFiles{2},'sweeps','a');
    %creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
    %M is dependent on sample rate
    %first data point in pCFile is at t=0
    noSweeps.RunOne=inputs.NumberofsweepsEditField.Value; %gives number of sweeps per file
    noSweeps.RunTwo=inputs.NumberofsweepsEditField_2.Value; %gives number of sweeps per file
    t = 1/samplerate:1/samplerate:size(output.ccFile.RunOne,1)/samplerate;
    
    %find dI
    sweepDelta.RunOne=zeros(noSweeps.RunOne,2);
    sweepDelta.RunOne(:,1)=1:1:noSweeps.RunOne;
    sweepDelta.RunOne(:,2)=(inputs.StartinginjectionpAEditField.Value:inputs.InjectiondeltapAEditField.Value:...
        (inputs.InjectiondeltapAEditField.Value*(noSweeps.RunOne-1)+inputs.StartinginjectionpAEditField.Value));
    sweepDelta.RunTwo=zeros(noSweeps.RunTwo,2);
    sweepDelta.RunTwo(:,1)=1:1:noSweeps.RunTwo;
    sweepDelta.RunTwo(:,2)=(inputs.StartinginjectionpAEditField_2.Value:inputs.InjectiondeltapAEditField_2.Value:...
        (inputs.InjectiondeltapAEditField_2.Value*(noSweeps.RunTwo-1)+inputs.StartinginjectionpAEditField_2.Value));
    [~,sweepTau]=min(abs(sweepDelta.RunTwo(:,2)+50)); %finds sweep index closest to -50pA sweep -- for tau-m
    
    %find AP Threshold
    [thVals.RunOne,thTimes.RunOne]=findSpkThresh(output.ccFile.RunOne,injstart,injend,samplerate,sweepDelta.RunOne(:,2),inputs.plotColor);
    thTimes.RunOne(thTimes.RunOne~=-1) = thTimes.RunOne(thTimes.RunOne~=-1) +1;
    rheoSweep.RunOne=find(thVals.RunOne(:,1)~=-1);
    if length(rheoSweep.RunOne) >1
        rheoSweep.RunOne(2:end)=[];
        output.rheo = sweepDelta.RunOne(rheoSweep.RunOne,2);
    elseif isempty(rheoSweep.RunOne)
        rheoSweep.RunOne=noSweeps.RunOne+1;
    end
    [thVals.RunTwo,thTimes.RunTwo]=findSpkThresh(output.ccFile.RunTwo,injstart,injend,samplerate,sweepDelta.RunTwo(:,2),inputs.plotColor);
    thTimes.RunTwo(thTimes.RunTwo~=-1) = thTimes.RunTwo(thTimes.RunTwo~=-1) +1;
    rheoSweep.RunTwo=find(thVals.RunTwo(:,1)~=-1);
    if length(rheoSweep.RunTwo) >1
        rheoSweep.RunTwo(2:end)=[];
        output.rheo = sweepDelta.RunTwo(rheoSweep.RunTwo,2);
    elseif isempty(rheoSweep.RunTwo)
        rheoSweep.RunTwo=noSweeps.RunTwo+1;
    end
end

%% Display Traces
ciFig = figure(1);
title('Voltage Responses to Current Injections')
if noRuns == 1
    ciFig.Position = [560 475 300 325];
    
    %get color range
    colorShift = 1-max(inputs.plotColot);
    for ii = 1:3
        cRange(:,ii) = inputs.plotColor(ii):colorShift/noSweeps:inputs.plotColor(ii)+colorShift;
    end
    
    %plot
    subplot(3,1,[1 2])
    hold on
    for ii = 1:size(output.ccFile,3)
        plot(t,output.ccFile(:,1,ii),'color',cRange(ii,:))
    end
    xlim([.5 1.5])
    ylabel('Vm (mV)')
    setAx(gca);
    subplot(3,1,3)
    hold on
    for ii = 1:size(output.ccFile,3)
        plot(t,output.ccFile(:,2,ii),'color',cRange(ii,:))
    end
    xlim([.5 1.5])
    xlabel('time (s)')
    ylabel('injected current (pA)')
    setAx(gca);
elseif noRuns == 2
    ciFig.Position = [560 475 550 325];
    
    %get color range
    colorShift = 1-max(inputs.plotColor);
    for ii = 1:3
        cRange.RunOne(:,ii) = inputs.plotColor(ii):colorShift/noSweeps.RunOne:inputs.plotColor(ii)+colorShift;
        cRange.RunTwo(:,ii) = inputs.plotColor(ii):colorShift/noSweeps.RunTwo:inputs.plotColor(ii)+colorShift;
    end
    
    subplot(3,2,[1 3])
    hold on
    for ii = 1:size(output.ccFile.RunOne,3)
        plot(t,output.ccFile.RunOne(:,1,ii),'color',cRange.RunOne(ii,:),'linewidth',1.25)
    end
    xlim([.5 1.5])
    ylabel('Vm (mV)')
    setAx(gca);
    subplot(3,2,5)
    hold on
    for ii = 1:size(output.ccFile.RunOne,3)
        plot(t,output.ccFile.RunOne(:,2,ii),'color',cRange.RunOne(ii,:),'linewidth',1.25)
    end
    xlim([.5 1.5])
    xlabel('time (s)')
    ylabel('injected current (pA)')
    setAx(gca);
    subplot(3,2,[2 4])
    hold on
    for ii = 1:size(output.ccFile.RunTwo,3)
        plot(t,output.ccFile.RunTwo(:,1,ii),'color',cRange.RunTwo(ii,:),'linewidth',1.25)
    end
    xlim([.5 1.5])
    setAx(gca);
    subplot(3,2,6)
    hold on
    for ii = 1:size(output.ccFile.RunTwo,3)
        plot(t,output.ccFile.RunTwo(:,2,ii),'color',cRange.RunTwo(ii,:),'linewidth',1.25)
    end
    xlim([.5 1.5])
    xlabel('time (s)')
    setAx(gca);
end


%% GET MEMBRANE PROPERTIES
%membrane resistance, get with low delta experiment if two runs
if noRuns == 1
    output.rm=getRm(output.ccFile,blstart,blend,injstart,samplerate,rheoSweep,sweepDelta(:,2),cRange);
elseif noRuns == 2
    output.rm=getRm(output.ccFile.RunOne,blstart,blend,injstart,samplerate,rheoSweep.RunOne,sweepDelta.RunOne(:,2),cRange.RunOne);
end

%tau, taken at sweep closest to -50 in run with largest steps
if noRuns == 1
    output.mtau=getTau(output.ccFile,sweepTau,injstart,samplerate,cRange);
elseif noRuns == 2
    output.mtau=getTau(output.ccFile.RunTwo,sweepTau,injstart,samplerate,cRange.RunTwo);
end

%voltage sag and rebound spikes, taken at most negative sweep
if noRuns == 1
    output.vsag = getSag(output.ccFile,1,injstart,injend,blstart,blend,samplerate,cRange);
    [~, ~, output.nRbnd] = getRbndspks(output.ccFile,1,injend+1,injend+.5*samplerate,cRange,samplerate);
elseif noRuns == 2
    output.vsag = getSag(output.ccFile.RunTwo,1,injstart,injend,blstart,blend,samplerate,cRange.RunTwo);
    [~, ~, output.nRbnd] = getRbndspks(output.ccFile.RunTwo,1,injend+1,injend+.5*samplerate,cRange.RunTwo,samplerate);
end

%Max firing rate, defined as inverse of ISI during first 200ms of last current injection before reduction in AP firing occurs
if noRuns == 1
    output.maxFR = getMaxFiring(output.ccFile,noSweeps,rheoSweep,thTimes,samplerate,injstart,cRange);
elseif noRuns == 2
    output.maxFR = getMaxFiring(output.ccFile.RunTwo,noSweeps.RunTwo,rheoSweep.RunTwo,thTimes.RunTwo,samplerate,injstart,cRange.RunTwo);
end

%F-I plot
if noRuns == 1
    output.FIcruve = getFICurve(output.ccFile,noSweeps,rheoSweep,sweepDelta,thTimes,samplerate,injstart,cRange);
elseif noRuns == 2
    output.FIcurve = getFICurve(output.ccFile.RunTwo,noSweeps.RunTwo,rheoSweep.RunTwo,sweepDelta.RunTwo,thTimes.RunTwo,samplerate,injstart,cRange.RunTwo);
end

%rheobase data
%set up figure
rheoFig = figure(7);
rheoFig.Position = [130 220 675 180];

%spike train
subplot(1,8,1:4)
hold on
if noRuns == 1
    plot(t,output.ccFile(:,1,rheoSweep),'linewidth',2,'color',cRange(rheoSweep,:))
elseif noRuns == 2
    plot(t,output.ccFile.RunTwo(:,1,rheoSweep.RunTwo),'linewidth',2,'color',cRange.RunTwo(rheoSweep.RunTwo,:))
end
setAx(gca)
xlim([injstart./samplerate-.05 injend./samplerate+.1])
xlabel('time (s)')
ylabel('Vm (mV)')
title('rheobase action potentials')

%spike overlays
subplot(1,8,5:6)
hold on
if noRuns == 1
    tRheoSpikes = thTimes(rheoSweep,thTimes(rheoSweep,:)>-1);
    output.rWaves=findAPwaveforms(output.ccFile,tRheoSpikes,rheoSweep,samplerate); %aligned to threshold
    for ii = 1:size(output.rWaves,1)
        plot(1000/samplerate.*(1:size(output.rWaves,2)),output.rWaves(ii,:),'color',cRange(rheoSweep,:),'linewidth',2)
    end
elseif noRuns == 2
    tRheoSpikes = thTimes.RunTwo(rheoSweep.RunTwo,thTimes.RunTwo(rheoSweep.RunTwo,:)>-1);
    output.rWaves=findAPwaveforms(output.ccFile.RunTwo,tRheoSpikes,rheoSweep.RunTwo,samplerate); %aligned to threshold
    for ii = 1:size(output.rWaves,1)
        plot(1000/samplerate.*(1:size(output.rWaves,2)),output.rWaves(ii,:),'color',cRange.RunTwo(rheoSweep.RunTwo,:),'linewidth',2)
    end
end
setAx(gca)
xlabel('time (ms)')
title('overlaid spikes')

%phase plot
subplot(1,8,7:8)
hold on
for ii = 1:size(output.rWaves,1)
    output.dVdt(ii,:)=diff(output.rWaves(ii,:)).*10; %gives derivative in mV/ms
    if noRuns == 1
        plot(output.rWaves(ii,2:end),output.dVdt(ii,:),'color',cRange(rheoSweep,:),'linewidth',2)
    elseif noRuns == 2
        plot(output.rWaves(ii,2:end),output.dVdt(ii,:),'color',cRange.RunTwo(rheoSweep.RunTwo,:),'linewidth',2)
    end
end
xlabel('Vm (mV)')
ylabel('dV/dt (mV/ms)')
title('phase plot')
phaseAx = gca;
setAx(phaseAx);
phaseAx.YAxisLocation= 'origin'; phaseAx.XAxisLocation = 'origin';

%find ap peak data for rheobase sweep
if noRuns == 1
    [output.valRheoPeaks, output.tRheoPeaks, output.noRheoSpks] = getPeaks(output.ccFile,thVals,thTimes,rheoSweep);
elseif noRuns == 2
    [output.valRheoPeaks, output.tRheoPeaks, output.noRheoSpks] = getPeaks(output.ccFile.RunTwo,thVals.RunTwo,thTimes.RunTwo,rheoSweep.RunTwo);
end

%latency to first AP; defined as time from current injection to first AP peak
output.spkLat = (output.tRheoPeaks(1) - injstart)*(1000/samplerate); %gives first spike latency in ms

%threshold at rheobase
if noRuns == 1
    output.rheoThreshold = thVals(rheoSweep,1:output.noRheoSpks);
elseif noRuns == 2
    output.rheoThreshold = thVals.RunTwo(rheoSweep.RunTwo,1:output.noRheoSpks);
end

subplot(1,8,5:6)
if noRuns == 1
    for ii = 1:output.noRheoSpks
        scatter((find(output.rWaves(ii,:)==output.rheoThreshold(ii)))*(1000/samplerate),output.rheoThreshold(ii),50,cRange(end,:),'filled')
    end
elseif noRuns == 2
    for ii = 1:output.noRheoSpks
        scatter((find(output.rWaves(ii,:)==output.rheoThreshold(ii)))*(1000/samplerate),output.rheoThreshold(ii),50,cRange.RunTwo(end,:),'filled')
    end
end

subplot(1,8,7:8)
if noRuns == 1
    for ii = 1:output.noRheoSpks
        scatter(output.rheoThreshold(ii),output.dVdt(ii,find(output.rWaves(ii,:) == output.rheoThreshold(ii))-1),50,cRange(end,:),'filled')
    end
elseif noRuns == 2
    for ii = 1:output.noRheoSpks
        scatter(output.rheoThreshold(ii),output.dVdt(ii,find(output.rWaves(ii,:) == output.rheoThreshold(ii))-1),50,cRange.RunTwo(end,:),'filled')
    end
end

%spk amplitude at rheobase
if noRuns == 1
    output.rheoAmp = output.valRheoPeaks-thVals(rheoSweep,1:output.noRheoSpks);
elseif noRuns == 2
    output.rheoAmp = output.valRheoPeaks-thVals.RunTwo(rheoSweep.RunTwo,1:output.noRheoSpks);
end

%AHP
if noRuns == 1
    [output.AHP, output.tAHP] = getAHP(output.ccFile,output.tRheoPeaks,thTimes(rheoSweep,1:output.noRheoSpks),injend,rheoSweep,samplerate);
elseif noRuns == 2
    [output.AHP, output.idxAHP] = getAHP(output.ccFile.RunTwo,output.tRheoPeaks,thTimes.RunTwo(rheoSweep.RunTwo,1:output.noRheoSpks),injend,rheoSweep.RunTwo,samplerate);
end

output.rheoAHPlat=zeros(length(output.AHP),1);
output.rheoAHPval=zeros(length(output.AHP),1);
if noRuns == 1
    for ii = 1:length(output.AHP)
        output.rheoAHPlat(ii)=(output.idxAHP(ii)-thTimes(rheoSweep,ii))/(samplerate/1000);
        output.rheoAHPval(ii)=thVals(rheoSweep,ii)-output.AHP(ii);
    end
elseif noRuns == 2
    for ii = 1:length(output.AHP)
        output.rheoAHPlat(ii)=(output.idxAHP(ii)-thTimes.RunTwo(rheoSweep.RunTwo,ii))/(samplerate/1000);
        output.rheoAHPval(ii)=thVals.RunTwo(rheoSweep.RunTwo,ii)-output.AHP(ii);
    end
end

%halfwidth
if noRuns == 1
    output.halfwidth = getHalfwidth(output.ccFile,rheoSweep,samplerate,thTimes(rheoSweep,1:output.noRheoSpks),output.tRheoPeaks,output.rheoAmp);
elseif noRuns == 2
    output.halfwidth = getHalfwidth(output.ccFile.RunTwo,rheoSweep.RunTwo,samplerate,thTimes.RunTwo(rheoSweep.RunTwo,1:output.noRheoSpks),output.tRheoPeaks,output.rheoAmp);
end

%rheobase + 2sweeps data
%set up figure
rheop2Fig = figure(9);
rheop2Fig.Position = [130 1 675 133];

%spike train
hold on
if noRuns == 1
    plot(t,output.ccFile(:,1,rheoSweep+2),'linewidth',2,'color',cRange(rheoSweep+2,:))
    title({'rheobase + ' num2str(inputs.InjectiondeltapAEditField.Value*2) ' pA action potentials'})
elseif noRuns == 2
    plot(t,output.ccFile.RunTwo(:,1,rheoSweep.RunTwo+2),'linewidth',2,'color',cRange.RunTwo(rheoSweep.RunTwo+2,:))
    title({'rheobase + ' num2str(inputs.InjectiondeltapAEditField_2.Value*2) ' pA action potentials'})
end
setAx(gca)
xlim([injstart./samplerate-.05 injend./samplerate+.1])
xlabel('time (s)')
ylabel('Vm (mV)')

%find ap peak data for rheobase +2 sweep
if noRuns == 1
    [output.valRheop2Peaks, output.tRheop2Peaks, output.noRheop2Spks] = getPeaks(output.ccFile,thVals,thTimes,rheoSweep+2);
elseif noRuns == 2
    [output.valRheop2Peaks, output.tRheop2Peaks, output.noRheop2Spks] = getPeaks(output.ccFile.RunTwo,thVals.RunTwo,thTimes.RunTwo,rheoSweep.RunTwo+2);
end

%FR adaptation, ratio of first first ISI to last 3 ISI (<1 slowing down; >1 speeding up)
if noRuns == 1
    output.frRatio = getRatioFR(thTimes,rheoSweep+2,output.noRheop2Spks);
elseif noRuns == 2
    output.frRatio = getRatioFR(thTimes.RunTwo,rheoSweep.RunTwo+2,output.noRheop2Spks);
end

%amplitude adaptation, ratio of last three AP amplitudes to first AP amplitude (<1 getting smaller; >1 getting bigger)
if noRuns == 1
    output.rheop2Amp = output.valRheop2Peaks-thVals(rheoSweep+2,1:output.noRheop2Spks);
elseif noRuns == 2
    output.rheop2Amp = output.valRheop2Peaks-thVals.RunTwo(rheoSweep.RunTwo+2,1:output.noRheop2Spks);
end
output.ampRatio= mean([output.rheop2Amp(end) output.rheop2Amp(end-1) output.rheop2Amp(end-2)])/output.rheop2Amp(1);

%AP broadening, defined in Keshavarzi et al 2014 as the change in hw from AP1 to AP2
%redefined here as the ratio between the two (<1 gets faster; >1 gets slower)
if noRuns == 1
    output.halfwidthRP2 = getHalfwidth(output.ccFile,rheoSweep+2,samplerate,thTimes(rheoSweep+2,1:output.noRheop2Spks),output.tRheop2Peaks,output.rheop2Amp);
elseif noRuns == 2
    output.halfwidthRP2 = getHalfwidth(output.ccFile.RunTwo,rheoSweep.RunTwo+2,samplerate,thTimes.RunTwo(rheoSweep.RunTwo+2,1:output.noRheop2Spks),output.tRheop2Peaks,output.rheop2Amp);
end
output.broadeningRatio = output.halfwidthRP2(2)/output.halfwidthRP2(1);

%deltaAHP
if noRuns == 1
    [output.AHPRP2, output.tAHPRP2] = getAHP(output.ccFile,output.tRheop2Peaks,thTimes(rheoSweep+2,1:output.noRheop2Spks),injend,rheoSweep+2,samplerate);
elseif noRuns == 2
    [output.AHPRP2, output.idxAHPRP2] = getAHP(output.ccFile.RunTwo,output.tRheop2Peaks,thTimes.RunTwo(rheoSweep.RunTwo+2,1:output.noRheop2Spks),injend,rheoSweep.RunTwo+2,samplerate);
end
output.dAHP = output.AHPRP2(end)-output.AHPRP2(1);

output
end

function [thVals,thTimes] = findSpkThresh(file,injectionStart,injectionEnd,fs,dI,pColor)

%init vars
noSweeps=length(dI);
thVals=[];
thTimes=[];
dVm = diff(file); %takes approx derivative (lose element n = 1)
dVm(:,2,:) = [];
dVm = dVm.*(fs/1000); %converts to V/s
file(1,:,:) = []; %remove element n = 1 from Vm data
t = fs^-1:fs^-1:size(file,1)/fs;
dVmThresh = 20; %manually set to 20V/s

for ii = 1:noSweeps
    tThresh=[];
    oThresh=[];
    if dI(ii) > 0
        %find spikes
        noSpkSweep=1;
        oThresh=find(dVm((.002*fs+injectionStart):injectionEnd,1,ii)>dVmThresh); %find points above the threshold of dVm during current injection (2ms buffer from start of step)
        for jj = 1:length(oThresh)
            if jj == 1
                %get point before and after threshold
                p0 = oThresh(jj)+injectionStart+.002*fs-1; %after
                p1 = oThresh(jj)+injectionStart+.002*fs-2; %before
                afterVal = abs(dVmThresh-dVm(p0,1,ii));
                beforeVal = abs(dVmThresh-dVm(p1,1,ii));
                if afterVal <= beforeVal %if value after threshold is closer to threshold...
                    tThresh(noSpkSweep) = p0;
                else %if value before threshold is closer to threshold...
                    tThresh(noSpkSweep) = p1;
                end
                noSpkSweep=noSpkSweep+1;
            else
                if oThresh(jj)>.001*fs+tThresh(noSpkSweep-1)-injectionStart %if time over threshold is at least 1ms past the last AP (refractory period)
                    p0 = oThresh(jj)+injectionStart+.002*fs-1; %after
                    p1 = oThresh(jj)+injectionStart+.002*fs-2; %before
                    afterVal = abs(dVmThresh-dVm(p0,1,ii));
                    beforeVal = abs(dVmThresh-dVm(p1,1,ii));
                    if afterVal <= beforeVal %if value after threshold is closer to threshold...
                        tThresh(noSpkSweep) = p0;
                    else %if value before threshold is closer to threshold...
                        tThresh(noSpkSweep) = p1;
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
        plot(t,file(:,1,ii),'color',pColor.*1.5)
        ylim([min(min(file(:,1,:)))-0.1*max(range(file(:,1,:))) max(max(file(:,1,:)))+0.1*max(range(file(:,1,:)))])
        ylabel('Membrane Potential (mV)')
        xlabel('Time (ms)')
        xlim([.5 1.75])
        if isempty(tThresh) ~= 1
            scatter(tThresh/fs,file(tThresh,1,ii),'markeredgecolor',pColor)
        end
        legend('Vm','threshold');
        thisLegend = legend;
        thisLegend.Box = 'off';
        thisAx = gca;
        setAx(thisAx);
        subplot(1,2,2)
        hold on
        plot(t,dVm(:,1,ii),'color',pColor.*1.5)
        line([t(1) t(end)],[dVmThresh dVmThresh],'color','k','linewidth',1.25);
        ylabel('dV/dt (V/sec)')
        xlabel('Time (ms)')
        xlim([.5 1.75])
        legend('dVm/dt','threshold');
        thisLegend = legend;
        thisLegend.Box = 'off';
        thisAx = gca;
        setAx(thisAx);
    end
    
    %Fill Matrix
    %sweep x spk number; -1s designate no spike
    if isempty(tThresh)~=1
        if length(tThresh)==size(thVals,2) %if noSpks(thisSweep) == maxSpks(1stSweep:thisSweep-1)
            thVals(ii,:)=file(tThresh,1,ii);
            thTimes(ii,:)=tThresh;
        elseif length(tThresh)<size(thVals,2) %if noSpks(thisSweep) < maxSpks(1stSweep:thisSweep-1)
            %adds -1 to complete row; script will ignore -1s
            noRows = size(thVals,2)-length(tThresh);
            addrows=-1.*ones(1,noRows);
            thVals(ii,:)=vertcat(file(tThresh,1,ii),addrows');
            tThresh=horzcat(tThresh,addrows);
            thTimes(ii,:)=tThresh;
        elseif length(tThresh)>size(thVals,2) %if noSpks(thisSweep) > maxSpks(1stSweep:thisSweep-1)
            %adds -1 to complete all other rows; script will ignore -1s
            noCol = length(tThresh)-size(thVals,2);
            addcol=-1.*ones(size(thVals,1),noCol);
            thVals=horzcat(thVals,addcol);
            thTimes=horzcat(thTimes,addcol);
            thVals(ii,:)=file(tThresh,1,ii);
            thTimes(ii,:)=tThresh;
        end
    elseif isempty(tThresh)==1
        thVals(ii,:)=-1;
        thTimes(ii,:)=-1;
    end
    clf
end
close(figure(1))

end

function Rm =getRm(file,blstart,blend,injstart,fs,rheobase,dI,cRange)
%uses all sweeps prior to rheobase, first 100ms of injection compared to 100ms pre injection
if rheobase > 3 %if you start with negative current this shouldn't be a problem
    
    %IV plot to get membrane resistance
    %get dI and dV
    dI(rheobase:end) = [];
    dV=-1.*ones(1,rheobase-1);
    for ii = 1:rheobase-1
        %take voltage from 50-150ms post injection, 50ms for membrane charging & go for 100ms before potential HCN activation (HCN would introduce non-linearities)
        dV(ii) = mean(file(injstart+0.05*fs:(injstart+0.15*fs-1),1,ii))-mean(file(blstart:blend,1,ii));
    end
    
    %get linear fit of IV plot
    IVfit=nlinfit(dI',dV,@dr_fitline,[0 1]);
    plotIVfit=IVfit(1).*dI+IVfit(2);
    Rm=1000*IVfit(1); %gives Rm in MOhms
    
    %Plot IV
    rmFig = figure(2);
    rmFig.Position = [1135 600 275 210];
    hold on
    scatter(dI,dV,50,cRange(1,:),'fill')
    plot(dI,plotIVfit,'linewidth',2,'color',cRange(end,:))
    ylim([-40 20])
    ylabel('Vm (mV)')
    xlabel('current (pA)')
    title({'IV plot; Rm = ' num2str(Rm) ' M\Omega'})
    ivAx = gca;
    setAx(ivAx);
    ivAx.XAxisLocation = 'origin'; ivAx.YAxisLocation = 'origin';
else
    Rm = NaN;
end
end

function [mtau] = getTau(data,sweep,injstart,samplerate,cRange)
%set up single exp fit
xvals=(injstart:injstart+.25*samplerate)'; %first 250ms of injection
endFallingV=find(data(xvals,1,sweep)==min(data(xvals,1,sweep)))+xvals(1); %when voltage reaches most negative Vm
yvals=(data(xvals(1):endFallingV(1),1,sweep));
yPrime = data(xvals,1,sweep);
newXVals=(xvals(1):endFallingV(1))';
dt=newXVals-newXVals(1);
dtPrime = xvals-xvals(1);
dV=yvals-yvals(end);
dVPrime = yPrime - yvals(end);
betaBest=nlinfit(dt,dV,@exp1fit,[dV(1) mean(dt)]);

mtau=betaBest(2)/10; %gets tau_m in milliseconds

%plot
tauFig=figure(3);
tauFig.Position = [1135 300 275 210];
hold on
plot(dtPrime,dVPrime,'color',cRange(sweep,:),'linewidth',2)
plot(dt,exp1fit(betaBest,dt),'color',cRange(end,:),'linewidth',2,'linestyle','-.')
title({'Membrane \tau = ' mtau ' ms'})
tauAx = gca;
setAx(tauAx);
tauAx.XTick = 0:.05*samplerate:.25*samplerate;
tauAx.XTickLabel = {'0' '50' '100' '150' '200' '250'};
xlabel('time (ms)')
ylabel('relative Vm (mV)')


end

function vsag = getSag(data,sweep,injstart,injend,blstart,blend,samplerate,cRange)
%returns voltage sag
% 100*(minV-steadystateV)/(minV-baselineV)
% minV is the mean value of the minumum point in the hyperpolarizing injection
% steadystateV is the mean value over the last 200ms of the hyperpolarizing injection
% baselineV is the mean value over the first 500ms of the sweep

t = 1/samplerate:1/samplerate:size(data,1)/samplerate;
minV=min(data(injstart:injend,1,sweep));
ssV=mean(data((injend-0.2*samplerate):injend,1,sweep));
blV=mean(data(blstart:blend,1,sweep));

vsag=100*((minV-ssV)/(minV-blV));

%plot
sagFig=figure(4);
sagFig.Position = [1135 1 275 210];
hold on
plot(t,data(:,1,sweep),'color',cRange(1,:),'linewidth',2)
line([injstart./samplerate-.05 injstart./samplerate-.05],[blV minV],'linewidth',2','color',cRange(round(.67*size(cRange,1)),:))
line([injstart./samplerate-.025 injstart./samplerate-.025],[blV ssV],'linewidth',2','color',cRange(end,:))
text(injstart./samplerate+.025, blV,'maximum voltage response','color',cRange(round(.67*size(cRange,1)),:),'FontWeight','bold')
text(injstart./samplerate+.025, blV-5,'steady state voltage response','color',cRange(end,:),'FontWeight','bold')
xlim([injstart./samplerate-.1 injend./samplerate+.15])
title({'Voltage sag = ' num2str(vsag) '%'})
sagAx = gca;
setAx(sagAx);
xlabel('time (s)')
ylabel('Vm (mV)')
end

function [peaks, times, number] = getRbndspks(data,sweep,t1,stop,cRange,samplerate)
%gets rebound spikes, define as voltage peaks above -15mV

t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

[allpeaks,alltimes]=findpeaks(data(t1:stop,1,sweep));

%plot
rbndFig = figure(5);
rbndFig.Position = [825 260 285 130];
hold on
plot(t,data(:,1,sweep),'color',cRange(sweep,:),'linewidth',2)
scatter(alltimes./samplerate,allpeaks,'markeredgecolor',cRange(end,:))
xlim([t1./samplerate t1./samplerate+1])
setAx(gca);
xlabel('time (s)')
ylabel('Vm (mv)')

peaks=allpeaks(allpeaks>-15);
times=alltimes(allpeaks>-15)-1+t1;
number=length(peaks);
title({num2str(number) ' rebound spikes'})

end

function [maxFiring] = getMaxFiring(data,noSweeps,rheoSweep,spkMatrix,samplerate,injStart,cRange)
t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

%find last current injection before AP reduction
nospks=zeros(noSweeps,1);
for ii = rheoSweep:noSweeps
    nospks(ii)=sum(spkMatrix(ii,:)~=-1);
end
dspks=zeros(length(rheoSweep+2:noSweeps),1);
for ii = rheoSweep+2:noSweeps
    dspks(ii)=nospks(ii)-nospks(ii-1);
end
negsweeps=find(dspks<0);
if isempty(negsweeps) ~= 1
    maxsweep=negsweeps(1)-1;
else
    maxsweep=noSweeps;
end

%get isi during first 200ms
spks=(spkMatrix(maxsweep,:)>=injStart)&(spkMatrix(maxsweep,:)<=(injStart+.2*samplerate));
tspks=spkMatrix(maxsweep,spks); %gets time for spk threshold

%get ISI in sec
for ii = 1:length(tspks)
    if ii ~=length(tspks)
        ISI(ii)=(tspks(ii+1)-tspks(ii))/samplerate;
    else
        if ii ~= length(spkMatrix(maxsweep,:))
            ISI(ii)=(spkMatrix(maxsweep,ii+1)-tspks(ii))/samplerate;
        end
    end
end
meanISI=mean(ISI);

%get FR in Hz
maxFiring=1/meanISI;

%plot
maxfrFig = figure(6);
maxfrFig.Position = [825 5 285 175];
hold on
plot(t,data(:,1,maxsweep),'linewidth',2,'color',cRange(maxsweep,:))
xlim([injStart./samplerate-.1 injStart./samplerate+.9])
title({'Max Firing Rate = ' num2str(maxFiring) ' Hz'})
xlabel('time (s)')
ylabel('Vm (mV)')
setAx(gca);

end
function [curve] = getFICurve(data,noSweeps,rheoSweep,sweepDelta,spkMatrix,samplerate,injStart,cRange)
%gets and plots current injection/firing rate relationship
curve = -1.*ones(noSweeps,2);
t = 1/samplerate:1/samplerate:size(data,1)/samplerate;

for ii = 1:noSweeps
    if ii < rheoSweep
        curve(ii,1) = 0;
        curve(ii,2) = sweepDelta(ii,2);
    elseif ii >= rheoSweep
        spks = []; tspks = [];
        %get isi during first 200ms
        spks=(spkMatrix(ii,:)>=injStart)&(spkMatrix(ii,:)<=(injStart+.2*samplerate));
        tspks=spkMatrix(ii,spks); %gets time for spk threshold
        
        %get isi in sec
        isi = [];
        for jj = 1:length(tspks)
            if jj ~= length(tspks)
                isi(jj) = (tspks(jj+1)-tspks(jj))/samplerate;
            else
                if jj ~= length(spkMatrix(ii,:))
                    if spkMatrix(ii,jj+1) ~= -1
                        isi(jj) = (spkMatrix(ii,jj+1)-tspks(jj))/samplerate;
                    end
                end
            end
        end
        if isempty(isi) ~= 1
            mISI = mean(isi);
            %get FR in Hz
            curve(ii,1) = mISI^(-1);
        else
            curve(ii,1) = 0;
        end
        curve(ii,2) = sweepDelta(ii,2);
    end
end

deleteThese = find(curve(:,2)<0);
curve(deleteThese,:) = [];


%plot
fiFig = figure(8);
fiFig.Position = [235 500 285 175];
hold on
plot(curve(:,2),curve(:,1),'linewidth',2,'color',cRange(1,:))
title({'FR plot'})
xlabel('current injection (pA)')
ylabel('Firing Rate (Hz)')
setAx(gca);

end

function [peakval, tpeak, nopeaks] = getPeaks(data,threshVals,threshTimes,sweep)
%detects action potential peaks

nopeaks=sum(threshVals(sweep,:)~=-1);
peakval=zeros(nopeaks,1);
tpeak=zeros(nopeaks,1);
if nopeaks~=0 %if APs in current sweep
    for ii=1:nopeaks
        if ii ~= nopeaks
            [peakval(ii),tpeak(ii)]=max(data(threshTimes(sweep,ii):threshTimes(sweep,ii+1),1,sweep)); %finds max point between threshold_n and threshold_n+1
        elseif ii == nopeaks
            [peakval(ii),tpeak(ii)]=max(data(threshTimes(sweep,ii):threshTimes(sweep,ii)+100,1,sweep)); %finds max point between threshold_n and threshold_n+1
        end
        tpeak(ii)=tpeak(ii)+threshTimes(sweep,ii)-1; %fixes time
    end
end
end

function [AHP, tAHP] = getAHP(data,peaktimes,threshtimes,injend,sweep,samplerate)
%detects AHP, defined as the largest hyperpolarization deflection within the first 100ms after AP
AHP=zeros(length(peaktimes),1);
tAHP=zeros(length(peaktimes),1);
for ii = 1:length(peaktimes)
    if ii ~=length(peaktimes)
        deltaspk=threshtimes(ii+1)-threshtimes(ii);
        if deltaspk <= 0.1*samplerate
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii+1),1,sweep));
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.1*samplerate,1,sweep));
        end
        tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
    else
        deltaend=injend-threshtimes(ii);
        if deltaend <= 0.1*samplerate
            if deltaend <=0
                AHP(ii)=[];
                tAHP(ii)=[];
            else
                if (injend-peaktimes(ii)) <= 0.02*samplerate
                    AHP(ii)=[];
                    tAHP(ii)=[];
                else
                    [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):injend,1,sweep));
                    tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
                end
            end
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.1*samplerate,1,sweep));
            tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
        end
    end
end
end

function [hw,halfamp1,tHalf1,halfamp2,tHalf2] = getHalfwidth(data,sweep,samplerate,threshtimes,peaktime,amp)
%defined as the width (in msec) of the AP waveform at the half-amplitude
halfamp1=zeros(length(threshtimes),1);
halfamp2=zeros(length(threshtimes),1);
hw=zeros(length(threshtimes),1);
tHalf1=zeros(length(threshtimes),1);
tHalf2=zeros(length(threshtimes),1);
for ii = 1:length(threshtimes)
    dfrom_halfamp1=[];
    dfrom_halfamp2=[]; %#ok<*NASGU>
    halfamp=amp(ii)/2;
    dfrom_halfamp1=abs(data(threshtimes(ii):peaktime(ii)-1,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist1]=min(dfrom_halfamp1);
    dfrom_halfamp2=abs(data(peaktime(ii):peaktime(ii)+0.0025*samplerate,1,sweep)-(data(peaktime(ii),1,sweep)-halfamp));
    [~,t_mindist2]=min(dfrom_halfamp2);
    tHalf1(ii)=threshtimes(ii)+t_mindist1-1;
    halfamp1(ii)=data(tHalf1(ii),1,sweep);
    tHalf2(ii)=peaktime(ii)+t_mindist2-1;
    halfamp2(ii)=data(tHalf2(ii),1,sweep);
    hw(ii)=(tHalf2(ii)-tHalf1(ii))/(samplerate/1000);
end
end

function [ratio] = getRatioFR(threshtimes,sweep,nospks)
%ratio of last 3 ISI to first ISI in sweep

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
