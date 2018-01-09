%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biophys Properties from Current Injection %
%%%%%%%%%%% Created: 09-28-2015 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 07-09-2017 %%%%%%%%%%%%%
%%%%%%%%%%% Renamed on 11-23-15 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% INIT VARS
close all
clearvars

%initiate vars -- set these to what is appropriate
sample_rate = 10000; %10kHz corresponds to each data point being 0.1 ms long
injstart = 5469; %injection start
injend = 11468; %injection end
blstart = injstart-.1*sample_rate-1; %baseline holding current start
blend = injstart-1; %baseline holding current end

% Choose Color GUI
cc = inputdlg('Which Current Injections? small (-100 to +100, d10), large (-200 to +400, d25), both');
if strcmp(cc{1},'small')
    clear cc; cc = 0;
elseif strcmp(cc{1},'large')
    clear cc; cc = 1;
elseif strcmp(cc{1},'both')
    clear cc; cc = 2;
end

%% LOAD DATA
%Load files
ci_dir=uigetdir;
cd(ci_dir);

%data
contents = dir('*.abf');
filenames = {contents.name}';
ci_files = fullfile(cd,filenames);

%% -100 TO +100 D10 DATA
if cc ~= 1
    pCFile = abfload(ci_files{1},'sweeps','a');
    %creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
    %M is dependent on sample rate
    %first data point in pCFile is at t=0
    no_sweeps=size(pCFile,3); %gives number of sweeps per file
    
    %find dI
    sweep_dI=zeros(no_sweeps,2);
    sweep_dI(:,1)=1:1:no_sweeps;
    sweep_dI(:,2)=-100:10:100;
    [~,sweep_tau]=min(abs(sweep_dI(:,2)+50)); %finds sweep index closest to -50pA sweep -- for tau-m
    
    %find AP Threshold
    [thresh_vals_smallsteps,thresh_times_smallsteps]=find_apthresh(pCFile,blstart,blend,injstart,injend,sample_rate,sweep_dI(:,2));
    thresh_sweep=find(thresh_vals_smallsteps(:,1)~=0);
    if length(thresh_sweep) >1
        thresh_sweep(2:end)=[];
    elseif isempty(thresh_sweep)
        thresh_sweep=no_sweeps+1;
    end
    
    % Set up CS figure
    for ii = 1:no_sweeps
        spkFig=figure(1);
        spkFig.Position=[0 60 850 475];
        subplot(2,1,1)
        hold on
        
        %plot voltage response to current step
        plot(0:1/sample_rate:(size(pCFile,1)-1)/sample_rate,pCFile(:,1,ii))
        ylim([min(min(pCFile(:,1,:)))-0.1*max(range(pCFile(:,1,:))) max(max(pCFile(:,1,:)))+0.1*max(range(pCFile(:,1,:)))])
        set(gca,'xtick',0:.5:3);
        set(gca,'xticklabel',{'0','500','1000','1500','2000','2500','3000'}) %sets labels to milliseconds -- assumes 3 second long sweep
        ylabel('Membrane Potential (mV)','fontweight','bold')
        xlabel('Time (ms)','fontweight','bold')
        
        %plot current injection
        figure(1)
        subplot(2,1,2)
        hold on
        plot(0:1/sample_rate:(size(pCFile,1)-1)/sample_rate,pCFile(:,2,ii))
        ylim([min(min(pCFile(:,2,:)))-0.1*max(range(pCFile(:,2,:))) max(max(pCFile(:,2,:)))+0.1*max(range(pCFile(:,2,:)))])
        set(gca,'xtick',0:.5:3);
        set(gca,'xticklabel',{'0','500','1000','1500','2000','2500','3000'}) %sets labels to milliseconds -- assumes 3 second long sweep
        ylabel('Holding Current (pA)','fontweight','bold')
        xlabel('Time (ms)','fontweight','bold')
        text(1.75,min(min(pCFile(:,2,:)))+10,['Current Step: ' num2str(sweep_dI(ii,2)) ' pA'],'fontweight','bold')
        clf
    end
    close(spkFig)
    
    %% Get membrane resistance
    [Rm.value,Rm.dp]=MemRes(pCFile,blstart,blend,injend,sample_rate,thresh_sweep);
    
    if cc == 0
        %% get tau-m
        tau_sweep = find(sweep_dI(:,2) == -100); %take tau at -100pA sweep;
        [mtau,~]=findmtau(pCFile,tau_sweep,injstart,injend,blstart,blend);
    end
    smallData = pCFile;
end

%% -200 TO +400 D25 DATA
if cc ~= 0
    pCFile=[];
    if cc == 1
        pCFile = abfload(ci_files{1},'sweeps','a');
    elseif cc == 2
        pCFile = abfload(ci_files{2},'sweeps','a');
    end
    %creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
    %M is dependent on sample rate
    no_sweeps=size(pCFile,3); %gives number of sweeps per file
    
    %find dI
    sweep_dI=zeros(no_sweeps,2);
    sweep_dI(:,1)=1:1:no_sweeps;
    sweep_dI(:,2)=-200:25:400;
    
    %find AP Threshold
    % thresh_sweep=[];
    [thresh_vals_largesteps,thresh_times_largesteps]=find_apthresh(pCFile,blstart,blend,injstart,injend,sample_rate,sweep_dI(:,2));
    
    %fix threshold times
    thresh_sweep=find(thresh_vals_largesteps(:,1)~=0);
    if length(thresh_sweep) >1
        thresh_sweep(2:end)=[];
    end
    no_sweepswspks=no_sweeps-thresh_sweep-1; %no of sweeps with spikes.
    
    %% Get Voltage Sag
    % Vsag is calculated as the percent voltage sag for the -200pA step
    % 100*(minV-steadystateV)/(minV-baselineV)
    % minV is the mean value of the minumum point in the hyperpolarizing injection
    % steadystateV is the mean value over the last 200ms of the hyperpolarizing injection
    % baselineV is the mean value over the first 500ms of the sweep
    
    %get value
    vsag=findvsag(pCFile,1,injstart,injend,blstart,blend,sample_rate);
    
    %% Get rebound spikes
    [rspk_vals,rspk_times,no_rspks]=find_rspks(pCFile,1,injend+1,injend+.5*sample_rate);
    
    %% Display sag and rebound spikes (if any)
    vsagFig=figure(1);
    vsagFig.Position=[0 655 675 150];
    hold on
    plot(pCFile(:,1,1))
    scatter(rspk_times,rspk_vals);
    text(2*sample_rate,-55,['Voltage Sag: ' num2str(vsag) '%'],'fontweight','bold')
    text(2*sample_rate,-60,['No. Rebound Spks: ' num2str(no_rspks)],'fontweight','bold')
    set(gca,'xtick',0:.5*sample_rate:3*sample_rate);
    set(gca,'xticklabel',{'0','500','1000','1500','2000','2500','3000'}) %sets labels to milliseconds -- assumes 3 second long sweep
    ylabel('Membrane Potential (mV)','fontweight','bold')
    xlabel('Time (ms)','fontweight','bold')
    
    %% Get MaxFR for cell
    %Defined as inverse of ISI during first 200ms of last current injection before reduction in AP firing occurs
    
    %find last current injection before AP reduction
    nospks=zeros(no_sweeps,1);
    for ii = thresh_sweep:no_sweeps
        nospks(ii)=sum(thresh_vals_largesteps(ii,:)~=0);
    end
    dspks=zeros(length(thresh_sweep+2:no_sweeps),1);
    for ii = thresh_sweep+2:no_sweeps
        dspks(ii)=nospks(ii)-nospks(ii-1);
    end
    negsweeps=find(dspks<0);
    if isempty(negsweeps) ~= 1
        maxsweep=negsweeps(1)-1;
    else
        maxsweep=no_sweeps;
    end
    
    %get maxFR
    maxFR=findFR(thresh_times_largesteps,maxsweep,injstart,injstart+.2*sample_rate,sample_rate);
    
    %plot figure of maxFR sweep
    maxFRfig=figure(3);
    maxFRfig.Position=[0 205 675 150];
    hold on
    t_ms=1000/sample_rate:1000/sample_rate:3000;
    plot(t_ms,pCFile(:,1,maxsweep))
    text(2000,-40,['Max FR: ' num2str(maxFR) 'Hz'],'fontweight','bold')
    ylim([min(min(pCFile(:,1,maxsweep)))-0.1*max(range(pCFile(:,1,maxsweep))) max(max(pCFile(:,1,maxsweep)))+0.1*max(range(pCFile(:,1,maxsweep)))])
    ylabel('Membrane Potential (mV)','fontweight','bold')
    xlabel('Time (ms)','fontweight','bold')
    
    %% Get rheobase data -- AP waveform data
    rheo=sweep_dI(thresh_sweep,2);
    
    %plot rheobase sweep
    rFig=figure(4);
    rFig.Position=[675 590 765 215];
    hold on
    plot(pCFile(:,1,thresh_sweep))
    text(20000,-40,['Rheobase: ' num2str(rheo) ' pA'],'fontweight','bold')
    ylim([min(min(pCFile(:,1,thresh_sweep)))-0.1*max(range(pCFile(:,1,thresh_sweep))) max(max(pCFile(:,1,thresh_sweep)))+0.1*max(range(pCFile(:,1,thresh_sweep)))])
    ylabel('Membrane Potential (mV)','fontweight','bold')
    xlabel('Time (ms)','fontweight','bold')
    
    %% Get ADF (calculated as in Cea-del Rio et al 2010)
    rADF = findADF(pCFile,blstart,blend,injend,thresh_sweep,sample_rate);
    text(20000,-30,['ADF: ' num2str(rADF) ' mV'],'fontweight','bold')
    
    %% Get AP latency
    %ap latency as defined as time from current step init to 1st ap peak at rheobase
    [rSpk_peaks,t_rSpks,n_rSpks]=findspkpeaks(pCFile,thresh_vals_largesteps,thresh_times_largesteps,thresh_sweep);
    rheo_APlat=(t_rSpks(1)-injstart)/(sample_rate/1000); %gives ap latency in msec
    text(20000,-20,['AP latency: ' num2str(rheo_APlat) ' msec'],'fontweight','bold')
    
    %% Get AP threshold
    rThresh=thresh_vals_largesteps(thresh_sweep,1:n_rSpks);
    text(20000,-10,['AP threshold: ' num2str(mean(rThresh)) ' mV'],'fontweight','bold')
    
    %% Get AP amplitude
    rAmp=rSpk_peaks-thresh_vals_largesteps(thresh_sweep,1:n_rSpks)';
    text(20000,0,['AP Amplitude: ' num2str(mean(rAmp)) ' mV'],'fontweight','bold')
    
    %plot thresh and peaks
    rFig;
    scatter(thresh_times_largesteps(thresh_sweep,1:n_rSpks),thresh_vals_largesteps(thresh_sweep,1:n_rSpks))
    scatter(t_rSpks,rSpk_peaks);
    
    %% Get AHP
    [rAHP,rawAHPtime]=findAHP(pCFile,t_rSpks,thresh_times_largesteps(thresh_sweep,1:n_rSpks),injend,thresh_sweep,sample_rate);
    rAHPlat=zeros(length(rAHP),1);
    rAHPval=zeros(length(rAHP),1);
    for ii = 1:length(rAHP)
        rAHPlat(ii)=(rawAHPtime(ii)-thresh_times_largesteps(thresh_sweep,ii))/(sample_rate/1000);
        rAHPval(ii)=thresh_vals_largesteps(thresh_sweep,ii)-rAHP(ii);
    end
    
    %plot AHP
    rFig;
    scatter(rawAHPtime,rAHP)
    text(20000,10,['AHP: ' num2str(mean(rAHPval)) ' mV, ' num2str(mean(rAHPlat)) ' msec'],'fontweight','bold')
    
    %% Get AP halfwidth
    [hw,halfamp1,t_half1,halfamp2,t_half2] = findAPhalfwidth(pCFile,thresh_sweep,sample_rate,thresh_times_largesteps(thresh_sweep,1:n_rSpks),t_rSpks,rAmp);
    
    %plot
    rFig;
    scatter(t_half1,halfamp1);
    scatter(t_half2,halfamp2);
    text(20000,20,['AP Halfwidth: ' num2str(mean(hw)) ' msec'],'fontweight','bold')
    
    %% Plot AP waveform from rheobase
    apFig=figure(5);
    apFig.Position=[675 60 250 200];
    hold on
    rWaves=findAPwaveforms(pCFile,t_rSpks,thresh_sweep,sample_rate);
    for ii = 1:size(rWaves,1)
        plot(rWaves(ii,:))
    end
    
    %% Get data from rheobase + 50pA sweep -- broadening, AP accomodation, and FR accomodation
    r50sweep=thresh_sweep+2;
    %plot rheobase +50pA sweep
    r50Fig=figure(6);
    r50Fig.Position=[675 310 765 215];
    hold on
    plot(pCFile(:,1,r50sweep))
    text(20000,-40,['Current Injection: ' num2str(rheo+50) 'pA'],'fontweight','bold')
    ylim([min(min(pCFile(:,1,r50sweep)))-0.1*max(range(pCFile(:,1,r50sweep))) max(max(pCFile(:,1,r50sweep)))+0.1*max(range(pCFile(:,1,r50sweep)))])
    ylabel('Membrane Potential (mV)','fontweight','bold')
    xlabel('Time (ms)','fontweight','bold')
    %find peaks
    [spk_peaks,t_spk,r50spks]=findspkpeaks(pCFile,thresh_vals_largesteps,thresh_times_largesteps,r50sweep);
    %plot peaks and amplitudes
    r50Fig;
    scatter(thresh_times_largesteps(r50sweep,1:r50spks),thresh_vals_largesteps(r50sweep,1:r50spks))
    scatter(t_spk,spk_peaks);
    
    %% Get FR adaptation ratio
    arFR=findFRadapt(thresh_times_largesteps,r50sweep,r50spks);
    text(20000,0,['FR adaptation ratio: ' num2str(arFR)],'fontweight','bold')
    
    %% Get Amp adaptation ratio
    r50Amp=spk_peaks-thresh_vals_largesteps(r50sweep,1:r50spks)';
    arAmp=findAmpadapt(r50Amp);
    text(20000,-10,['AP adaptation ratio: ' num2str(arAmp)],'fontweight','bold')
    
    %% Get AP broadening ratio
    %defined in Keshavarzi et al 2014 as the change in hw from AP1 to AP2
    %redefined here as the ratio between the two
    r50Amp=spk_peaks-thresh_vals_largesteps(r50sweep,1:r50spks)';
    [r50hw,~,~,~,~] = findAPhalfwidth(pCFile,r50sweep,sample_rate,thresh_times_largesteps(r50sweep,1:r50spks),t_spk,r50Amp);
    brAP=r50hw(2)/r50hw(1);
    text(20000,-20,['AP broadening ratio: ' num2str(brAP)],'fontweight','bold')
    
    %% Get deltaAHP
    [r50AHP,r50AHPt]=findAHP(pCFile,t_spk,thresh_times_largesteps(r50sweep,1:r50spks),injend,r50sweep,sample_rate);
    dAHP=finddeltAHP(r50AHP);
    
    %plot
    r50Fig;
    scatter(r50AHPt(1),r50AHP(1))
    scatter(r50AHPt(end),r50AHP(end))
    text(20000,-30,['\DeltaAHP: ' num2str(dAHP) ' mV'],'fontweight','bold')
    
    %% Plot AP waveform from rheobase
    ap50Fig=figure(7);
    ap50Fig.Position=[900 60 250 200];
    hold on
    r50Waves=findAPwaveforms(pCFile,t_spk,r50sweep,sample_rate);
    for ii = 1:size(r50Waves,1)
        plot(r50Waves(ii,:))
    end
    
    if cc == 1
        %% Get membrane resistance
        [Rm.value,Rm.dp]=MemRes(pCFile,blstart,blend,injend,sample_rate,thresh_sweep);
    end
    
    %% get tau-m
    tau_sweep = find(sweep_dI(:,2) == -100); %take tau at -100pA sweep;
    [mtau,~]=findmtau(pCFile,tau_sweep,injstart,injend,blstart,blend,sample_rate);
    if mtau > 100 %if mtau is over 100ms, it's likely not correct. Redo analysis with small current injections and use that data.
        mtau = [];
        [mtau,~]=findmtau(smallData,tau_sweep,injstart,injend,blstart,blend,sample_rate);
    end
    figure(1);
    text(2*sample_rate,-65,['R_{m}: ' num2str(Rm.value)],'fontweight','bold')
    text(2*sample_rate,-70,['\tau_{m}: ' num2str(mtau)],'fontweight','bold')
end
%% Save data
button = questdlg('Save data');
if strcmp(button,'Yes')
    [save_file,path_save]=uiputfile;
    full_save_file=fullfile(path_save,save_file);
    save(full_save_file);
end