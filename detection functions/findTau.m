function [tau] = findTau(data,tStart,samplerate,SoD)
%returns tau as single or double exponential fit for a given current.
%create 04-10-18
%last edit 10-18-18

%find current peak
[~,peak] = max(abs(data(tStart:end)));
tPeak = peak - 1 + tStart;

%rescale data
thisData=data(tPeak:end);

%get decay
xvals=(1:length(thisData))';
dt=xvals-xvals(1);
dI=thisData-thisData(end);

%plot to decide if EPSC/IPSC single/double exp
decayFig=figure;
decayFig.Position=[50 600 400 200];
plot(data)

%get tau
if nargin == 3
    SoD = inputdlg('Tau from Single or Double Exponential?');
    if strcmp(SoD,'Single')
        beta_best=nlinfit(dt,dI,@exp1fit,[dI(1) mean(dt)]);
        tau=beta_best(2)*(1000/samplerate); %gets tau_m in milliseconds
        
        %plot on figure
        clf
        hold on
        plot(dI)
        plot(dt,exp1fit(beta_best,dt))
        close(decayFig);
    elseif strcmp(SoD,'Double')
        beta_best=nlinfit(dt,dI,@exp2fit,[dI(1) mean(dt) dI(1) mean(dt)]); %beta(1) == A1, beta(3) == A2, beta(2) == tau-d1, beta(4) == tau-d2
        %weighted decay constant (based on Li & Huntsman 2014, Neuroscience)
        wdecay = ((beta_best(1)*beta_best(2))+(beta_best(3)*beta_best(4)))/(beta_best(1)+beta_best(3)); %gives weighted decay constant in data points
        tau = wdecay*(1000/samplerate); %gives weighted decay in msec
        
        %plot on figure
        clf
        hold on
        plot(dI)
        plot(dt,exp2fit(beta_best,dt))
        close(decayFig);
    end
elseif nargin == 4
    if strcmp(SoD,'Single')
        beta_best=nlinfit(dt,dI,@exp1fit,[dI(1) mean(dt)]);
        tau=beta_best(2)*(1000/samplerate); %gets tau_m in milliseconds
        
        %plot on figure
        clf
        hold on
        plot(dI)
        plot(dt,exp1fit(beta_best,dt))
        close(decayFig);
    elseif strcmp(SoD,'Double')
        beta_best=nlinfit(dt,dI,@exp2fit,[dI(1) mean(dt) dI(1) mean(dt)]); %beta(1) == A1, beta(3) == A2, beta(2) == tau-d1, beta(4) == tau-d2
        %weighted decay constant (based on Li & Huntsman 2014, Neuroscience)
        wdecay = ((beta_best(1)*beta_best(2))+(beta_best(3)*beta_best(4)))/(beta_best(1)+beta_best(3)); %gives weighted decay constant in data points
        tau = wdecay*(1000/samplerate); %gives weighted decay in msec
        
        %plot on figure
        clf
        hold on
        plot(dI)
        plot(dt,exp2fit(beta_best,dt))
        close(decayFig);
    end
end 
end



