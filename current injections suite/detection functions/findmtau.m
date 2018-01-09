function [mtau,vdrop]=findmtau(data,sweep,injstart,injend,blstart,blend,samplerate)
%finds membrane tau. only calc if dV <=10mV

%gives voltage drop to sweep
vdrop=min(data(injstart:injend,1,sweep))-mean(data(blstart:blend,1,sweep));

%set up single exp fit
xvals=(injstart:injstart+.25*samplerate)'; %first 250ms of injection
end_fallingV=find(data(xvals,1,sweep)==min(data(xvals,1,sweep)))+xvals(1); %when voltage reaches most negative Vm
yvals=(data(xvals(1):end_fallingV(1),1,sweep));
new_xvals=(xvals(1):end_fallingV(1))';
dt=new_xvals-new_xvals(1);
dV_in=yvals-yvals(end);
beta_best=nlinfit(dt,dV_in,@emg_exp1fit,[dV_in(1) mean(dt)]);

%plot
tauFig=figure(99);
hold on
plot(dt,dV_in)
plot(dt,emg_exp1fit(beta_best,dt))
title('Single Exp Fit for falling voltage response to -100pA current step')

mtau=beta_best(2)/10; %gets tau_m in milliseconds

close(tauFig)

% if mtau > 100 %if mtau is over 100ms, it's likely not correct. Redo analysis at -125pA current injection and use this measure.
%     % Report as tau was calculated at -75pA if accurate exponential fit was unable to be calculated at -100pA current injection due to synaptic activity.
%     clear mtau vdrop xvals end_fallingV yvals new_xvals dt dV_in beta_best
%     
%     sweep = sweep -1;
%     
%     %gives voltage drop to sweep
%     vdrop=min(data(injstart:injend,1,sweep))-mean(data(blstart:blend,1,sweep));
%     
%     %set up single exp fit
%     xvals=(injstart:injstart+.25*samplerate)'; %first 250ms of injection
%     end_fallingV=find(data(xvals,1,sweep)==min(data(xvals,1,sweep)))+xvals(1); %when voltage reaches most negative Vm
%     yvals=(data(xvals(1):end_fallingV(1),1,sweep));
%     new_xvals=(xvals(1):end_fallingV(1))';
%     dt=new_xvals-new_xvals(1);
%     dV_in=yvals-yvals(end);
%     beta_best=nlinfit(dt,dV_in,@emg_exp1fit,[dV_in(1) mean(dt)]);
%     
%     %plot
%     tauFig=figure(99);
%     hold on
%     plot(dt,dV_in)
%     plot(dt,emg_exp1fit(beta_best,dt))
%     title('Single Exp Fit for falling voltage response to -100pA current step')
%     
%     mtau=beta_best(2)/10; %gets tau_m in milliseconds
% end

% close(tauFig)

end