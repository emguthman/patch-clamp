function [dI] = exp2fit(beta,dt)
%functino for double exponential decay, for nlinfit, created 08-13-2016,
%modified 08-13-2016

dI=beta(1)*exp(-dt/beta(2))+beta(3)*exp(-dt/beta(4));
end