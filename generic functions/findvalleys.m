%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%$%%%%% Calculates Local Minima %%%%%%%
%%%%%%%%%% Created: 11-05-2015 %%%%%%%%%%
%%%%%%%%%%% Edited: 7-19-2016 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [min_val,loc] = findvalleys(trace)
trace_inv = -1.*trace;
[inv_val,inv_loc] = findpeaks(trace_inv);
min_val=-1*inv_val; loc=inv_loc;

