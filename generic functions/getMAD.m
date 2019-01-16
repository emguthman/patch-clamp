function [mad] = getMAD(data)
%used to determine median absolute deviation of electrical noise for detecting synaptic events
%in my opinion, better than standard dev bc median abs dev allows for less influence of large deflections (ie spontaneous events) in determining the deviation of the background
mad = median(abs(data-median(data)));

end