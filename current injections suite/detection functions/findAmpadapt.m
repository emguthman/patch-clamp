function [ratio] = findAmpadapt(peakvals)
%ratio of last three AP amplitudes to first AP amplitude
ratio=mean([peakvals(end) peakvals(end-1) peakvals(end-2)])/peakvals(1);

end