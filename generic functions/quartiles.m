function [lowerBound, upperBound] = quartiles(data)
%function to determine the standard error of the mean for a data set
%created 05-03-18
lowerBound = quantile(data,.25);
upperBound = quantile(data,.75);
end