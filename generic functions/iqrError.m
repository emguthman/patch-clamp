function [lowerBound, upperBound] = iqrError(data,dim)
%function to determine the standard error of the mean for a data set
%created 05-03-18
lowerBound = median(data,dim) - quantile(data,.25,dim);
upperBound = quantile(data,.75,dim) - median(data,dim);
end