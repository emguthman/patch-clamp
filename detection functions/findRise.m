function [rt, meanRT, tTwentyPercent] = findRise(data,tStart,tPeak,samplerate)
%returns 20-80% risetime and tTwentyPercent for use as latwentycy.
%requires data martix in time x trial format and for data to be baseline subtracted
%create 04-10-18
%last edit 04-26-18

%inits
rt = -1.*ones(size(data,2),1);
tTwentyPercent = -1.*ones(size(data,2),1);
for ii = 1:size(data,2)
    %inits
    eightyVal = [];
    eightyPts = [];
    eightyDist = [];
    eightyPt = [];
    twentyVal = [];
    twentyPts = [];
    twentyDist = [];
    twentyPt = [];
    peak = [];
    peakval = [];
    thisPeak = [];
    thisData = [];
    
    
    %find current peak
    peakval = abs(data(tPeak,ii));
    if length(tPeak) == 1
        thisPeak = tPeak;
    elseif length(tPeak) == size(data,2)
        thisPeak = tPeak(ii);
    end
    
    %rescale data
    thisData=data(tStart:thisPeak,ii);
    
    %get 20% and 80% points
    % 80% point
    eightyVal = 0.8*peakval;
    thesePts = find(abs(thisData)>eightyVal);
    if thesePts(1) == 1 %if this is the case, there's likely a spontaneous event potentially obscuring the rise
        %look for 80 pt from halfway from start to peak
        %20 pt section has a second check for spon event and will return Nan if the spon event prevents calculating rise
        clear thesePts
        thisStart = round(length(thisData)/2);
        thesePts = find(abs(thisData(thisStart:end))>eightyVal);
        thesePts = thesePts + thisStart -1;
        clear thisStart
    end
    if thesePts(1) == 1 %if still the case, likely spontaneous event, so return Nan
        rt(ii) = NaN;
        tTwentyPercent(ii) = NaN;
    else %continue
        eightyPts = [thesePts(1)-1; thesePts(1)];
        eightyDist = abs(abs(thisData(eightyPts))-eightyVal);
        [~,thisPt] = min(eightyDist);
        eightyPt = eightyPts(thisPt);
        thesePts = []; thisPt = [];
        
        % 20% point
        twentyVal = 0.2*peakval;
        thesePts = find(abs(thisData(1:eightyPt-1))<twentyVal);
        if isempty(thesePts) ~= 1
            twentyPts = [thesePts(end) thesePts(end)+1];
            twentyDist = abs(abs(thisData(twentyPts))-twentyVal);
            [~,thisPt] = min(twentyDist);
            twentyPt = twentyPts(thisPt);
            tTwentyPercent(ii) = twentyPt;
            %calculate 20-80 risetime in milliseconds
            rt(ii) = (eightyPt-twentyPt)*(1000/samplerate);
        else %if no data below 20% pt, spontaneous activity is obscuring rise, do not use for rise or latwentycy calculation
            rt(ii) = NaN;
            tTwentyPercent(ii) = NaN;
        end
        thesePts = []; thisPt = [];
    end
    
    
    %plot to validate
    riseFig=figure;
    riseFig.Position=[2555 300 800 500];
    hold on
    plot(thisData)
    scatter(eightyPt,thisData(eightyPt))
    scatter(twentyPt,thisData(twentyPt))
    close(riseFig);
    
    
end
meanRT = mean(rt(~isnan(rt)));
end
