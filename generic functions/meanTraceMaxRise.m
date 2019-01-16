%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% mean current, aligned by max slope %%%%%
%%%%%%%%%%% Created: 03-27-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 04-18-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alignedTraces,meanTrace,tMaxSlope, tMaxSlopeMeanTrace] = meanTraceMaxRise(data,tPeak,samplerate,windowDur,direction) %requires data in t x trial format
%align by max slope, five data point approximation using the 5pt midpoint derivative approximation method to determine slope
%ie get slope by 5pt deriv. method at t-2 to t+2, and get mean slope to get slope at point t
slope0 = -1.*ones(size(data,1)-6,size(data,2));
slope1 = -1.*ones(size(data,1)-6,size(data,2));
slope2 = -1.*ones(size(data,1)-6,size(data,2));
slope3 = -1.*ones(size(data,1)-6,size(data,2));
slope4 = -1.*ones(size(data,1)-6,size(data,2));
slope = -1.*ones(size(data,1)-6,size(data,2));

%five point midpoint derivative approx. (f'(x) = 1/12*dt * (f(x-2*dt) - 8*f(x-dt) + 8f(x+dt) - f(x+2dt))
for ii = 1:size(slope,2)
    for jj = 5:size(data,1)-4
        slope0(jj-2,ii) = 1/12*(data(jj-2,ii) - 8*data(jj-1,ii) + 8*data(jj+1,ii) - data(jj+2,ii)); %indexing is off by 2 dp relative to original data
        slope1(jj-2,ii) = 1/12*(data(jj-2+1,ii) - 8*data(jj-1+1,ii) + 8*data(jj+1+1,ii) - data(jj+2+1,ii));
        slope2(jj-2,ii) = 1/12*(data(jj-2+2,ii) - 8*data(jj-1+2,ii) + 8*data(jj+1+2,ii) - data(jj+2+2,ii));
        slope3(jj-2,ii) = 1/12*(data(jj-2-1,ii) - 8*data(jj-1-1,ii) + 8*data(jj+1-1,ii) - data(jj+2-1,ii));
        slope4(jj-2,ii) = 1/12*(data(jj-2-2,ii) - 8*data(jj-1-2,ii) + 8*data(jj+1-2,ii) - data(jj+2-2,ii));
    end
    slope(:,ii) = (slope0(:,ii) + slope1(:,ii) + slope2(:,ii) + slope3(:,ii) + slope4(:,ii))./5;
end


%pick how far back from peak to look for max slope
thisTime = 1000.*(1/samplerate:1/samplerate:size(data,1)/samplerate);
stimArtefactFig = figure;
plot(thisTime,data)
line([mean(tPeak)*(1000/samplerate) mean(tPeak)*(1000/samplerate)],[-20 20],'color','k','linewidth',2)
xlim([mean(tPeak)*(1000/samplerate)-25 mean(tPeak)*(1000/samplerate)+10]) %15ms window around stimulus
artBuffer = inputdlg('window time (ms) for finding peak slope (back from tPeak)')
artBuffer = (samplerate/1000)*str2double(artBuffer);
close(stimArtefactFig)

%find most negative slope
tMaxSlope = -1.*ones(size(slope,2),1);
maxSlope = -1.*ones(size(slope,2),1);
for ii = 1:size(slope,2)
    if length(tPeak) > 1
        thisPeak = tPeak(ii);
    elseif length(tPeak) == 1
        thisPeak = tPeak;
    end
    if direction == -1
        if (thisPeak-artBuffer-2) < 1
            [maxSlope(ii),tMaxSlope(ii)] = min(slope(1:thisPeak-2,ii));
        else
            [maxSlope(ii),tMaxSlope(ii)] = min(slope(thisPeak-artBuffer-2:thisPeak-2,ii)); %find min between 4ms pre peak and peak
            tMaxSlope(ii) = tMaxSlope(ii) + thisPeak - artBuffer -3 ; %correct to real time (in datapoints): 4ms pre-peak - 1
        end
    elseif direction == 1
        if (thisPeak-artBuffer-2) < 1
            [maxSlope(ii),tMaxSlope(ii)] = max(slope(1:thisPeak-2,ii));
        else
            [maxSlope(ii),tMaxSlope(ii)] = max(slope(thisPeak-artBuffer-2:thisPeak-2,ii));
            tMaxSlope(ii) = tMaxSlope(ii) + thisPeak - artBuffer -3 ; %correct to real time (in datapoints): 10ms pre-peak - 1
        end
    end
end

for ii = 1:size(slope,2)
    if length(tPeak) > 1
        thisPeak = tPeak(ii);
    elseif length(tPeak) == 1
        thisPeak = tPeak;
    end
%     %plot to check
%     tempFig = figure(99);
%     tempFig.Position;
%     subplot(2,1,1)
%     plot(data(1:thisPeak,ii))
%     if direction == -1
%         xlim([thisPeak-artBuffer-2 thisPeak])
%     elseif direction == 1
%         xlim([thisPeak-artBuffer-2 thisPeak])
%     end
%     
%     subplot(2,1,2)
%     plot(slope(1:thisPeak,ii))
%     if direction == -1
%         xlim([thisPeak-artBuffer-2 thisPeak])
%     elseif direction == 1
%         xlim([thisPeak-artBuffer-2 thisPeak])
%     end
end
% close (figure(99))

%align based on max slope
if direction == -1
    if (thisPeak-artBuffer-2) < 1
        alignedTraces = -1.*ones(size(slope));
        disp('trace too small for alignment//only use tMaxRise outputs for the data')
    else
        alignedTraces = -1.*ones(windowDur*samplerate+1,size(data,2));
        for ii = 1:size(data,2)
            if (min(tMaxSlope)-.035*samplerate) > 0
                alignedTraces(:,ii) = data(tMaxSlope(ii)-.035*samplerate:tMaxSlope(ii)+windowDur*samplerate-.035*samplerate,ii);
            else
                alignedTraces(:,ii) = data(tMaxSlope(ii)-artBuffer:tMaxSlope(ii)+windowDur*samplerate-artBuffer,ii);
            end
        end
    end
elseif direction == 1
    if (thisPeak-artBuffer-2) < 1
        alignedTraces = -1.*ones(size(slope));
        disp('trace too small for alignment//only use tMaxRise outputs for the data')
    else
        alignedTraces = -1.*ones(windowDur*samplerate+1,size(data,2));
        for ii = 1:size(data,2)
            if (min(tMaxSlope)-.035*samplerate) > 0
                alignedTraces(:,ii) = data(tMaxSlope(ii)-.035*samplerate:tMaxSlope(ii)+windowDur*samplerate-.035*samplerate,ii);
            else
                alignedTraces(:,ii) = data(tMaxSlope(ii)-artBuffer:tMaxSlope(ii)+windowDur*samplerate-artBuffer,ii);
            end
        end
    end
end

%get mean
meanTrace = mean(alignedTraces,2);
if direction == -1
    if (tMaxSlope(ii)-.035*samplerate) > 0
        tMaxSlopeMeanTrace = .035*samplerate +1;
    else
        tMaxSlopeMeanTrace = artBuffer +1;
    end
elseif direction == 1
    if (tMaxSlope(ii)-.035*samplerate) > 0
        tMaxSlopeMeanTrace = .035*samplerate +1;
    else
        tMaxSlopeMeanTrace = artBuffer +1;
    end
end

end