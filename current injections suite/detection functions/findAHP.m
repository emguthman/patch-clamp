function [AHP, tAHP] = findAHP(data,peaktimes,threshtimes,injend,sweep,sample_rate)
%detects AHP, defined as the largest hyperpolarization deflection within the first 100ms after AP
AHP=zeros(length(peaktimes),1);
tAHP=zeros(length(peaktimes),1);
for ii = 1:length(peaktimes)
    if ii ~=length(peaktimes)
        deltaspk=threshtimes(ii+1)-threshtimes(ii);
        if deltaspk <= 0.1*sample_rate
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii+1),1,sweep));
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.1*sample_rate,1,sweep));
        end
        tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
    else
        deltaend=injend-threshtimes(ii);
        if deltaend <= 0.1*sample_rate
            if deltaend <=0
                AHP(ii)=[];
                tAHP(ii)=[];
            else
                if (injend-peaktimes(ii)) <= 0.02*sample_rate
                    AHP(ii)=[];
                    tAHP(ii)=[];
                else
                    [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):injend,1,sweep));
                    tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
                end
            end
        else
            [AHP(ii), tAHP(ii)]=min(data(peaktimes(ii):threshtimes(ii)+0.1*sample_rate,1,sweep));
            tAHP(ii)=tAHP(ii)-1+peaktimes(ii); %corrects AHP time
        end
    end
end
end