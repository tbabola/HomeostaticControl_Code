function [ROIstats, SGNdFoF] = SGN_ROIsignals(SGNsingles, baseline, drugs1, drugs2)
    switch nargin
        case 2
            drugs1 = 0;
            drugs2 = 0;
        case 3
            drugs2 = 0;
        case 1
            baseline = int16([1 300]);
            drugs1 = 0;
            drugs2 = 0;
    end
    
    numCellsForEvent = 3;
    
    %correct baseline
    Fo = [];
    SGNdFoF = [];
    [t,n] = size(SGNsingles);
     for i = 1:n
         meanSGN = mean(SGNsingles(:,i));
         SGNsingles(:,i) = msbackadj([1:t]',SGNsingles(:,i),'WindowSize',100,'StepSize',100) + meanSGN;%     
     end   

%plot signals
    offset = [];
    for i=1:n
        if max(SGNsingles(:,1)) > 256
            mult = 1000;
        else
            mult = 10;
        end
    offset(:,i) = SGNsingles(:,i) - (i-1)*mult;
    end
    figure; plot(offset);

%     offset = [];
%     for i=1:n
%     offset(:,i) = dFoF(:,i) - (i-1)*1;
%     end
%     figure; plot(offset);

    SGNbl = SGNsingles;
    SGNmean = mean(SGNbl,2);
    
    
    Fo = prctile(SGNbl,10); %bottom tenth percentile for baseline
    SGNdFoF = (SGNbl - Fo)./(repmat(Fo,t,1));

    %plot signals
    offset = [];
    for i=1:n
    offset(:,i) = SGNdFoF(:,i) - (i-1)*1;
    end
    figure; plot(offset);
    
    dontInc = 200;
    stdDFOF = std(SGNdFoF(1:end-dontInc,:),1,1);
    meanDFOF = mean(SGNdFoF(1:end-dontInc,:),1);

    %plot(meanDFOF,'Color','g');
    thr = meanDFOF + 2.5*stdDFOF;
    thrDFOF = zeros(t,n);
    for i=1:n
        thrDFOF(:,i) = SGNdFoF(:,i) > thr(i);
        %check for cells that come up for seconds and remove those signals
%         labeled = bwlabel(thrDFOF(1:end-dontInc,i));
%         numEvents = max(labeled);
%         for j=1:numEvents
%             if size(find(labeled == j),1) > 10
%                 labeled(labeled == j) = 0;
%                 thrDFOF(:,i) = zeros(size(thrDFOF(:,i)));
%                 disp('Gotcha');
%             end
%         end
        %thrDFOF(:,i) = thrDFOF(:,i) .* (labeled > 0);    
    end
    sumThr = sum(thrDFOF,2);
    
    [baselinePks] = findpeaks(sumThr(baseline(1):baseline(2)),'MinPeakHeight',numCellsForEvent,'MinPeakDistance',2);
    ROIstats.baselineFreq = size(baselinePks,1)/((baseline(2)-baseline(1))/60);
    ROIstats.baselinePks = baselinePks;
    if drugs1 ~= 0
        [drug1Pks] = findpeaks(sumThr(drugs1(1):drugs1(2)),'MinPeakHeight',numCellsForEvent,'MinPeakDistance',2);
        ROIstats.drugs1Freq = size(drug1Pks,1)/((drugs1(2)-drugs1(1))/60); 
        ROIstats.drug1Pks = drug1Pks;
        if drugs2 ~= 0
            [drug2Pks] = findpeaks(sumThr(drugs2(1):drugs2(2)),'MinPeakHeight',numCellsForEvent,'MinPeakDistance',2);
            ROIstats.drugs2Freq = size(drug2Pks,1)/((drugs2(2)-drugs2(1))/60);
            ROIstats.drug2Pks = drug2Pks;
        end    
    end
    
    
    figure; subplot(2,1,1); plot(sum(thrDFOF,2))
    subplot(2,1,2); plotSpikeRaster(logical(thrDFOF)','PlotType','vertline');


    [chi, chiShuf, inputShuf] = synchronyIndex(SGNdFoF(baseline(1):baseline(2),:));
    figure; imagesc([corr(SGNdFoF(1:500,:)) corr(inputShuf)]);
    caxis([0 1]);
  
    ROIstats.chi = chi;
    ROIstats.chiShuf = chiShuf;

end

function [chi, chiShuf, inputShuf] = synchronyIndex(input)
%synchrony index computed as outlined @ http://www.scholarpedia.org/article/Neuronal_synchrony_measures
    input = input;
    [t,n] = size(input);
    Vt = mean(input,2);

    %shuffle data;
    R = randi(t,n,1);
    inputShuf = zeros(t,n);
    for i = 1:n
        inputShuf(:,i) = [input(R(i):end,i); input(1:R(i)-1,i)];   
    end
    VtShuf = mean(inputShuf,2);

    time = [1:t]';
    sig = trapz(time,Vt.^2)/time(end) - (trapz(time,Vt)/time(end)).^2;
    sigSing = trapz(time,input.^2)/time(end) - (trapz(time,input)/time(end)).^2;
    sigSing = mean(sigSing);
    chi2 = sig/sigSing;
    chi = sqrt(chi2)

    sigShuf = trapz(time,VtShuf.^2)/time(end) - (trapz(time,VtShuf)/time(end)).^2;
    sigSingShuf = trapz(time,inputShuf.^2)/time(end) - (trapz(time,inputShuf)/time(end)).^2;
    sigSingShuf = mean(sigSingShuf);
    chi2shuf = sigShuf/sigSingShuf;
    chiShuf = sqrt(chi2shuf)

    figure; plot(Vt);

end

