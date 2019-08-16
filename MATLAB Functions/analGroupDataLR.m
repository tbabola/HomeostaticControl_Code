function [groupHistoL,groupHistoR,totalEvents,individualStats,histoLsingle, histoRsingle] = analGroupDataLR(paths,plotFlag,normType)
    if nargin < 3 || isempty(normType), normType = 'count'; end
    if nargin < 2 || isempty(plotFlag), plotFlag = 1; end;
    %paths = paths(:,3);
    groupHistoL = [];
    groupHistoR = [];
    totalEvents = [];
    individualStats = [];
    hwtL = [];
    hwtR = [];
    hwxL = [];
    hwxR = [];
    ampL = [];
    ampR = [];
    
    for i=1:size(paths)
        load(paths{i},'-regexp','^(?!LICmov|RICmov)...');
%         load(paths{i});
        [peaksBinaryL] = getPeaks_dFoF(smLIC,0,0.02);
         [peaksBinaryR] = getPeaks_dFoF(smRIC,1,0.02);
%         savefile = '\ICmovs_peaks.mat';
%         [path, fn, ext ] = fileparts(paths{i});
%         save([path savefile],'LICmov','RICmov','smLIC','smRIC','peaksBinaryR','peaksBinaryL');
        %size(peaksBinaryR)
        %size(peaksBinaryL)
        [peakStat, eventStats] = peakStats_dFoF(smLIC, peaksBinaryL, smRIC, peaksBinaryR);
       
        %eventStat = table(eventLabel, eventClassification, leftOrRightDom, numPeaks, domAmp, maxLAmp, maxRAmp,xloc, tloc, hwt, hwx);
        peaksL = peakStat{1};
        peaksR = peakStat{2};
        
        groupHistoL = [groupHistoL; peaksL(:,2)];
        groupHistoR = [groupHistoR; peaksR(:,2)];
        histoLat = [peaksL(:,2); abs(peaksR(:,2)-125)];
        histoLat = histoLat(histoLat > 50);
        histoLat = histoLat(histoLat < 100);
        totalEvents = [totalEvents; eventStats];
        %individualStats(i,:) = [max(eventStats.eventLabel) median(eventStats.domAmp) median(eventStats.hwt) nanmedian(eventStats.hwx) size((find(eventStats.numPeaks > 1)),1)/max(eventStats.eventLabel)*100 size(histoLat,1)];
        
        %individual histograms
        edges = 1:7:125;
        leftPeaks = peaksL(:,2);
        leftPeaks = leftPeaks(leftPeaks > 20 & leftPeaks < 100);
        rightPeaks = abs(peaksR(:,2)-130);
        rightPeaks = rightPeaks(rightPeaks > 20 & rightPeaks < 100);
        hl = histogram(leftPeaks,edges,'Normalization',normType);
        histoLsingle(:,i) = hl.Values';
        hr = histogram(rightPeaks,edges,'Normalization',normType);
        histoRsingle(:,i) = hr.Values';
        
        
        %left right compariason
        hwtLdist = eventStats.hwt(strcmp('Left', cellstr(eventStats.leftOrRightDom)));
        hwtRdist = eventStats.hwt(strcmp('Right', cellstr(eventStats.leftOrRightDom)));
        hwtL = [hwtL; hwtLdist];
        hwtR = [hwtR; hwtRdist];
        hwtLavg = median(hwtLdist);
        hwtRavg = median(hwtRdist);
        
        hwxLdist = eventStats.hwx(strcmp('Left', cellstr(eventStats.leftOrRightDom)));
        hwxRdist = eventStats.hwx(strcmp('Right', cellstr(eventStats.leftOrRightDom)));
        hwxL = [hwxL; hwxLdist];
        hwxR = [hwxR; hwxRdist];
        hwxLavg = nanmedian(hwxLdist);
        hwxRavg = nanmedian(hwxRdist);
        
        ampLdist = eventStats.domAmp(strcmp('Left', cellstr(eventStats.leftOrRightDom)));
        ampRdist = eventStats.domAmp(strcmp('Right', cellstr(eventStats.leftOrRightDom)));
        ampL = [ampL; ampLdist];
        ampR = [ampR; ampRdist];
        ampLavg = median(ampLdist);
        ampRavg = median(ampRdist);
        if ~isnan(max(eventStats.eventLabel))
            individualStats(i,:) = [max(eventStats.eventLabel) median(eventStats.domAmp) median(eventStats.hwt) nanmedian(eventStats.hwx) size((find(eventStats.numPeaks > 1)),1)/max(eventStats.eventLabel)*100 size(histoLat,1)...
                hwtLavg hwtRavg hwxLavg hwxRavg ampLavg ampRavg 0 0 0];
            %kstest2(hwtLdist,hwtRdist) kstest2(hwxLdist,hwxRdist) kstest2(ampLdist,ampRdist)];
        else
            individualStats(i,:) = zeros(1,15);
        end
        
    
    end
    
    if plotFlag
          %plotTimeSeries(smLIC, smRIC, peaksBinaryL, peaksBinaryR, peakStat);
          compareFigLR([hwtL hwxL ampL], [hwtR hwxR ampR],groupHistoL, groupHistoR);
    end
    
end

function [groupHistoL,groupHistoR,totalEvents,individualStats, individualPks] = analGroupData(paths,plotFlag)
    groupHistoL = [];
    groupHistoR = [];
    totalEvents = [];
    individualStats = [];

    for i=1:size(paths)
         %load(paths{i},'-regexp','^(?!LICmov|RICmov)...');
         load(paths{i});
         [peaksBinaryL] = getPeaks_dFoF(smLIC,0,0.02);
         [peaksBinaryR] = getPeaks_dFoF(smRIC,1,0.02);
         savefile = '\ICmovs_peaks_dFoF.mat';
         [path, fn, ext ] = fileparts(paths{i});
         save([path savefile],'smLIC','smRIC','peaksBinaryR','peaksBinaryL');
        %size(peaksBinaryR)
        %size(peaksBinaryL)
        
        [peakStat, eventStats] = peakStats_dFoF(smLIC, peaksBinaryL, smRIC, peaksBinaryR);
       
        %eventStat = table(eventLabel, eventClassification, leftOrRightDom, numPeaks, domAmp, maxLAmp, maxRAmp,xloc, tloc, hwt, hwx);
        peaksL = peakStat{1};
        peaksR = peakStat{2};
        groupHistoL = [groupHistoL; peaksL(:,2)];
        groupHistoR = [groupHistoR; peaksR(:,2)];
        histoLat = [peaksL(:,2); abs(peaksR(:,2)-125)];
        histoLat = histoLat(histoLat > 30);
        histoLat
        histoLat = histoLat(histoLat < 50);
        histoMed = [peaksL(:,2); abs(peaksR(:,2)-125)];
        histoMed = histoMed(histoMed > 50);
        histoMed = histoMed(histoMed < 70);
        totalEvents = [totalEvents; eventStats];
        [max(eventStats.eventLabel) median(eventStats.domAmp) median(eventStats.hwt) nanmedian(eventStats.hwx) size((find(eventStats.numPeaks > 1)),1)/max(eventStats.eventLabel)*100 size(histoLat,1) median(eventStats.integral)]
        individualStats(i,:) = [max(eventStats.eventLabel) median(eventStats.domAmp) median(eventStats.hwt) nanmedian(eventStats.hwx) size((find(eventStats.numPeaks > 1)),1)/max(eventStats.eventLabel)*100 size(histoLat,1) median(eventStats.integral) size(histoMed,1)];
        individualPks{i} = [peaksL(:,2); abs(peaksR(:,2)-125)];
        
        if(plotFlag)
          plotTimeSeries_dFoF(smLIC, smRIC, peaksBinaryL, peaksBinaryR, peakStat);
        end
    
    end
    
end

function compareFigLR(lEvents, rEvents,lHisto, rHisto)
    figure;
    cond1 = 'Left';
    cond2 = 'Right';
    
    subplot(3,1,1);
    cdfplot(lEvents(:,2)*7.45);
    hold on;
    cdfplot(rEvents(:,2)*7.45);
    xlabel('Spatial Spread (\mum)');
    legend(cond1,cond2);
    [h,p]=kstest2(lEvents(:,2),rEvents(:,2))
    text(80,0.5, strcat('p=',num2str(p)))
    
    subplot(3,1,2);
    cdfplot(lEvents(:,1)/10);
    hold on;
    cdfplot(rEvents(:,1)/10);
    xlabel('Half-width (s)');
    legend(cond1,cond2);
    [h,p]=kstest2(lEvents(:,1),rEvents(:,1))
    text(30,0.5,  strcat('p=',num2str(p)))
    
    subplot(3,1,3);
    cdfplot(lEvents(:,3));
    hold on;
    cdfplot(rEvents(:,3));
    xlabel('Amplitude (AIU)');
    legend(cond1,cond2);
    [h,p]=kstest2(lEvents(:,3),rEvents(:,3))
    text(35,0.5,  strcat('p=',num2str(p)))
    
    %%more plot
    figure;
    histogram([lHisto],25,'faceColor','black','edgeColor','black','faceAlpha',1)
    xlim([0,125]);
    hold on;
    histogram([rHisto],25,'faceColor','red','edgeColor','black','faceAlpha',1)
    
end

function compareStats(cond1, indStats, cond2, indStats2, plotFlag)
    m = size(indStats,1);
    m2 = size(indStats2,1);
    
    if(plotFlag)
        figure;
        %frequency
        subplot(1,3,1);
            scatter(ones(m,1), indStats(:,1)/10)
            hold on;
            scatter(2*ones(m2,1), indStats2(:,1)/10)
            ylim([0,25]);
            xlim([0 3]);
            xticklabels({'',cond1,cond2,''});
            ylabel('Frequency (events/min)');
        %amplitude
        subplot(1,3,2);
            scatter(ones(m,1), indStats(:,2))
            hold on;
            scatter(2*ones(m2,1), indStats2(:,2))
            ylim([0,25]);
            xlim([0 3]);
            xticklabels({'',cond1,cond2,''});
            ylabel('Amplitude (AIU)');
        subplot(1,3,3);
            scatter(ones(m,1), indStats(:,4)*7.45)
            hold on;
            scatter(2*ones(m2,1), indStats2(:,4)*7.45)
            ylim([0,800]);
            xlim([0 3]);
            xticklabels({'',cond1,cond2,''});
            ylabel('Spatial Spread (\mum)');
        
    end
end