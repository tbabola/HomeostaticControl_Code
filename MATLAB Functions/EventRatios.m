lt_org = [255, 166 , 38]/255;
dk_org = [255, 120, 0]/255;
lt_blue = [50, 175, 242]/255;
dk_blue = [0, 13, 242]/255;


TMIEFiles = loadFileList('..\..\TMIE\Fig. 8 - TMIE KO\Data\TMIE KO\*\*ICinfo16_dFoF.mat');
TMIEstats = flipFlopStats(TMIEFiles,0);

wtFiles = loadFileList('.\Data\WILDTYPE\*\*ICinfo16_dFoF.mat');
WTstats = flipFlopStats(wtFiles,0);

lablFiles = loadFileList('.\Data\Left Ablation\*\*ICinfo16_dFoF.mat');
Lablstats = flipFlopStats(lablFiles,0);

rablFiles = loadFileList('.\Data\Right Ablation\*\*ICinfo16_dFoF.mat');
Rablstats = flipFlopStats(rablFiles,1);

VG3Files = loadFileList('..\Figure 5 - VG3 KO\Data\VG3 KO\*\*ICinfo16_dFoF.mat');
VG3stats = flipFlopStats(VG3Files,0);

c1 =plotThings(WTstats.means',[Lablstats.means Rablstats.means]', VG3stats.means', TMIEstats.means','Mean Ratio', [0 1]);
c2 =plotThings(WTstats.cvs',[Lablstats.cvs Rablstats.cvs]', VG3stats.cvs', TMIEstats.cvs','Coefficient of Variation',[0 100]);
%plotThings(WTstats.normality',Lablstats.normality', Rablstats.normality', VG3stats.normality', TMIEstats.normality','Mean Ratio', [0 1]);
%plotThings(WTstats.means',Lablstats.means', Rablstats.means', VG3stats.means', TMIEstats.means','Mean Ratio');

function c = plotThings(thing1, thing2, thing3, thing4, ylbl, ylimits)
    figure;
    plot(ones(size(thing1,1),1),thing1,'.','Color',[0.7 0.7 0.7]); hold on;
    errorbar(1,mean(thing1,1),sterr(thing1,1),'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10);
    plot(2*ones(size(thing2,1),1),thing2,'.','Color',[0.7 0.7 0.7]);
    errorbar(2,mean(thing2,1),sterr(thing2,1),'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10);
    plot(3*ones(size(thing3,1),1),thing3,'.','Color',[0.7 0.7 0.7]);
    errorbar(3,mean(thing3,1),sterr(thing3,1),'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10);
    plot(4*ones(size(thing4,1),1),thing4,'.','Color',[0.7 0.7 0.7]);
    errorbar(4, mean(thing4,1),sterr(thing4,1),'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10);

    xlim([0 5]);
    ylim(ylimits);
    ylabel(ylbl);
    xticks([1:1:4]);
    xticklabels({'WT', 'Uni ablations', 'VG3 KO', 'TMIE KO'});
    box off;
    figQuality(gcf,gca,[3 2]);
    export_fig(['.\EPS Panels\FlipFlop\' ylbl '.eps'], '-eps', '-nocrop')
    [p,anovatab, stats] = anova1([thing1; thing2; thing3; thing4],[ones(size(thing1,1),1); 2*ones(size(thing2,1),1); 3*ones(size(thing3,1),1); 4*ones(size(thing4,1),1)])
    c = multcompare(stats);
    
end



function [statsG, delta] = flipFlopStats(files, plotFlag)
    m = size(files,1)
    paths = files;
    if plotFlag
        h = figure('Position',[500 -500 1200 1200]); 
        h2 = figure;
    end
    means = [];
    stds = [];
    cvs = [];
    histoBinCounts = [];
    statsG = struct();
    for i = 1:m
             thr = .01
             load(paths{i});
             disp(paths{i});
             [filePath name ext] = fileparts(paths{i});
             [stats, pkData] = findICpeaksdFoF(ICsignal,filePath,'dFoF',0);
             temp = pkData(pkData(:,7)==1 & pkData(:,2) > thr,:);
             deltaCalcLDOM = [temp(:,2) temp(:,4)];
             time = temp(:,1);
             delL = deltaCalcLDOM(:,2)./deltaCalcLDOM(:,1);

             temp = pkData(pkData(:,7)==2 & pkData(:,4) > thr,:);
             deltaCalcRDOM = [temp(:,4) temp(:,2)];
             delR = deltaCalcRDOM(:,2)./deltaCalcRDOM(:,1);
             time2 = temp(:,3);

             delta = [deltaCalcRDOM; deltaCalcLDOM];
             delta = delta(:,2)./delta(:,1);
             timeTot = [time; time2];
             [timeTot, o] = sort(timeTot,'ascend');
             delta = delta(o);
             if plotFlag
                 figure(h);
                 subplot(m,1,i);
                 plot(timeTot,abs(smooth(delta,1)),'.-','Color',[0.7 0.7 0.7]); hold on;
                 ylim([0 1]);
                 xlim([0 6000]);
                 box off;

                 figure(h2);
                 subplot(m,1,i);
                 histo = histogram(delta,[0:.10:1]);
             end
             means(i) = mean(delta);
             stds(i) = std(delta,[],1);
             cvs(i) = stds(i)/means(i)*100;
             normality(i) = kstest(delta);
             histoBinCounts(i,:) = histcounts(delta,[0:.10:1]);


    %          figure; scatter(deltaCalcLDOM(:,1),delL, 15, lt_org); hold on;
    %          scatter(deltaCalcRDOM(:,1),delR,15, lt_blue);
    %          ylim([0 1]);
    %          xlim([0 inf]);
    %          ylabel('Ratio');
    %          xlabel('Amplitude of Dominant Side');
    %          corrcoef(deltaCalcLDOM(:,1), delL)
    %          corrcoef(deltaCalcRDOM(:,1), delR)
    %          
    %          figure; scatter(deltaCalcLDOM(:,1),deltaCalcLDOM(:,2),10, lt_org); hold on;
    %          scatter(deltaCalcRDOM(:,1),deltaCalcRDOM(:,2),10, lt_blue);
    %          xlim([0 .35]);
    %          ylim([0 .35]);

    %          figure; plot(time,smooth(abs(deltaCalcLDOM(:,1)-deltaCalcLDOM(:,2)),5),'-','MarkerSize',10,'Color',lt_org); hold on;
    %          %scatter(time,abs(deltaCalcLDOM(:,2)./deltaCalcLDOM(:,1)),abs(deltaCalcLDOM(:,1))*100, lt_org);
    %          plot(time2,smooth(abs(deltaCalcRDOM(:,1)-deltaCalcRDOM(:,2)),5),'-','MarkerSize',10,'Color',lt_blue);
    %          %scatter(time2,abs(deltaCalcRDOM(:,2)./deltaCalcRDOM(:,1)),abs(deltaCalcRDOM(:,1))*100, lt_blue);
    %          ylim([0 1]);
    %          ylabel('Ratio');
    %          xlabel('Time (s)');
    end
    statsG.means = means;
    statsG.stds = stds;
    statsG.cvs = cvs;
    statsG.normality = normality;
    statsG.histoBinCounts = histoBinCounts;
end
     
%close all;
