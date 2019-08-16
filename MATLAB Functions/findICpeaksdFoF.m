function [stats, pkData] = findICpeaksdFoF(ICsignal,filePath,analysisName,plotFlag)
    if ~exist('analysisName','var')
        analysisName = 'stock';
    end
    if ~exist('plotFlag','var')
        plotFlag = 1;
    end
    
    time = [1:1:size(ICsignal,1)]';
    LIC = smooth(ICsignal(:,1));
    RIC = smooth(ICsignal(:,2));
    ctx = smooth(ICsignal(:,3));
    LIC = msbackadj(time,LIC,'WindowSize',500,'StepSize',500);
    RIC = msbackadj(time,RIC,'WindowSize',500,'StepSize',500);
    ctx = msbackadj(time,ctx,'WindowSize',500,'StepSize',500);
    
    %parameters
    pkThreshold = .01;
    pkMinHeight = .01;
    pkDistance = 5; %in frame, 10 = 1s
    [pks,locs,w] = findpeaks(LIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkMinHeight,'MinPeakDistance',pkDistance,'Annotate','extents');
    LIC_itk = ctx_PkRemoval(ctx, LIC, locs);
    pksLIC = pks(LIC_itk);
    locsLIC = locs(LIC_itk);
    wLIC = w(LIC_itk);
    LICinfo = [pksLIC locsLIC wLIC];
    
    
    [pks,locs,w] = findpeaks(RIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkMinHeight,'MinPeakDistance',pkDistance,'Annotate','extents');
    RIC_itk = ctx_PkRemoval(ctx, RIC, locs);
    pksRIC = pks(RIC_itk);
    locsRIC = locs(RIC_itk);
    wRIC = w(RIC_itk);
    RICinfo = [pksRIC locsRIC wRIC];
    pkData = ICcompare2(LIC, RIC, pksLIC,locsLIC,pksRIC,locsRIC);
    
    if(plotFlag)
        figure;
        subplot(3,1,1);
            plot(time/10,LIC);
            hold on;
            plot(time/10,RIC);
            plot(time/10,ctx);
            plot(locsRIC/10,pksRIC,'o');
            plot(locsLIC/10,pksLIC,'o');
            ylabel('Amplitude (AIU)')
            xlabel('Time (s)');
        subplot(3,1,2);
            plot(time/10,ctx);
            ylabel('Amplitude (AIU)')
            xlabel('Time (s)');
        subplot(3,1,3);
            plot(time/10,ctx);
            hold on;
            plot(time/10,ctx);
        %disp('plot active');
        graphEvents(pkData,[filePath,'timevert_',analysisName]);
    end
    %savefig([filePath,'dfof_bg_',analysisName]);
    
    assignin('base','pkdata',pkData);
    
    %get stats about events
    totLpks = size([pkData(pkData(:,6)==1); pkData(pkData(:,6)==2)],1);
    totRpks = size([pkData(pkData(:,6)==1); pkData(pkData(:,6)==3)],1);
    totMatchedPks = size(pkData(pkData(:,6)==1),1);
    totPks = totLpks + totRpks - totMatchedPks;
    Ldom = size(pkData(pkData(:,7)==1),1); % dominant LIC peaks
    Rdom = size(pkData(pkData(:,7)==2),1) %dominant RIC peaks
    meanR = mean(pkData(:,4));
    meanL = mean(pkData(:,2));
    width = mean([RICinfo(:,3); LICinfo(:,3)],1)
    %calculate event ratios
    tempPkData = pkData(pkData(:,7)==1,:)
    ratio = tempPkData(:,4)./tempPkData(:,2);
    tempPkData = pkData(pkData(:,7)==2,:);
    ratio = [ratio; tempPkData(:,2)./tempPkData(:,4)];
    
    disp(['LIC: ',num2str(totLpks),' peaks total, ',num2str(totMatchedPks),' corresponding RIC peaks.']);
    disp(['RIC: ',num2str(totRpks),' peaks total, ',num2str(totMatchedPks),' corresponding LIC peaks.']);
    disp(['Mean RIC amplitude: ', num2str(meanR)]);
    disp(['Mean LIC amplitude: ', num2str(meanL)]);
    stats = table(totLpks, totRpks, totMatchedPks, meanL, meanR, Ldom/totPks, totLpks+totRpks-totMatchedPks, Ldom, Rdom, width, mean(ratio));
    assignin('base','stats1',stats);
    save([filePath,['\ICinfo16_',analysisName]],'LICinfo','RICinfo','ICsignal','filePath','stats','pkData');
    disp([filePath,['\ICinfo16_',analysisName]]);
end

function [pkData] = ICcompare2(LIC, RIC, master_pks,master_locs,pks,locs)
    windowforhit = 10;
    
    pkData = []; %col 1: LIC loc, 2: LIC pk, 3:RIC loc, 4:RIC pk, 5:delta, 6: peak type (1=matched, 2=LIC only, 3=RIC only) 7: which is bigger (1=LIC, 2=RIC)
    locs = [locs(:) zeros(size(locs,1),1)];
    for i=1:size(master_pks,1)
        match = 0;
        j = 1;
        while match == 0
            if abs(master_locs(i)-locs(j,1)) <= windowforhit && locs(j,2) ~= 1
                pkData(i, 1:6) = [master_locs(i) master_pks(i) locs(j,1) pks(j) abs(pks(j)-master_pks(i)) 1];
                locs(j,2) = 1;
                match = 1;
            else
                j = j+1;
            end
            
            if j > size(pks,1)
                match = 1;
                pkData(i, 1:6) = [master_locs(i) master_pks(i) master_locs(i) RIC(master_locs(i)) abs(master_pks(i)-RIC(master_locs(i))) 2];
                %pkData(i, 1:5) = [master_locs(i) master_pks(i) master_locs(i) 0 pks(i)];
            end
        end
    end
   
    %insert unmatched pks
    unmatchedPks = pks(locs(:,2)==0);
    unmatchedLocs = locs(locs(:,2) == 0);
    for i=1:size(unmatchedPks,1)
        loc_under = 1;
        j=1;
        while loc_under == 1
          if j > size(pkData,1)
              temp_pkdata = [unmatchedLocs(i) LIC(unmatchedLocs(i)) unmatchedLocs(i) unmatchedPks(i) abs(unmatchedPks(i)-LIC(unmatchedLocs(i))) 3];
              pkData = [pkData; temp_pkdata];
              loc_under = 0;
          end          
          if pkData(j, 3) < unmatchedLocs(i)
              j = j + 1;
          else
              temp_pkdata = [unmatchedLocs(i) LIC(unmatchedLocs(i)) unmatchedLocs(i) unmatchedPks(i) abs(unmatchedPks(i)-LIC(unmatchedLocs(i))) 3];
              %temp_pkdata = [unmatchedLocs(i) 0 unmatchedLocs(i) unmatchedPks(i) unmatchedPks(i)];
              pkData = [pkData(1:j-1,:); temp_pkdata; pkData(j:end,:)];
              loc_under = 0;
          end   
        end
    end
    
    %compute which side is bigger
    for i=1:size(pkData,1)
        if pkData(i,2) > pkData(i,4)
            pkData(i,7) = 1;
        else
            pkData(i,7) = 2;
        end
    end 
end

function [indexToKeep] = ctx_PkRemoval(ctxSignal, ICsignal, IClocs)
    time = [1:1:size(ctxSignal,1)]';
    
    %[pks,locs,w] = findpeaks(msbackadj(time,ctxSignal),'MinPeakHeight',0.03,'Annotate','extents');
    
    ctxBright = ctxSignal > ICsignal;
    
    indexToKeep = [];
    for i=1:size(IClocs)
        if ctxBright(IClocs(i)) == 0
            indexToKeep = [indexToKeep i];
        end      
    end
end

function [events] = graphLocsOnly(pks1,locs1,pks2,locs2)
    windowforhit = 10;
    
    result_index = 1;
    events = []; %small 1, small 2, big 1, big 2
    while ~isempty(locs1) & ~isempty(locs2)
        if isempty(locs1)
           events(result_index,4) = pks2(1);
           events(result_index,6) = 1;
           events(result_index,1) = -1;
           locs2(1) = [];
           pks2(1) = [];
           result_index = result_index + 1;
        elseif isempty(locs2)
           events(result_index,3) = pks1(1);
           events(result_index,5) = 1;
           events(result_index,2) = -1;
           locs1(1) = [];
           pks1(1) = [];
           result_index = result_index + 1;
        else
            if locs1(1) < locs2(1)
                %check if match
                match = 0;
                j = 1;
                while match == 0
                    if abs(locs1(1)-locs2(j)) <= windowforhit
                        %events(result_index,2) = pks2(j);
                        if pks2(j) > pks1(1)
                            events(result_index,4)=pks2(j);
                            events(result_index,1)=pks1(1);
                            events(result_index,6) = 1;
                        else
                            events(result_index,2)=pks2(j);
                            events(result_index,3)=pks1(1);
                            events(result_index,5) = 1;
                        end
                        locs2(j) = [];
                        pks2(j) = [];
                        match = 1;
                    else
                        j = j+1;
                    end
                    if j > size(locs2,1)
                        match = 1;
                        events(result_index,3)=pks1(1);
                        events(result_index,5) = 1;
                        events(result_index,2)=-1;
                    end
                end
                locs1(1) = [];
                pks1(1) = [];
            else
                %check if match
                match = 0;
                j = 1;
                while match == 0
                    if abs(locs1(j)-locs2(1)) <= windowforhit
                        if pks2(j) > pks1(1)
                            events(result_index,4)=pks2(1);
                            events(result_index,1)=pks1(j);
                            events(result_index,6) = 1;
                        else
                            events(result_index,2)=pks2(1);
                            events(result_index,3)=pks1(j);
                            events(result_index,5) = 1;
                        end
                        locs1(j) = [];
                        pks1(j) = [];
                        match = 1;
                    else
                        j = j+1;
                    end
                    if j > size(locs1,1)
                        match = 1;
                        events(result_index,4)=pks2(1);
                        events(result_index,6) = 1;
                        events(result_index,1)=-1;
                    end
                end
                locs2(1) = [];
                pks2(1) = [];
            end  
            result_index = result_index + 1;
        end
    end
   
    events(events==0)=NaN;
    events(events==-1)=0;
end

function graphEvents(pkData, filePath)
    LIC_big = pkData(pkData(:,7)==1,:);
    RIC_big = pkData(pkData(:,7)==2,:);
    lt_org = [255, 166 , 38]/255;
    dk_org = [255, 120, 0]/255;
    lt_blue = [50, 175, 242]/255;
    dk_blue = [0, 13, 242]/255;
    %figure;
        %plot(-LIC_big(:,2),-LIC_big(:,1),'o','Color','red','MarkerFaceColor','red');
        %hold on;
        %plot(LIC_big(:,4),-LIC_big(:,3),'o','Color','black');
        %plot(-RIC_big(:,2),-RIC_big(:,1),'o','Color','red');
        %plot(RIC_big(:,4),-RIC_big(:,3),'o','Color','black','MarkerFaceColor','black');
        %all points with lines
       % line([-pkData(:,2)'; zeros(size(pkData,1),1)'],[-pkData(:,1)'; -pkData(:,1)'],'Color','red');
        %line([zeros(size(pkData,1),1)'; pkData(:,4)'; ],[-pkData(:,3)'; -pkData(:,3)'],'Color','black');
        %only bigger events with lines
        %line([-LIC_big(:,2)'; zeros(size(LIC_big,1),1)'],[-LIC_big(:,1)'; -LIC_big(:,1)'],'Color','red');
        %line([zeros(size(RIC_big,1),1)'; RIC_big(:,4)'; ],[-RIC_big(:,3)'; -RIC_big(:,3)'],'Color','black');
   % savefig(filePath);
    
    %h = figure(3);
    h = figure;
    ax = gca;
        line([(-pkData(:,2))'; zeros(size(pkData,1),1)'],[-pkData(:,1)'; -pkData(:,1)'],'Color',lt_org);
        hold on;
        line([zeros(size(pkData,1),1)'; (pkData(:,4))'; ],[-pkData(:,3)'; -pkData(:,3)'],'Color',lt_blue);
        line([0 0]',[-6050 50]','Color','black');
        scatter(pkData((pkData(:,7)==2),4),-pkData((pkData(:,7)==2),3),pkData(pkData(:,7)==2,5)*200,dk_blue,'MarkerFaceColor', dk_blue);
        scatter(-pkData((pkData(:,7)==1),2),-pkData((pkData(:,7)==1),1),pkData(pkData(:,7)==1,5)*200,dk_org,'MarkerFaceColor', dk_org);
        set(h,'Position',[200,0,350,500]);
        ax.XLim = [-0.4 0.4];
        ax.YLim = [-6050 50];
        %plot(pkData((pkData(:,7)==2),4),-pkData((pkData(:,7)==2),3),'o','Color','black');
        %plot(-pkData((pkData(:,7)==1),2),-pkData((pkData(:,7)==1),1),'o','Color','red');
        
        %plot(find(pkData(:,7)==2),0,pkData(pkData(:,7)==2,5),'o','Color','black','MarkerFaceColor','black');
        %plot(find(pkData(:,7)==1),0,'o','Color','red','MarkerFaceColor','red');
        
        Rpks = pkData((pkData(:,7)==2),4);
        Lpks = pkData((pkData(:,7)==1),2);
        hold off;
        binLim = [0 0.02:.05:.35 100];
        Rcounts= histcounts(Rpks,binLim);
        Lcounts =histcounts(Lpks,binLim);
        binY = [0 .0450:.05:.35];
        
        figh = figure;
            h=barh(binY(2:end),Rcounts(2:end),.9);
            hold on;
            barh(binY(2:end),-Lcounts(2:end),.9,'FaceColor',lt_org,'EdgeColor','none');
             line([0 0],[0 0.4],'LineWidth',0.75,'Color',[0.6 0.6 0.6]);
            %barh(Lbins,-Lcounts,.9,'FaceColor',lt_org,'EdgeColor',lt_org);
            h.FaceColor = lt_blue;
            h.EdgeColor = 'none';
           
           xax_lim = [-60 60];
            xtick = [-60 -30 0 30 60];
            xticklabel = {'60' '30' '0' '30' '60'};
             xlim(xax_lim);
            xticks(xtick);
            xticklabels(xticklabel);
            ylim([0 .40]);
            yticks([0 .2 .4]);
            yticklabels({'0' '20' '40'});
            box off;
            xlabel('# of Dominant Events','FontSize',8);
            ylabel('\DeltaF/F (%)','FontSize',8);
            axh = gca;
            figQuality(figh,axh,[1.8 1.2]);
        
end

