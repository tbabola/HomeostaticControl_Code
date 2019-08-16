function [avgWave, stErr, half_width, ampl] = waveform2LR(paths)
    m = size(paths,1);
    
    win_before = 20;
    win_after = 50;
    LICavgTraces = [];
    RICavgTraces = [];
    half_widthL = [];
    half_widthR = [];
    amplL = [];
    amplR = [];
    
    for i=1:m
        load(paths{i});
        ICsignal = [smooth(ICsignal(:,1)) smooth(ICsignal(:,2)) smooth(ICsignal(:,3))];
        t = size(ICsignal,1);
        l_locs = pkData(pkData(:,7)==1,1);
        l_m = size(l_locs);
        r_locs = pkData(pkData(:,7)==2,3);
        r_m = size(r_locs);
        k = 1;
        LIC_trace = [];
        half_widthR = [];
        amplR = [];
        half_widthR = [];
        amplR = [];
        
        for j=1:l_m
            loc = l_locs(j);
            if ~(loc <= win_before || (loc + win_after) > t)
                loc_l = loc - win_before;
                loc_r = loc + win_after;
                LIC_trace(:,k) = ICsignal(loc_l:loc_r,1);
                [pks,locs,w] = findpeaks(ICsignal(loc_l:loc_r,1),[loc_l:1:loc_r]);
                if size(find(locs == loc)) > 0
                    half_widthL(end+1) = w(find(locs == loc));
                    amplL(end+1) = pks(find(locs == loc));
                end
                 k = k+1;
            end    
        end
        
        k = 1;
        RIC_trace = [];
        for j=1:r_m
            loc = r_locs(j);
            if ~(loc <= win_before || (loc + win_after) > t)
                loc_l = loc - win_before;
                loc_r = loc + win_after;
                RIC_trace(:,k) = ICsignal(loc_l:loc_r,2);
                [pks,locs,w] = findpeaks(ICsignal(loc_l:loc_r,2),[loc_l:1:loc_r]);
                if size(find(locs == loc)) > 0
                    half_widthR(end+1) = w(find(locs == loc));
                    amplR(end+1) = pks(find(locs == loc));
                end
                k = k+1;
            end    
        end
       
        LICavgTraces(:,i) = mean(LIC_trace,2);
        RICavgTraces(:,i) = mean(RIC_trace,2);
        totavgTraces(:,i) = mean([LIC_trace RIC_trace],2);
        totavgTraces(:,i) = totavgTraces(:,i) - totavgTraces(1,i);
        [pks, locs, w] = findpeaks(totavgTraces(:,i),[1:1:win_before+win_after+1],'MinPeakHeight',0.02);
        half_width(i,:) = [mean(half_widthL), mean(half_widthR)];
        ampl(i,:) = [mean(amplL), mean(amplR)];
        
        
    end
    
    time = [1:1:win_before+win_after+1]'/10;
    lt_org = [255, 166 , 38]/255;
    dk_org = [255, 120, 0]/255;
    lt_blue = [50, 175, 242]/255;
    dk_blue = [0, 13, 242]/255;
    
    
    if m == 1
            totTraces = [LIC_trace, RIC_trace];
            f = figure; 
%             set(f,'renderer','painters');
%             shadedErrorBar(time, mean(totTraces,2), std(totTraces,1,2)/sqrt(size(totTraces,2)),{'-','color','r'});
%             xlim([0 6.5]);
%             ylim([-0.01 .15]);
            avgWave = mean(totTraces,2)-mean(totTraces(1,:),2);
            stErr = std(totTraces,1,2)/sqrt(i);
    else
%         f = figure; 
%         set(f,'renderer','painters');
%         shadedErrorBar(time, mean(LICavgTraces,2)-mean(LICavgTraces(1,:),2), std(LICavgTraces,1,2),{'-','color',dk_org}); hold on;
%         shadedErrorBar(time, mean(RICavgTraces,2)-mean(RICavgTraces(1,:),2), std(RICavgTraces,1,2),{'-','color',dk_blue});
        avgWave = [mean(LICavgTraces,2)-mean(LICavgTraces(1,:),2) mean(RICavgTraces,2)-mean(RICavgTraces(1,:),2)];
        stErr = [std(LICavgTraces,1,2)/sqrt(m) std(RICavgTraces,1,2)/sqrt(m)];
        xlim([0 6.5]);
            ylim([-0.01 .15]);
            
%          f = figure; 
%         set(f,'renderer','painters');
%         shadedErrorBar(time, mean(totavgTraces,2)-mean(totavgTraces(1,:),2), std(totavgTraces,1,2),{'-'});
%         xlim([0 6.5]);
%             ylim([-0.01 .15]);
    end
    
end