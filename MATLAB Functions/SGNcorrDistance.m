dname = '.\Data\In Vitro Experiments';
list = loadFileList([dname '\WT + NBQX\*\*V3.mat']);
csvlist = loadFileList([dname '\WT + NBQX\*\*CM.csv']);


for i = 1:size(csvlist,1)
    load(list{i});
    keepTrace = int8(find([eventStruct.SGN(:).keepTrace]));

    SGNsignals = [eventStruct.SGN(keepTrace).dFoF];
    SGNsignals = SGNsignals(1:600,:);
    corrMat = corrcoef(SGNsignals);
    
    csvR = csvread(csvlist{i},1,1);
    csvR = csvR(keepTrace',:);
    pts = csvR;
    
    distances = [];
    for j=1:size(csvR,1)
        pt = csvR(j,:); 
        distances(:,j) = sqrt((pts(:,1)-pt(:,1)).^2+(pts(:,2)-pt(:,2)).^2);
    end
    
    [b, ordered] = sort(distances,'ascend');
    
    meanNeighbors = 5;
    close = [];
    far = [];
    for j=1:size(csvR,1)
        close(:,j) = corrMat(j,ordered(2:meanNeighbors+1,j));
        far(:,j) = corrMat(j,ordered(end-meanNeighbors:end,j));
    end
    
    WTcloseGroup(i) = mean(mean(close));
    WTfarGroup(i) = mean(mean(far));
    
end

%% dname = 'M:\Bergles Lab Data\Projects\In Vivo Imaging and VG3 Paper\In Vitro Experiments';
list = loadFileList([dname '\VG3KO + NBQX + BBE\*\*V3.mat']);
csvlist = loadFileList([dname '\VG3KO + NBQX + BBE\*\*CM.csv']);


for i = 1:size(csvlist,1)
    load(list{i});
    keepTrace = int8(find([eventStruct.SGN(:).keepTrace]));

    SGNsignals = [eventStruct.SGN(keepTrace).dFoF];
    SGNsignals = SGNsignals(1:600,:);
    corrMat = corrcoef(SGNsignals);
    
    csvR = csvread(csvlist{i},1,1);
    csvR = csvR(keepTrace',:);
    pts = csvR;
    
    distances = [];
    for j=1:size(csvR,1)
        pt = csvR(j,:); 
        distances(:,j) = sqrt((pts(:,1)-pt(:,1)).^2+(pts(:,2)-pt(:,2)).^2);
    end
    
    [b, ordered] = sort(distances,'ascend');
    
    meanNeighbors = 5;
    close = [];
    far = [];
    for j=1:size(csvR,1)
        close(:,j) = corrMat(j,ordered(2:meanNeighbors+1,j));
        far(:,j) = corrMat(j,ordered(end-meanNeighbors:end,j));
    end
    
    VG3closeGroup(i) = mean(mean(close));
    VG3farGroup(i) = mean(mean(far));
    
end

%% plot
figure;
figure;
meanWTclose = mean(WTcloseGroup);
sterrWTclose = std(WTcloseGroup,[],2)/sqrt(size(WTcloseGroup,2));
meanWTfar = mean(WTfarGroup);
sterrWTfar = std(WTfarGroup,[],2)/sqrt(size(WTfarGroup,2));
meanVG3close = mean(VG3closeGroup);
sterrVG3close = std(VG3closeGroup,[],2)/sqrt(size(VG3closeGroup,2));
meanVG3far = mean(VG3farGroup);
sterrVG3far = std(VG3farGroup,[],2)/sqrt(size(VG3farGroup,2));
WTn = size(WTcloseGroup,2);
VG3n = size(VG3closeGroup,2);
errorbar(1, meanWTclose, sterrWTclose,'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10); hold on;
errorbar(1.5, meanWTfar, sterrWTfar,'LineStyle', 'none','Color','k','CapSize',0,'Marker','.','MarkerSize',10);
errorbar(2, meanVG3close, sterrVG3close,'LineStyle', 'none','Color','r','CapSize',0,'Marker','.','MarkerSize',10);
errorbar(2.5, meanVG3far, sterrVG3far,'LineStyle', 'none','Color','r','CapSize',0,'Marker','.','MarkerSize',10);
line([1*ones(1,WTn); 1.5*ones(1,WTn)], [WTcloseGroup;WTfarGroup],'Color',[0.7 0.7 0.7]);
line([2*ones(1,VG3n); 2.5*ones(1,VG3n)], [VG3closeGroup;VG3farGroup],'Color',[0.7 0.7 0.7]);
xlim([.5 3]);
ylim([0 1]);
