function [avg, se, fg] = LRpie(stats)
    lt_org = [255, 166 , 38]/255;
    dk_org = [255, 120, 0]/255;
    lt_blue = [50, 175, 242]/255;
    dk_blue = [0, 13, 242]/255;
 
    avg = mean(stats.Var6,1);
    se = std(stats.Var6,1,1);%/sqrt(size(stats.Var6,1));
    fg = figure;
    fg.Units = 'inches';
    fg.Position = [5 5 1 1];
    h = pie([avg 1-avg]);
    hp = findobj(h, 'Type', 'patch');
    hp(1).FaceColor = lt_org;
    hp(1).LineStyle = 'none';
    hp(2).FaceColor = lt_blue;
    hp(2).LineStyle = 'none';
    
end