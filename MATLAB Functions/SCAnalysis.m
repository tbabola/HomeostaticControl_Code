file = '.\Data\WILDTYPE\Experiment-265\Experiment-265-SC.tif'
img = loadTif(file,16);
dFoF = normalizeImg(img,10,1);

%%
[LSC, RSC] = ROIselectionSC(dFoF);
LSCmean = mean(LSC);
LSCstd = std(LSC,[],1);
thr =  LSCmean + 2*LSCstd;
LSCbin = LSC > thr;
imagesc(LSCbin);

RSCmean = mean(RSC);
RSCstd = std(RSC,[],1);
thr =  RSCmean + 2*RSCstd;
RSCbin = RSC > thr;
imagesc(RSCbin);


function [LSC, RSC] = ROIselectionSC(movie)
    show = movie(:,:,1);
    [m,n,t] = size(movie);
    movieRsp = reshape(movie,m*n,t);
    figure;
    h_im = imagesc(show);
    caxis([min(show(:)) max(show(:))]);
    %rotate images
    
    RSC = [];
    RSCpoly = impoly(gca,[255,33;255,168;400,168;350,75]);
    wait(RSCpoly);
    RSCmask = createMask(RSCpoly, h_im);
    [r] = find(RSCmask);
    RSC = movieRsp(r,:)';
    
    
     LSC = impoly(gca,[190,33;190,168;45,168;95,75]);
     wait(LSC);
     LSCmask = createMask(LSC, h_im);
     [r] = find(LSCmask);
     LSC = movieRsp(r,:)';
end