function data = corticalPeaksQuant(paths)
     ctxEvents = [];
     m = size(paths,1);
     for i = 1:m
         load(paths{i});
         disp(paths{i});
         [filePath name ext] = fileparts(paths{i});
         
            time = [1:1:size(ICsignal,1)]';
            LIC = smooth(ICsignal(:,1));
            RIC = smooth(ICsignal(:,2));
            ctx = smooth(ICsignal(:,3));
            LIC = msbackadj(time,LIC,'WindowSize',500,'StepSize',500);
            RIC = msbackadj(time,RIC,'WindowSize',500,'StepSize',500);
            ctx = msbackadj(time,ctx,'WindowSize',500,'StepSize',500);

            %parameters
            pkThreshold = .02;
            pkMinHeight = .02;
            pkDistance = 5; %in frame, 10 = 1s
            [pks,locs,w] = findpeaks(LIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkMinHeight,'MinPeakDistance',pkDistance,'Annotate','extents');
            [LIC_itk L_itkill] = ctx_PkRemoval(ctx, [LIC RIC], locs);
            sizeL = size(pks,1) - size(LIC_itk,2);
            locsLIC = locs(LIC_itk);
            killLocs = locs(L_itkill);



            [pks,locs,w] = findpeaks(RIC,'MinPeakProminence',pkThreshold,'MinPeakHeight',pkMinHeight,'MinPeakDistance',pkDistance,'Annotate','extents');
            [RIC_itk R_itkill] = ctx_PkRemoval(ctx, [LIC RIC], locs);
            sizeR = size(pks,1) - size(RIC_itk,2)
            locsRIC = locs(RIC_itk);

            data(i).time = time;
            data(i).signals = [LIC RIC ctx];
            data(i).ctxEvents = (sizeL+sizeR)/2;
            data(i).killLocs = killLocs;
            
     end
   
end
 
function [indexToKeep, indexToKill] = ctx_PkRemoval(ctxSignal, ICsignal, IClocs)
    time = [1:1:size(ctxSignal,1)]';
    
    %[pks,locs,w] = findpeaks(msbackadj(time,ctxSignal),'MinPeakHeight',0.03,'Annotate','extents');
    
    ctxBright = ctxSignal > ICsignal(:,1)*1.2 & ctxSignal > ICsignal(:,2)*1.2 & ctxSignal > 0.01;
    
    indexToKeep = [];
    indexToKill = [];
    for i=1:size(IClocs)
        if ctxBright(IClocs(i)) == 0
            indexToKeep = [indexToKeep i];
        else
            indexToKill = [indexToKill i];
        end      
    end
end