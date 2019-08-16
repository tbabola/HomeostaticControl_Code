function [d, time] = abfloadClean(fn)
%abfloadClean loads abf and cleans up file
    [d,si] = abfload(fn);
    d = squeeze(d);
    
    %remove clampex "holding period" that is 1/64th of trace
    pntsToRemove = round(size(d,1)/64);
    d = d(pntsToRemove:end,:);
    
    time = zeros(size(d));
    time = 0:si*10^-3:(size(time,1)-1)*si*10^-3;
    time = time';

end

