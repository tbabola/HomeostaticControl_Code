function [ img ] = bfLoadTif( filename )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
        bf = bfopen([filename]);
        [m,n] = size(bf{1}{1});
        t = size(bf{1},1);
        img = zeros(m,n,t,'int16');
        for i=1:t
            img(:,:,i) = bf{1}{i};
        end
end

