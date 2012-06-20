function [xc rng] = mXcorr(A, B)

ta = size(A, 2);
tb = size(B, 2);

if ta ~= tb
    error('different time lengths');
end

rng = -ta+1:ta-1;

xc = zeros(size(rng));
for i = 1:length(rng)
    del = rng(i);
    
    %C = zeros(size(B));
    
    if del <= 0
        xc (i) = sum(sum(B(:,-del+1:ta) .* A(:,1:ta+del)));
        %C(:,1:ta+del) = B(:,-del+1:ta);
    else
        xc(i) = sum(sum((B(:,1:ta-del) .* A(:,1+del:ta))));
        %C(:,1+del:ta) = B(:,1:ta-del);
    end   
    
    %xc(i) = corr(A(:), C(:));
    %xc(i) = C(:)' * D(:);
end

