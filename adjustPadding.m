function [nPadMask, nPhiLP] = adjustPadding(phi, padMask, phiLP, exTrs)

extL = 2;
nPadMask = padMask;
nPhiLP = phiLP;

for k = 1:size(phi,2)
    edges = phi(padMask(:,k) == 1, k);
    edgeSum = sum(abs(edges));
    trs = length(edges) * exTrs;
    
    if edgeSum > trs && phiLP(k) + 2 * extL <= size(phi,1)
        
        str = find(padMask(:,k), 1, 'first');
        lst = find(padMask(:,k), 1, 'last');
        
        inds = [str-extL:str-1 lst+1:lst+extL];
        inds = inds(inds >= 1 & inds <= size(phi,1));
        
        nPadMask(inds, k) = 1;
        nPhiLP(k) = nPhiLP(k) + 2 * extL;
    
%     elseif edgeSum < trs && phiLP(k) - 2 * extL > 0
%         
%         str = find(padMask(:,k), 1, 'first');
%         lst = find(padMask(:,k), 1, 'last');
%         
%         inds = [str:str+extL-1 lst-extL+1:lst];
%         inds = inds(inds >= 1 & inds <= size(phi,1));
%                 
%         nPadMask(inds, k) = 0;
%         nPhiLP(k) = nPhiLP(k) - 2 * extL;
    
    end
    
end