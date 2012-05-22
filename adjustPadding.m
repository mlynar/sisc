function [nPadMask, nPhiLP] = adjustPadding(phi, padMask, phiLP, exTrs)
    nKernels = size(phi,2);
    totL = size(phi,1);
    
    nPadMask = padMask;
    nPhiLP = phiLP;
    
    relL = 0.1;
    
    shrTrs = 0.5 * exTrs;
    
    for i = 1:nKernels        
        padInds = find(padMask(:,i));
        pad = phi(padInds, i);
        
        lpI = length(padInds) / 2;
        kS = padInds(lpI);
        kE = padInds(lpI + 1);
        endL = floor(relL * (kE - kS + 1));
       
        nds = phi([kS+1:kS+endL kE-endL:kE-1], i);
           
        mFirst = find(padMask(:,i), 1, 'first');
        mLast = find(padMask(:,i), 1, 'last');
        
        if sum(abs(pad)) > exTrs * sum(padMask(:,i))
            nPadL = floor(relL * phiLP(i)); 
            
            inds = [(mFirst-nPadL:mFirst-1) (mLast+1:mLast+nPadL)];
            inds = inds(inds > 0 & inds < totL);
            
            nPadMask(:,i) = 0;
            nPadMask(inds, i) = 1;
            nPhiLP(i) = phiLP(i) + 2 * nPadL;

         %TODO: add kernel shrinkage, remove arbitrary value
         %kernel shrinkage from both sides
         elseif  sum(abs(nds)) < shrTrs * 2 * endL & phiLP(i) > 91
             nPadL = max(2, floor(relL * phiLP(i))); 
             
             inds = [(mFirst:mFirst+nPadL-1) (mLast-nPadL+1:mLast)];
             inds = inds(inds > 0 & inds < totL);
             
             nPadMask(:,i) = 0;
             nPadMask(inds, i) = 1;
             
             nl = phiLP(i) - 2 * (nPadL - 1);           
             nPhiLP(i) = nl;
        end        
    end