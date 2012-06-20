function [fval dPhi] = phiObjective(phi, w, tData, mdP, phiLP)
%TODO: incorporate noise variance
%fval should be maximized
%and dPhi is a gradient of log posterior (a function to be maximized)

    rec = reconstructSignal(w', phi);
    res = tData - rec;
        
    fval = - sum(res.^2) / 2;
    
    if nargout > 1
            nKernels = size(phi, 2);
            totL = size(phi, 1);
        
            dPhi = zeros(size(phi));
            res = [zeros(1, totL) res zeros(1, totL)];
            
         
         for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);
            
            hLK = floor(phiLP(k) / 2);
            rngK = -hLK:hLK;
            kInds = mdP + rngK;
            
            for tp = 1:length(tPs)
                t = tPs(tp);         
                tInds = t + totL + rngK;    
                dPhi(kInds, k) = dPhi(kInds, k) + w(t,k) * res(tInds)';
            end
            
        end
     end