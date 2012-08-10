function [fval new_phi] = phi_ndim_cd(phi, w, tData, mdP, phiLP)
%TODO: 
%   -incorporate noise variance write EPSILON
%   -generalize to all numbers of channels (from one to N);
%       
%fval should be maximized
%and dPhi is a gradient of log posterior (a function to be maximized)

    rec = reconstructNdimSignal(w', phi);
    res = tData - rec;
    res = res(:);
    
    %fval = -sum(res.^2) / 2;
    fval = -sum(res.^2) / 2;
    
    eps = 1e-8;
    
    if nargout > 1
        nKernels = size(phi,3);
        totL = size(phi, 2);
        nChn = size(phi, 1);

        new_phi = phi;
        tData = [zeros(nChn, totL) tData zeros(nChn, totL)];
         
         for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);

            if ~isempty(tPs)                

                hLK = floor(phiLP(k) / 2);
                rngK = -hLK:hLK;
                kInds = mdP + rngK;

                phi_m = zeros(nChn, totL); 
                
                %new version
                resK = tData - [zeros(nChn, totL) reconstructNdimSignal(w', new_phi, k) zeros(nChn, totL)];
                %end new version

                for tp = 1:length(tPs)
                    t = tPs(tp);         
                    tInds = t + totL + rngK;    

                    phi_m(:,kInds) = phi_m(:,kInds) + resK(:,tInds) / w(t,k);                
                end
                
                new_phi(:, :, k) = phi_m ./ length(tPs); 
            end    
         end        
    end