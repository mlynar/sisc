function [fval new_phi] = phi_cd(phi, w, tData, mdP, phiLP)
%TODO: incorporate noise variance write EPSILON
%fval should be maximized
%and dPhi is a gradient of log posterior (a function to be maximized)

    rec = reconstructSignal(w', phi);
    res = tData - rec;
        
    fval = - sum(res.^2) / 2;
    
    eps = 1e-8;
    
    if nargout > 1
        nKernels = size(phi, 2);
        totL = size(phi, 1);

        new_phi = phi;
        tData = [zeros(1, totL) tData zeros(1, totL)];
         
         for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);

            if ~isempty(tPs)                

                hLK = floor(phiLP(k) / 2);
                rngK = -hLK:hLK;
                kInds = mdP + rngK;

                phi_m = zeros(size(phi,1),1);
                
                %new version
                resK = tData - [zeros(1, totL) reconstructSignal(w', new_phi, k) zeros(1, totL)];
                %end new version

                for tp = 1:length(tPs)
                    t = tPs(tp);         
                    tInds = t + totL + rngK;    

                    phi_m(kInds) = phi_m(kInds) + resK(tInds)' / w(t,k);                
                end
                
                new_phi(:, k) = phi_m ./ length(tPs); 
            end    
         end        
    end
            
        