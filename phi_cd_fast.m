function [fval new_phi] = phi_cd_fast(phi, w, tData, mdP, phiLP)
%TODO: incorporate noise variance write EPSILON
%      reconstruct only user parts of signal, not entire one
%fval should be maximized
%and dPhi is a gradient of log posterior (a function to be maximized)

    rec = reconstructSignal(w', phi);
    res = tData - rec;
        
    fval = - sum(res.^2) / 2;
    
    eps = 1e-8;
    
    if nargout > 1
        nKernels = size(phi, 2);
        totL = size(phi, 1);
        datL = length(w);

        new_phi = phi;
        tData = [zeros(1, totL) tData zeros(1, totL)];
        
        xHat = [zeros(1, totL) reconstructSignal(w', new_phi) zeros(1, totL)];
        resK = tData - xHat;
        
         for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);

            if ~isempty(tPs)                

                hLK = floor(phiLP(k) / 2);
                rngK = -hLK:hLK;
                kInds = mdP + rngK;

                %new version
                xHatK = convnfft(w(:,k)', new_phi(:,k)');
                xHatK = xHatK(mdP + (0:datL-1));
                %xHat(totL + (1:datL)) = xHat(totL + (1:datL)) - xHatK;
                resK(totL + (1:datL)) = resK(totL + (1:datL)) + xHatK;
                                
                %resK = tData - xHat;
                %end new version
                
                phi_m = zeros(size(phi,1),1);
                
                for tp = 1:length(tPs)
                    t = tPs(tp);         
                    tInds = t + totL + rngK;  

                    phi_m(kInds) = phi_m(kInds) + resK(tInds)' / w(t,k);                
                end
                
                new_phi(:, k) = phi_m ./ length(tPs); 
                
                xHatK = convnfft(w(:,k)', new_phi(:,k)');
                xHatK = xHatK(mdP + (0:datL-1));
                %xHat(totL + (1:datL)) = xHat(totL + (1:datL)) + xHatK;
                resK(totL + (1:datL)) = resK(totL + (1:datL)) - xHatK;
            end    
         end        
    end