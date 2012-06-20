function [fval new_phi] = phi_mstep(phi, w, tData, mdP, phiLP)
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

        new_phi = zeros(size(phi));
        tData = [zeros(1, totL) tData zeros(1, totL)];
         
         for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);

            hLK = floor(phiLP(k) / 2);
            rngK = -hLK:hLK;
            kInds = mdP + rngK;
            
            if ~isempty(tPs)
                %new version
                resK = tData - [zeros(1, totL) reconstructSignal(w', phi, k) zeros(1, totL)];
                %end new version

                for tp = 1:length(tPs)
                    t = tPs(tp);         
                    tInds = t + totL + rngK;    
    %                 new_phi(kInds, k) = new_phi(kInds, k) + (1 / w(t,k)) * tData(tInds)';

                    new_phi(kInds, k) = new_phi(kInds, k) + resK(tInds)' / w(t,k);                
                end
                
                new_phi(kInds, k) = new_phi(kInds, k) ./ length(tPs); 
                
            else
                new_phi(kInds, k) = phi(kInds, k);
            end
            
         end        
    end
            
            %new version
%             sC = sum(w(tPs, k));
%             if sC ~= 0
%                 new_phi(:, k) = new_phi(:, k) ./ sC;
%             end
            %end new version
        