%TODO:
%   -parametrization into one object X
%   -sparse objective into separate function (in pimping.m) X
%   -kernel shrinkage
%   -likelihood history, convergence X
%   -data sampling
%   -modularize code

%basis params
mdP = (p.totL - 1) / 2;
rng = -floor(p.kML/2):floor(p.kML/2);
padL = floor(0.1 * p.kML);
kL = p.kML + 2 * padL;

hL = floor(kL / 2);
rngP = repmat(-hL:hL, p.nKernels);

%basis initialization
phi = zeros(p.totL, p.nKernels);
phi(mdP+rng,:) = randn(p.kML, p.nKernels);
phi = normalize_matrix(phi);

phiLP = repmat(kL, p.nKernels, 1);

%TODO: change mask from binary inds to array inds
padMask = zeros(p.totL, p.nKernels);
padMask(mdP + rngP, :) = 1;
padMask(mdP + rng, :) = 0;

%learning history
lHist = zeros(1, p.niter);

%testing toy dataset params
fcoefs = MakeERBFilters(6000, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,256)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);


for i = 1:p.niter
    fprintf('ITERATION %d\n', i);
        
    dPhi = zeros(p.totL, p.nKernels);    
    
    for j = 1:p.nbatch(i)
        tData = sampleToyData(gammaTones(:,28:32), 3000, 3);
                
        [w r] = myTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16);%, 0, p.mp.nsnr(i));
        
        wp = mapCoeff(tData, phi, w, p.map.mapIter);
        r = reconstructSignal(wp', phi)';
        
        %TODO: why converges faster when residue = r?
        residue = tData - r';
        residue = [zeros(1, p.totL) residue zeros(1, p.totL)];

        lHist(i) = lHist(i) + norm(residue, 2);
            
        for k = 1:p.nKernels
            tPs = find(w(:,k) ~= 0);

            hLK = floor(phiLP(k) / 2);
            rngK = -hLK:hLK;
            
            for tp = 1:length(tPs)
                t = tPs(tp);         
                tInds = t + p.totL + rngK;
                kInds = mdP + rngK;
                
                dPhi(kInds, k) = dPhi(kInds, k) + w(t,k) * residue(tInds)';
            end            
            
        end
    end

    dPhi = dPhi / p.nbatch(i);
    phi = phi + p.lr * dPhi;
    
    phi = normalize_matrix(phi);
    
    %setup kernel lengths TODO: add kernel shrinkage
    %[padMask phiLP] = adjustPadding(phi, padMask, phiLP, p.exTrs);
    
    lHist(i) = lHist(i) / p.nbatch(i);
       
     if ~mod(i, 10)
         close all;
         plotbs(phi,0);
         pause(0.3);
     end
end
