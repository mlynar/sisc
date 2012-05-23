%TODO:
%   -parametrization into one object
%   -sparse objective into separate function (in pimping.m) X
%   -kernel shrinkage
%   -likelihood history, convergence
%   -data sampling
%   -modularize code

%basis params
totL = 513;
mdP = (totL - 1) / 2;
nKernels = 32;

kML = 101;
rng = -floor(kML/2):floor(kML/2);
padL = floor(0.1 * kML);
kL = kML + 2 * padL;

hL = floor(kL / 2);
rngP = repmat(-hL:hL, nKernels);

%basis initialization
phi = zeros(totL, nKernels);
phi(mdP+rng,:) = randn(kML, nKernels);
phi = normalize_matrix(phi);

phiLP = repmat(kL, nKernels, 1);

%TODO: change mask from binary inds to array inds
padMask = zeros(totL, nKernels);
padMask(mdP + rngP, :) = 1;
padMask(mdP + rng, :) = 0;

%learning params
niter = 2000;
lr = 0.01;
exTrs = 0.05;
noiseVar = 1;

nbatch = ceil(linspace(1, 20, niter));

%learning history
lHist = zeros(1, niter);

%MP params
nonneg = false;
mpIter = 5000;
nsnr = ceil(linspace(5, 25, niter));

%testing toy dataset params
fcoefs = MakeERBFilters(6000, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,256)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);

for i = 1:niter
    fprintf('ITERATION %d\n', i);
        
    dPhi = zeros(totL, nKernels);    
    

    for j = 1:nbatch(i)
        tData = sampleToyData(gammaTones, 10000, 3);
        [w r] = myTemporalMP(tData, phi, nonneg, mpIter, nsnr(i));

        lHist(i) = lHist(i) + sum(sqrt((tData - r').^2));
        
        %TODO: why converges faster when residue = r?
        residue = tData'- r;

        for k = 1:nKernels
            tPs = find(w(:,k) ~= 0);

            hLK = floor(phiLP(k) / 2);
            rngK = -hLK:hLK;

            for tp = 1:length(tPs)
                t = tPs(tp);         

                tInds = t + rngK;
                inds = tInds > 0 & tInds < size(w,1);
                tInds = tInds(inds);
                kInds = mdP + rngK(inds);

                dPhi(kInds, k) = dPhi(kInds, k) + w(t,k) * residue(tInds);
            end
        end
        
    end

    phi = phi + lr * dPhi;
    phi = normalize_matrix(phi);
    
    %setup kernel lengths TODO: add kernel shrinkage
    [padMask phiLP] = adjustPadding(phi, padMask, phiLP, exTrs);
    
    lHist(i) = lHist(i) / nbatch(i);
    
%     if ~mod(i, 10)
%         close all;
%         plotbs(phi,0);
%         pause(0.3);
%     end
end

save('TTT.mat');