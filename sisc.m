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
niter = 300;
nMapIter = 300;
mapIter = 100;
lr = 0.05;
aTrs = 1e-5;
exTrs = 0.05;
noiseVar = 1;

%learning history
lHist = zeros(1, niter);

%MP params
nonneg = false;
mpIter = 5000;
eps =  1e-4;

%testing toy dataset params
fcoefs = MakeERBFilters(6000, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,100)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);

for i = 1:niter
    fprintf('ITERATION %d\n', i);
        
    dPhi = zeros(totL, nKernels);    
    
    tData = sampleToyData(gammaTones, 3000, 20000); %TODO: write real version
    
    [w r] = temporalMP(tData, phi, nonneg, mpIter, eps);
    
    %lHist(i) = sum((r-tData').^2 / (noiseVar * length(tData)));
    %fprintf('Likelihood: %.5f\n', lHist(i));
    
    %MAP tuning; 
    if i >= nMapIter
        w = pimping(tData, phi, w, mapIter);
        r = reconstruct_signal(w, phi);
    end
        
    for k = 1:nKernels
        tPs = find(abs(w(:,k)) > aTrs);
        
        hLK = floor(phiLP(k) / 2);
        rngK = -hLK:hLK;
        
        for tp = 1:length(tPs)
            t = tPs(tp);         

            tInds = t + rngK;
            inds = tInds > 0 & tInds < size(w,1);
            tInds = tInds(inds);
            kInds = mdP + rngK(inds);

            dPhi(kInds, k) = dPhi(kInds, k) + w(t,k) * r(tInds);
        end
    end
    
    phi = phi + lr * dPhi;
    phi = normalize_matrix(phi);
    
    %setup kernel lengths TODO: add kernel shrinkage
    [padMask phiLP] = adjustPadding(phi, padMask, phiLP, exTrs);
    
    if ~mod(i, 1)
        close all;
        plotbs(phi,0);
        pause(0.3);
    end
end

%save('LEARNED_KERNELS.mat');