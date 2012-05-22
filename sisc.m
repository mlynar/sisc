%TODO:
%   -parametrization into one object
%   -sparse objective into separate function (in pimping.m)
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
niter = 300;
nMapIter = 300;
mapIter = 100;
lr = 0.05;
aTrs = 1e-5;
exTrs = 0.05;
noiseVar = 1;
nbatch = ceil(linspace(1, 50, niter));

%learning history
lHist = zeros(1, niter);

%MP params
nonneg = false;
mpIter = 2000;
eps =  1e-4;

%testing toy dataset params
fcoefs = MakeERBFilters(6000, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,100)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);

for i = 70:niter
    fprintf('ITERATION %d\n', i);
        
    dPhi = zeros(totL, nKernels);    
    
    for jB = 1:nbatch(i);
        tData = sampleToyData(gammaTones, 12000, 12001); %TODO: write real version

        %[w r] = temporalMP(tData, phi, nonneg, mpIter, eps);
        [w r] = myTemporalMP(tData, phi, nonneg, mpIter, 0.1);

        lHist(i) = lHist(i) + sum(sqrt((tData - r').^2));
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
    end

    phi = phi + lr * dPhi;
    phi = normalize_matrix(phi);
    
    %setup kernel lengths TODO: add kernel shrinkage
    [padMask phiLP] = adjustPadding(phi, padMask, phiLP, exTrs);
    
    lHist(i) = lHist(i) / nbatch(i);
    
    %if ~mod(i, 5)
    %    close all;
    %    plotbs(phi,0);
        %figure; plot(lHist, 'r','LineWidth',2); grid()
    %    pause(0.3);
    %end
end

%save('LEARNED_KERNELS.mat');