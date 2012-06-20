%TODO:
%   -parametrization into one object X
%   -sparse objective into separate function (in pimping.m) X
%   -kernel shrinkage
%   -switching off inactive kernels during learning
%   -likelihood history, convergence X
%   -batch learning (?)
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
phi(mdP+rng,:) = -0.01 + 0.02 * rand(p.kML, p.nKernels);
%phi = normalize_matrix(phi);

phiLP = repmat(kL, p.nKernels, 1);

%TODO: change mask from binary inds to array inds
padMask = zeros(p.totL, p.nKernels);
padMask(mdP + rngP, :) = 1;
padMask(mdP + rng, :) = 0;

%learning history
lHist = zeros(1, p.niter);

%rprop params
d = repmat(d0, size(phi));
dP = repmat(d0, size(phi));
dPhiOld = zeros(size(dP));
% 
% fcoefs = MakeERBFilters(6000, 32, 200);
% gammaTones = ERBFilterBank([1 zeros(1,256)], fcoefs)';
% gammaTones = normalize_matrix(gammaTones);
% freqs = ERBSpace(200, 6000/2, 32);

% tData = sampleToyData(gammaTones(:, p.td.kf), 6000 * p.td.ns, p.td.spk);
load('data/tData2');
tData = tData(1:30*6000);

for i = 1:p.niter
    
    %tData = sampleToyData(gammaTones(:, p.td.kf), 6000 * p.td.ns, p.td.spk);
    
    %tData = loadChunk(p.dataPath)';
    
    %matching pursuit
    %[w r] = myTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.nsnr(i));
    [w r] = myTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.nsnr(i), p.mp.spkTrs(i));
    
    %w = mapCoeff(tData, phi, w, p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
    
    
    [fval dPhi] = phiObjDesc(phi, w, tData, mdP, phiLP);
    lHist(i) = fval;
    
    %rprop
    for k = 1:p.nKernels
        for l = 1:p.kML
            ind = mdP + rng(l);
                        
            if dPhiOld(ind, k) * dPhi(ind, k) > 0
                d(ind, k) = min(d(ind, k) * etaP, dMax);
                dP(ind, k) = -sign(dPhi(ind, k)) * d(ind, k);
                phi(ind, k) = phi(ind, k) + dP(ind, k);
            
            elseif dPhiOld(ind, k) * dPhi(ind, k) < 0
                d(ind, k) = max(d(ind, k) * etaM, dMin);
                phi(ind, k) = phi(ind, k) - dP(ind, k);
                dPhi(ind, k) = 0;
                
            elseif dPhiOld(ind, k) * dPhi(ind, k) == 0
                dP(ind, k) = -sign(dPhi(ind, k)) * d(ind, k);
                phi(ind, k) = phi(ind, k) + dP(ind, k);
                
            end
        end
    end
    
    %padding
    %[padMask phiLP] = adjustPadding(phi, padMask, phiLP, p.exTrs);
    
    %normalize
    phi = normalize_matrix(phi);    
    
    %rewrite
    dPhiOld = dPhi;
    
    %give output
    fprintf('\n\nITERATION %d FVAL: %.5f\n', i, fval);
    
    %save
    if ~mod(i, p.saveEvery)
       fname = sprintf('rprop_%d_iter.mat', i);
       save(fname);
       fprintf('%s saved\n', fname);
    end
    
end
