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
mdP = ceil(p.totL / 2);
rng = -floor(p.kML/2):floor(p.kML/2);
padL = floor(0.1 * p.kML);
kL = p.kML + 2 * padL;

hL = floor(kL / 2);
rngP = repmat(-hL:hL, p.nKernels);

%basis initialization
phi = zeros(p.totL, p.nKernels);
phi(mdP+rng,:) = 0.1 * randn(p.kML, p.nKernels);
%phi = normalize_matrix(phi);

phiLP = repmat(kL, p.nKernels, 1);

%TODO: change mask from binary inds to array inds
padMask = zeros(p.totL, p.nKernels);
padMask(mdP + rngP, :) = 1;
padMask(mdP + rng, :) = 0;

%learning history
lHist = zeros(1, p.niter);
nHist = zeros(1, p.niter);
sHist = zeros(1, p.niter);

%testing toy dataset params
fcoefs = MakeERBFilters(6000, 32, 200);
gammaTones = ERBFilterBank([1 zeros(1,256)], fcoefs)';
gammaTones = normalize_matrix(gammaTones);
freqs = ERBSpace(200, 6000/2, 32);
tData = sampleToyData(gammaTones(:, 1:32), 6000 * 20, 3);
%load('data/tData2');
%tData = tData(1:10*6000);

for i = 1:p.niter
    %sample
    %tData = sampleToyData(gammaTones(:,32), 10000, 3);   
    %tData = loadChunk(p.dataPath)';
    
    %matching pursuit
    [w r] = myTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.nsnr(i), p.mp.spkTrs(i));
    
    %MAP
    %w = mapCoeff(tData, phi, w, p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
    
    %ml gradient step
    [fval dPhi] = phiObjective(phi, w, tData, mdP, phiLP);
    
    phiN = phi + p.lr * dPhi;
    phiN = normalize_matrix(phiN);
    
    fvalN = phiObjective(phiN, w, tData, mdP, phiLP);
    
    % pursue a constant change in angle
    angle_phi = acos(phiN(:)' * phi(:) / sqrt(sum(phiN(:).^2)) / sqrt(sum(phi(:).^2)));
    if angle_phi < p.target_angle
          p.lr = p.lr * 1.01;
    else
          p.lr = p.lr * 0.99;
    end

    %change only if increased
    if fvalN < fval
        fprintf('WARNING: gradient increase\n');
        p.lr = p.lr * 0.5;
    else
        phi = phiN;
    end
    
    lHist(i) = fval;
    nHist(i) = nnz(w);
    sHist(i) = sum(w(:));
    
    
    %normalize
%    phi = normalize_matrix(phi);
    
    %setup kernel lengths TODO: add kernel shrinkage should it be here or
    %two line later?
    %[padMask phiLP] = adjustPadding(phi, padMask, phiLP, p.exTrs);
    
    fprintf('\n\nITERATION %d LR: %.5f FVAL: %.5f \n', i, p.lr, fval);
    
    %save
    if ~mod(i, saveEvery)
       fname = sprintf('NspkTrs%.2f_%d_iter.mat', p.mp.spkTrs(i), i);
       save(fname);
       fprintf('%s saved\n', fname);
    end
        
    
end
