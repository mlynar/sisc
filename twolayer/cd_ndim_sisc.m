%basis params
mdP = ceil(p.totL / 2);
rng = -floor(p.kML/2):floor(p.kML/2);
padL = floor(0.1 * p.kML);
kL = p.kML + 2 * padL;

hL = floor(kL / 2);
rngP = repmat(-hL:hL, p.nKernels);

%basis initialization
phi = zeros(p.nCh, p.totL, p.nKernels);
phi(:, mdP+rng, :) = randn(p.nCh, p.kML, p.nKernels);
phi = normalize_ndim_matrix(phi);

phiLP = repmat(kL, p.nKernels, 1);

%TODO: change mask from binary inds to array inds
padMask = zeros(p.totL, p.nKernels);
padMask(mdP + rngP, :) = 1;
padMask(mdP + rng, :) = 0;

%learning history
lHist = zeros(1, p.niter);
nHist = zeros(1, p.niter);
sHist = zeros(1, p.niter);
snrHist = zeros(1, p.niter);
normHist = zeros(1, p.niter);

%load(p.dataPath);
%tData = tData(1:16000 * 30, 2)';

normPrev = norm(phi(:), 2);

for i = 1:p.niter
    
        
    %E-STEP
    [w r] = multiTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, p.mp.nsnr(i), p.mp.spkTrs(i));
        
    %MAP
    w = map_ndim_coeff(tData, phi, w', p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
    
    %M-STEP
    [fval new_phi] = phi_ndim_cd(phi, w, tData, mdP, phiLP);
    phi = new_phi;
    
    lHist(i) = fval;
    nHist(i) = nnz(w);
    sHist(i) = sum(w(:));
    noise = tData - reconstructNdimSignal(w', phi);
    snrHist(i) = snr(tData(:), noise(:));
    
    normCurr = norm(phi(:), 2);
    normHist(i) = normPrev - normCurr;
    normPrev = normCurr;
    
    %normalize
    phi = normalize_ndim_matrix(phi);
    
    %adjust kernel lengths
    %[padMask phiLP] = adjustPadding(phi, padMask, phiLP, p.exTrs);
    
    fprintf('\n\nITERATION %d FVAL: %.5f DNORM: %.5f\n', i, fval, normHist(i));
    
    %save
    if ~mod(i, saveEvery)
       fname = sprintf('CDspkTrs%.2f_%d_iter.mat', p.mp.spkTrs(i), i);
       save(fname);
       fprintf('%s saved\n', fname);
    end
    
end
