%TODO:
%   -parametrization into one object X
%   -sparse objective into separate function (in pimping.m) X
%   -kernel shrinkage X
%   -likelihood history, convergence X
%   -batch learning (?) X
%   -data sampling X
%   -modularize code X
%   -implement MAP fine tuning with a minimize() function
%   -check why it fails with short total kernel lengths X
%   -change trimming to both sides

%basis params
mdP = ceil(p.totL / 2);
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
nHist = zeros(1, p.niter);
sHist = zeros(1, p.niter);
snrHist = zeros(1, p.niter);
normHist = zeros(1, p.niter);

load(p.dataPath);
tData = tData(1:16000 * 30, 2)';

normPrev = norm(phi(:), 2);
for i = 1:p.niter
    
    %print
    if ~mod(i, p.printEvery)
        %convert -delay 20 -loop 200 phi_it*.png animation.gif
        if i < 10
            fname = sprintf('phi_it0%d.png', i);
        else
            fname = sprintf('phi_it%d.png', i);
        end
        drawKernels(phi, padMask, phiLP, fname);
        pause(0.3);
        close all;
    end
    
    %sample
    %tData = loadChunk(p.dataPath)';
    
    %E-STEP
    [w r] = myTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.nsnr(i), p.mp.spkTrs(i));
    if ~isreal(w)
        error('Imaginary sparse matrix\n');
    end
    
    %MAP
    %w = mapCoeff(tData, phi, w, p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
    %[w f] = minimize(w(:), 's_objective_mp', 10, tData, phi, p.map.lambda);
    %w = reshape(w, length(tData), p.nKernels);
    
    %M-STEP
    [fval new_phi] = phi_cd(phi, w, tData, mdP, phiLP);
    phi = new_phi;
    
    lHist(i) = fval;
    nHist(i) = nnz(w);
    sHist(i) = sum(w(:));
    snrHist(i) = snr(tData, tData - reconstructSignal(w', phi));
    
    normCurr = norm(phi(:), 2);
    normHist(i) = normPrev - normCurr;
    normPrev = normCurr;
    
    %normalize
    phi = normalize_matrix(phi);
    
    %adjust kernel lengths
    [padMask phiLP] = adjustPadding(phi, padMask, phiLP, p.exTrs);
    
    fprintf('\n\nITERATION %d FVAL: %.5f DNORM: %.5f\n', i, fval, normHist(i));
    
    %save
    if ~mod(i, saveEvery)
       fname = sprintf('CDspkTrs%.2f_%d_iter.mat', p.mp.spkTrs(i), i);
       save(fname);
       fprintf('%s saved\n', fname);
    end
    
end
