%basis params
mdP = ceil(p.totL / 2);
rng = -floor(p.kML/2):floor(p.kML/2);
padL = floor(0.1 * p.kML);
kL = p.kML + 2 * padL;

hL = floor(kL / 2);
rngP = repmat(-hL:hL, p.nKernels);

%stereo basis initialization
phi = zeros(p.totL, p.nKernels, 2);
phi(mdP+rng,:,:) = 0.1 * randn(p.kML, p.nKernels,2);
%phi = normalize_matrix(phi);

%kernel lengths
phiLP = repmat(kL, [p.nKernels, 2]);

%TODO: change mask from binary inds to array inds
padMask = zeros(p.totL, p.nKernels, 2);
padMask(mdP + rngP, :, :) = 1;
padMask(mdP + rng, :, :) = 0;

%learning history
lHist = zeros(p.niter, 2);
nHist = zeros(p.niter, 2);
sHist = zeros(p.niter, 2);
snrHist = zeros(p.niter, 2);

tData = wavread(p.dataPath);
% load(p.dataPath);
% tData = tData(1:16000 * 5, :);
% if size(tData,2) < 2
%     error('Non stereo sound');
% end
tData = tData';

for i = 1:p.niter
    
    %alternate between left and right channels
    for chn = 1:2
                
        %E-STEP
        [w r] = myTemporalMP(tData(chn,:), phi(:,:,chn), p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.nsnr(i), p.mp.spkTrs(i));
        w = mapCoeff(tData(chn,:), phi(:,:,chn), w, p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
        
        %M-STEP
        [fval new_phi] = phi_cd(phi(:,:,chn), w, tData(chn,:), mdP, phiLP(:,chn));
        phi(:,:,chn) = new_phi;

        lHist(i, chn) = fval;
        nHist(i, chn) = nnz(w);
        sHist(i, chn) = sum(w(:));
        snrHist(i, chn) = snr(tData(chn,:), tData(chn,:) - reconstructSignal(w', phi(:,:,chn)));

        %normalize
        phi(:,:,chn) = normalize_matrix(phi(:,:,chn));

        %adjust kernel lengths
        [padMask(:,:,chn) phiLP(:,chn)] = adjustPadding(phi(:,:,chn), padMask(:,:,chn), phiLP(:,chn), p.exTrs);

        fprintf('\n\nITERATION %d CHANNEL: %d FVAL: %.5f \n', i, chn, fval);
        
    end
    
    %save
    if ~mod(i, saveEvery)
       fname = sprintf('CDspkTrs%.2f_%d_iter_STEREO.mat', p.mp.spkTrs(i), i);
       save(fname);
       fprintf('%s saved\n', fname);
    end
    
end
