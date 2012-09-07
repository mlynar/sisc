%basis params
mdP = ceil(p.kML / 2);

%stereo basis initialization
phi = rand(2, p.kML, p.nKernels);
phi = normalize_ndim_matrix(phi);

%kernel lengths
phiLP = repmat(p.kML, p.nKernels, 1);

%learning history
lHist = zeros(p.niter);
nHist = zeros(p.niter);
sHist = zeros(p.niter);
snrHist = zeros(p.niter);

tData = wavread(p.dataPath)';
% load(p.dataPath);
% tData = tData(1:16000 * 120, :);
% if size(tData,2) < 2
%     error('Non stereo sound');
% end
% tData = tData';

for i = 1:p.niter
       
    %E-STEP
    [w r] = N_temporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, 0, p.mp.spkTrs(i));
    %multiTemporalMP(tData, phi, p.mp.nonneg, p.mp.mpIter(i), 1e-16, p.mp.nsnr(i), p.mp.spkTrs(i));
    %w = map_ndim_coeff(tData, phi, w', p.map.mapIter, p.map.lr, p.map.lambda, p.map.noiseVar, p.map.verbose);
    
    %M-STEP
    tic
    [fval phi] = phi_ndim_cd_fast(phi, w, tData, mdP, phiLP);
    toc
    phi = normalize_ndim_matrix(phi);
    
    lHist(i) = fval;
    nHist(i) = nnz(w);
    
    fprintf('ITERATION %d LOG L: %.5f\n\n', i, lHist(i));    
    
    
    %save
    if ~mod(i, saveEvery)
       %fname = sprintf('res/joint_stereo/CDspkTrs%.2f_%d_iter_STEREO_SI.mat', p.mp.spkTrs(i), i);
       fname = [p.savePath sprintf('jnt_st_180s_%d_kernels_%d_iter_STEREO_SI.mat', p.nKernels, i)];
       save(fname);
       fprintf('%s saved\n', fname);
    end
    
    %print
    if ~mode(i, p.printEvery)
        fname = sprintf('stereo_%d.png', i);
        drawStereoKernels(phi, fname);
    end
    
end
