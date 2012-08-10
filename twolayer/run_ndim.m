%clear all
close all
clc

addpath(genpath('~/Desktop/code'))
addpath(genpath('~/Desktop/sisc'))

%basis params
p.totL = 9;
p.nKernels = 3;
p.kML = 5;
p.nCh = 3;
p.exTrs = 0.05; %0.02;

%learning params
p.niter = 20;
p.noiseVar = 1;

%MAP params
p.map.mapIter = 200;
p.map.lr = 0.01;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(50000, 1, p.niter);
p.mp.nsnr = repmat(120, 1, p.niter);
p.mp.spkTrs = repmat(0.1, 1, p.niter);
    
%display params
showEvery = 100;
saveEvery = 100;
p.printEvery = 100;

%data params
%p.dataPath = '~/Desktop/sounds/training/stereo_sounds.mat';

%run
cd_ndim_sisc

%%
figure;
subplot(2,3,1);
imagesc(phi(:,5-3:5+2,1));colorbar(); colormap('Gray');
subplot(2,3,2);
imagesc(phi(:,5-3:5+2,2)); colorbar(); colormap('Gray');
subplot(2,3,3)
imagesc(phi(:,5-3:5+2,3)); colorbar(); colormap('Gray');

subplot(2,3,4);
imagesc(B(:,:,1)); colorbar(); colormap('Gray');
subplot(2,3,5);
imagesc(B(:,:,2)); colorbar(); colormap('Gray');
subplot(2,3,6)
imagesc(B(:,:,3)); colorbar(); colormap('Gray');
