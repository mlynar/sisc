%clear all
close all
clc

addpath(genpath('../code'))
addpath(genpath('.'))

%basis params
p.totL = 513;
p.nKernels = 32;
p.kML = 101;
p.exTrs = 0.05; %0.02;

%learning params
p.niter = 20;
p.noiseVar = 1;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.01;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(50000, 1, p.niter);
p.mp.nsnr = repmat(120, 1, p.niter);
p.mp.spkTrs = repmat(0.3, 1, p.niter);
    
%display params
showEvery = 100;
saveEvery = 1;
p.printEvery = 100;

%data params
p.dataPath = '~/Desktop/sounds/training/stereo_sounds.mat';

%run
stereo_cd_sisc
