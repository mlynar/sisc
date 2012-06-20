%clear all
close all
clc

addpath(genpath('../code'))

%basis params
p.totL = 513;
p.nKernels = 2;
p.kML = 101;
p.exTrs = 0.05;

%learning params
p.niter = 3000;
p.lr = 0.001; %0.001;
p.noiseVar = 1;
p.nbatch = ceil(linspace(1, 1, p.niter));
p.target_angle = 0.02;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.1;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(50000, 1, p.niter);
p.mp.nsnr = repmat(500, 1, p.niter);
p.mp.spkTrs = repmat(0.3, 1, p.niter);

%display params
showEvery = 100;
saveEvery = 50;

%data params
p.dataPath = '~/Desktop/sounds/training/env/';

%test data params
% p.td.ns = 15;
% p.td.spk = 3;
% p.td.kf = 30:32;

%run
sisc
