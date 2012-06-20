%clear all
close all
clc

addpath(genpath('../code'))

%basis params
p.totL = 221; %513;
p.nKernels = 32;
p.kML = 101;
p.exTrs = 0.02;

%learning params
p.niter = 50;
p.noiseVar = 1;
p.nbatch = ceil(linspace(1, 1, p.niter));

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
p.mp.spkTrs = repmat(0.5, 1, p.niter);
    
%display params
showEvery = 100;
saveEvery = 25;
p.printEvery = 1;

%data params
p.dataPath = '~/Desktop/sounds/training/env/';

%test data params
% p.td.ns = 15;
% p.td.spk = 3;
% p.td.kf = 30:32;

%run
cd_sisc
