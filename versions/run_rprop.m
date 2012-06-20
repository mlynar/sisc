%clear all
close all
clc

addpath(genpath('../code'))

%basis params
p.totL = 513;
p.nKernels = 16;
p.kML = 101;
p.exTrs = 0.02;

%learning params
p.niter = 200;

d0 = 0.001;
dMax = 1;
dMin = 1e-12;
ksi = 0.1;
etaM = 0.5  ;
etaP = 1.2;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.01;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;


%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(50000, 1, p.niter);
p.mp.nsnr = repmat(800, 1, p.niter);
p.mp.spkTrs = repmat(0.4, 1, p.niter);

%data params
p.dataPath = '~/Desktop/sounds/training/env/';

%test data params
% p.td.ns = 20;
% p.td.spk = 3;
% p.td.kf = 17:32;

%save params
p.saveEvery = 200;

%run
rprop_sisc
