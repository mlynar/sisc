%clear all
close all
clc

addpath(genpath('../code'))
addpath(genpath('.'))

%basis params
p.totL = 513;
p.nKernels = 32; %64; %32;
p.kML = 121; %121
p.exTrs = 0.05; %0.02;

%learning params
p.niter = 50;
p.noiseVar = 1;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.03;
p.map.lambda = 0.05;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(1e6, 1, p.niter);
p.mp.nsnr = repmat(120, 1, p.niter);
p.mp.spkTrs = repmat(0.1, 1, p.niter);
    
%display params
showEvery = 100;
saveEvery = 10;
p.printEvery = 1;

%data params
p.dataPath = '~/Desktop/sounds/training/stereo_sounds_180s.wav';
p.savePath = '~/Desktop/sisc/res/joint_stereo/180s/';
%p.dataPath = '~/Desktop/sounds/training/gmSample1.wav';


%run
%stereo_cd_sisc
%cd_sisc
cd_joint_stereo_sisc
