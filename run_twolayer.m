clear all
close all
clc

addpath(genpath('../code'))
addpath(genpath('.'))

% ===== learning params =====

p.niter = 80;
p.flIter = 50;
p.noiseVar = 1;

% ===== 1st layer params =====

%basis params
p.totL = 513;
p.nKernels = 32;
p.kML = 121;
p.exTrs = 0.05;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.01;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = repmat(1e6, 1, p.niter);
p.mp.nsnr = repmat(120, 1, p.niter);
p.mp.spkTrs = repmat(0.1, 1, p.niter);
p.mp.deadzone = repmat(0, 1, p.niter);


% ===== 2nd layer params =====

%basis params
p.layertwo.totL = 17; 
p.layertwo.nKernels = 32;

%MAP params
p.layertwo.map.mapIter = 100;
p.layertwo.map.lr = 0.01;
p.layertwo.map.lambda = 0.1;
p.layertwo.map.noiseVar = 1;
p.layertwo.map.verbose = false;

%MP params
p.layertwo.mp.nonneg = false;
p.layertwo.mp.mpIter = repmat(floor(0.1 * 16000 * 30), 1, p.niter);
%p.layertwo.mp.nsnr = repmat(120, 1, p.niter);
p.layertwo.mp.spkTrs = repmat(0.1, 1, p.niter);


% ===== display params =====

showEvery = 1;
saveEvery = 5;
p.printEvery = 100;

% ===== data params =====
p.dataPath = '~/Desktop/sounds/training/stereo_sounds_180s.wav';    
p.savePath = '~/Desktop/sisc/res/two_layer/180s/';


% ===== run =====
fprintf('STARTING TWO LAYER TRAINING:\n');
twolayer_cd_sisc
