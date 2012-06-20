clear all
close all

addpath(genpath('../code'))

%basis params
p.totL = 513;
p.nKernels = 32;
p.kML = 201;
p.exTrs = 0.05;

%learning params
p.niter = 200;
%p.lr = 
p.noiseVar = 1;
p.nbatch = ceil(linspace(1, 1, p.niter));
p.target_angle = 0.05;

%MAP params
p.map.mapIter = 100;
p.map.lr = 0.01;
p.map.lambda = 0.1;
p.map.noiseVar = 1;
p.map.verbose = false;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = ceil(linspace(10000, 10000, p.niter)); %100
p.mp.nsnr = ceil(linspace(10, 10, p.niter));

%display params
showEvery = p.niter + 1;
saveEvery = p.niter;

%test data params
p.td.ns = 20;
p.td.spk = 3;
p.td.kf = 1:32;

%run
lrT = [5 1 0.5 0.1 0.05 0.01 0.005];

for q = 1:length(lrT)
    p.lr = repmat(lrT(q), 1, p.niter);
    
    sisc
    
    fname = sprintf('./res/lr_%.4f.mat', lrT(q));
    save(fname);
end