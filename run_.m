clear all

%basis params
p.totL = 513;
p.nKernels = 5;
p.kML = 201;
p.exTrs = 0.05;

%learning params
p.niter = 500;
p.lr = 0.01;
p.noiseVar = 1;
p.nbatch = ceil(linspace(1, 50, p.niter));

%MAP params
p.map.mapIter = 100;

%MP params
p.mp.nonneg = false;
p.mp.mpIter = ceil(linspace(100, 500, p.niter)); %100
p.mp.nsnr = ceil(linspace(5, 15, p.niter));

%run
sisc
