function smp = sampleToyData(phi, L, spkTrs, kP)
%generates time series of random lenght l \in [minL maxL]
%being sum of kernels phi (in columns) TODO: plus noise
%kP - nTau = floor(kP * L)

if nargin < 4
    kP = 0.01;
end

smp = zeros(1, L);
tauL = floor(kP * L);

kL = size(phi, 1);
hL = floor(kL/2);
rng = -hL:hL;
nKernels = size(phi, 2);

nspikes = 0;

for i = 1:nKernels
    tau = randi(L, tauL, 1);
    for j = 1:tauL
        t = tau(j);
        s = rndlap();
        if abs(s) > spkTrs
            nspikes = nspikes + 1;
            %s = abs(s);
            
            inds = t + rng;
            tInds = inds > 0 & inds < L;
            inds = inds(tInds);
            smp(inds) = smp(inds) + s * phi(tInds,i)';
        end
    end
end
    
smp = smp ./ (max(abs(smp)) + 1e-10);

fprintf('Spike no: %d\n', nspikes);