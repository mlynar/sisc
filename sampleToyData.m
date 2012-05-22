function smp = sampleToyData(phi, minL, maxL, nonneg)
%generates time series of random lenght l \in [minL maxL]
%being sum of kernels phi (in columns) plus noise

if nargin < 4
    nonneg = false;
end
 
l = minL + randi(maxL-minL);
smp = zeros(1,l);

nParts = 100; %randi(floor(l/5));
kL = size(phi, 1);
hL = floor(kL/2);
rng = -hL:hL;
nKernels = size(phi, 2);

for i = 1:nParts
    t = kL + randi(l - 2 * kL);
    kernInd = nKernels - mod(i,nKernels);
    smp(t+rng) = smp(t+rng) + (0.3 + 0.7 * rand()) * phi(:, kernInd)';
end

if nonneg
    smp = smp + abs(min(smp));
end

smp = smp ./ max(abs(smp));
