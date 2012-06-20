function rec = reconstructSignal(w, phi, k)
% w - nKernels * nTimePoints
% phi - nKernelTimePoints * nKernels

if nargin < 3
    k = -1;
end

sL = size(w, 2);
L = size(phi, 1);
hL = floor(L / 2);

rec = zeros(1, size(w,2));
for i = 1:size(phi,2)
    if i == k
        continue;
    end
    si = convnfft(w(i,:), phi(:,i)');
    rec = rec + si(hL + (1:sL));
end