function rec = reconstructSignal(w, phi)
% w - nKernels * nTimePoints
% phi - nKernelTimePoints * nKernels

sL = size(w, 2);
L = size(phi, 1);
hL = floor(L / 2);

rec = zeros(1, size(w,2));
for i = 1:size(phi,2)
    si = convnfft(w(i,:), phi(:,i)');
    rec = rec + si(hL + (1:sL));
end