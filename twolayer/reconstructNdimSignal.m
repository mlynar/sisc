function rec = reconstructNdimSignal(w, phi, k)
% w - nKernels * nTimePoints
% phi - nChannels * nKernelTimePoints * nKernels

if nargin < 3
    k = -1;
end

sL = size(w, 2);
isMultChn = ndims(phi);

if isMultChn == 2
    nK = size(phi, 2);
    L = size(phi, 1);
    nChn = 1;
elseif isMultChn > 2
    nK = size(phi,3);
    L = size(phi, 2);
    nChn = size(phi, 1);
end

hL = floor(L / 2);

rec = zeros(nChn, size(w,2));
for i = 1:nK
    if i == k
        continue;
    end
        
    if isMultChn > 2
        for c = 1:nChn
            si = convnfft(w(i,:), phi(c,:,i));
            rec(c,:) = rec(c,:) + si(hL + (1:sL));
        end
    else       
        si = convnfft(w(i,:), phi(:,i)');
        rec = rec + si(hL + (1:sL));
    end
end