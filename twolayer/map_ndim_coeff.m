function ws = map_ndim_coeff(sig, phi, coefM, niter, lr, lambda, noiseVar, verbose)
%perform MAP inference of activation coefficients based on coefM initialization
%	sig - signal vector
%	phi - kernel dictionary
%	coefM - coefficient matrix, initialized with MP for instance
%	niter - iteration count
%   lr - learning rate
%   lambda - sparseness coefficient
%   noiseVar - noise variance

%learning params
if nargin < 5
    lr = 0.03;
    lambda = 0.05;
    noiseVar = 1;
    verbose = false; %true;
end

stopCrt = 1;

%basis and signal params
sL = length(sig);
%limitation - basis functions of only odd length
L = size(phi,2);
hL = floor(L / 2); 
nB = size(phi,3);
rngK = -hL:hL;
nCh = size(phi, 1);

%MSR stopping criterion
oldMSR = 1e12;
newMSR = 0;

ws = coefM;
wsOld = ws;
for n = 1:niter    
    %compute signal reconstruction
    rec = reconstructNdimSignal(ws, phi);
    
    %compute residue
    r = sig - rec;
    r = [zeros(nCh, L) r zeros(nCh, L)];
    
    if verbose
        fprintf('Iteration %d nnz %d\n', n, nnz(ws));
        fprintf('MSR = %.3f Kurt = %.3f Sum = %.3f\n\n', newMSR, kurtosis(ws(:)), sum(abs(ws(:))));
    end
    
    newMSR = sum(sum(r.^2));
    if newMSR > oldMSR && stopCrt
        ws = wsOld;
        fprintf('MSR increased to %.4f. Stopping\n', newMSR);
        break
    end
            
    oldMSR = newMSR;
    wsOld = ws;
    
    for i = 1:nB
        tpoints = find(ws(i,:) ~= 0);  
        %TODO: this is not optimal - don't keep stuff in a full array
        ds = zeros(1,sL); 

        for tp = 1:length(tpoints)
            t = tpoints(tp);
            rng = t + L + rngK;

            
            ds(t) = sum(sum(phi(:, :, i) .* r(:,rng) ./ noiseVar));
        end

        ds = ds - lambda * sparseLog(ws(i,:));            
        ws(i,:) = ws(i,:) + lr * ds;
    end
    
end

ws = ws';
