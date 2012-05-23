function ws = mapCoeff(sig, phi, coefM, niter, lr, lambda, noiseVar)
%perform MAP inference of activation coefficients based on coefM initialization
%	sig - signal vector
%	phi - kernel dictionary
%	coefM - coefficient matrix, initialized with MP for instance
%	niter - iteration count
%   lr - learning rate
%   lambda - sparseness coefficient
%   noiseVar - noise variance

%learning params
if nargin < 4
    lr = 0.03;
    lambda = 0.05;
    noiseVar = 1;
end

phi = phi';
if size(sig, 2) == 1;
    sig = sig';
end

%basis and signal params
sL = length(sig);
%limitation - basis functions of only odd length
L = size(phi,2);
hL = floor(L / 2); 
nB = size(phi,1);

%MSR stopping criterion
oldMSR = 1e12;
newMSR = 0;

ws = coefM';
wsOld = ws;
for n = 1:niter
    fprintf('Iteration %d nnz %d\n', n, nnz(ws));
    
    %compute signal reconstruction
    rec = reconstructSignal(ws, phi');
    
    %compute residue
    r = sig - rec;

    newMSR = sum(r.^2);
    if newMSR > oldMSR
        ws = wsOld;
        fprintf('MSR increased to %.4f. Stopping\n', newMSR);
        break
    end
    
    fprintf('MSR = %.3f Kurt = %.3f Sum = %.3f\n\n', newMSR, kurt(ws(:)), sum(abs(ws(:))));
    
    oldMSR = newMSR;
    wsOld = ws;
    
    for i = 1:nB
        tpoints = find(ws ~= 0);
        
        if ~isempty(tpoints)             
            ds = zeros(1,sL);
            for tp = 1:length(tpoints)
                t = tpoints(tp);
                rng = t-hL:t+hL; 
                ind = rng > 0 & rng <= sL;
                rng = rng(ind);
                ds(t) = phi(i, ind) * r(rng)' ./ noiseVar;
            end
            
            %"differential" sparse prior
            ds = ds - lambda * sparseLog(ws(i,:));
            
            ws(i,:) = ws(i,:) + lr * ds;
        end
    end
end

ws = ws';
