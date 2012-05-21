function ws = pimping(sig, phi, coefM, niter)

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

%learning params
%TODO: put as arguments
lambda = 0.05; %0.02;
noiseVar = 1;
lr = 0.03;
trs = 1e-5;

%MSR stopping criterion
oldMSR = 1e12;
newMSR = 0;

ws = coefM';
wsOld = ws;
for n = 1:niter
    fprintf('Iteration %d nnz %d\n', n, nnz(ws));
    
    %compute signal reconstruction
    %TODO: replace with reconstruct_signal function call
    rec = zeros(1, length(sig));
    for i = 1:nB
        si = convnfft(ws(i,:), phi(i,:));
        rec = rec + si(hL + (1:sL));
    end
    
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
        tpoints = find(abs(ws(i,:)) >= trs);
        
        if ~isempty(tpoints)             
            ds = zeros(1,sL);
            for tp = 1:length(tpoints)
                t = tpoints(tp);
                rng = t-hL:t+hL; 
                ind = rng > 0 & rng <= sL;
                rng = rng(ind);
                ds(t) = phi(i, ind) * r(rng)' ./ noiseVar;
            end
            
            %laplace prior 
            %ds = ds + lambda * (-sign(ws(i,:)));
            %"differential" sparse prior
            ds = ds - lambda * 2 * ws(i,:) ./ (1 + ws(i,:).^2);
            ws(i,:) = ws(i,:) + lr * ds;
        end
    end
end

ws = ws';


%for niter
%   compute reconstruction
%   residue = sig - reconstruction
%   for basis_function
%       ds_n = zeros(size(coef(n,:)));
%       tp = find(coef(n,:) ~= 0);
%       for t in tp
%           rng = t-hL:t+hL;
%           ds_n(rng) = basis_function_n .* residue(rng);
%       end
%       ds_n = ds_n + sparse_penalty(ds_n)
%       s_n = s_n - lr * ds_n
%   end
%end