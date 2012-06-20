function [f, g] = s_objective(s, x, phi, lambda)

L = length(x);
EI = zeros(L,1);
B = size(phi,2);

s = reshape(s, B, L);

hL = floor(size(phi,1) / 2);

for i = 1:size(phi,2)
    EI = EI + conv(s(i,:), phi(:,i), 'same')';
end

E = x' - EI;

f_residue = 0.5 * sum(E.^2);
f_sparse = lambda * sum(abs(s(:)));

f = f_residue + f_sparse;

ds = zeros(size(s));
for t = 1:L
    srt = t - hL; 
    fin = t + hL; 
    
    tind = srt:fin; 
    ind = 1:2*hL+1; ind = ind(tind > 0 & tind <= L);
    tind = tind(tind > 0 & tind <= L);
    
    ds(:,t) = ds(:,t) - (phi(ind,:)' * E(tind));
end

ds = ds + lambda * sign(s);
%ds = ds + lambda * 2 * ds ./ log(1 + ds.^2);

g = ds(:);


if 0
    clf
    figure(3); 
    subplot(2,1,1);
    imagesc(ds, [-0.1 0.1]); colorbar();
    drawnow;
    
    subplot(2,1,2); 
    plot(x); hold on;
    plot(EI, 'r');
    plot(E, 'k');
    legend('signal', 'reconstruction', 'error');
    drawnow
end
    
    